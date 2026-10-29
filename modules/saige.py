# Copyright (C) 2026 Zachary Kenneth Stewart

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os, re, sys, glob, json, shutil, statistics, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from statistics import NormalDist

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _version import __version__
from plink import PLINK, PLINKTool, parse_sample_info
from subprocess_runners import run_rscript
from gwas_reporting import parse_lead_variants, parse_saige_leads, inter_chromosomal_r2, LeadInspector, \
    LeadGeneReporter
from annotarium_importers import import_annotarium_gff3, import_annotarium_annotable
from traits import BinaryTrait, get_trait

SAIGE_SUBMODULE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "external", "SAIGE")
GWAS_MANIFEST_SUFFIX = ".gwas.json"

class SAIGEPhenotype:
    '''
    Builds a SAIGE phenotype file (for step1_fitNULLGLMM.R --phenoFile) from a
    fileset's sample information and a PLINK2 PCA, with the chosen PCs, and
    optionally the ALT allele dosages of variants, as covariates.

    Samples are kept in sample information order, and are excluded if their
    phenotype is missing or they are absent from the PCA. Phenotypes are recoded
    by a trait handler (see traits.py): a binary phenotype in PLINK coding
    (1 = control, 2 = case; 0 or -9 = missing) is recoded to 0/1, while a
    quantitative phenotype (-9 = missing) is kept as it is. The phenotype column
    is named by the trait (e.g., 'y_binary').

    Example:
        pheno = SAIGEPhenotype("prefix.psam", "prefix.pca.eigenvec", BinaryTrait())
        covariates = pheno.write("prefix.pheno", [1, 2])
    '''
    ID_COLUMNS = {"FID", "IID", "SID", "PAT", "MAT", "SEX"}

    def __init__(self, sampleFile, eigenvecFile, trait=None):
        self.sampleFile = sampleFile
        self.eigenvecFile = eigenvecFile
        self.trait = trait if trait is not None else BinaryTrait()
        self.phenoColumn = None
        self.counts = {}

    @property
    def sampleFile(self):
        return self._sampleFile

    @sampleFile.setter
    def sampleFile(self, value):
        '''
        A .psam or .fam file of sample phenotypes.
        '''
        if not isinstance(value, str):
            raise TypeError("sampleFile must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"Sample file '{value}' is not a file")
        self._sampleFile = os.path.abspath(value)
        self.sampleOrder, self.sampleDict = parse_sample_info(self._sampleFile)

    @property
    def eigenvecFile(self):
        return self._eigenvecFile

    @eigenvecFile.setter
    def eigenvecFile(self, value):
        if not isinstance(value, str):
            raise TypeError("eigenvecFile must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"PCA eigenvector file '{value}' is not a file")
        self._eigenvecFile = os.path.abspath(value)
        self.pcaDict, self.numPCs = SAIGEPhenotype.parse_eigenvec(self._eigenvecFile)

    @property
    def phenoColumn(self):
        return self._phenoColumn

    @phenoColumn.setter
    def phenoColumn(self, value):
        '''
        The sample information column holding the phenotype, or None for the
        first non-ID column (PHENO1 in a variantopia .psam).
        '''
        candidates = [ c for c in self.sampleDict[self.sampleOrder[0]] if c not in SAIGEPhenotype.ID_COLUMNS ] \
                     if self.sampleOrder != [] else []
        if value is None:
            if candidates == []:
                raise ValueError(f"No phenotype column found in '{self.sampleFile}'")
            value = candidates[0]
        elif value not in candidates:
            raise ValueError(f"Phenotype column '{value}' not found in '{self.sampleFile}'; options are {candidates}")
        self._phenoColumn = value

    @property
    def phenoName(self):
        '''
        The phenotype's column name in the SAIGE phenotype file.
        '''
        return self.trait.PHENO_NAME

    @staticmethod
    def parse_eigenvec(eigenvecFile):
        '''
        Parameters:
            eigenvecFile -- a string pointing to a PLINK2 .eigenvec file
        Returns:
            pcaDict -- a dictionary with structure like:
                       {
                           iid: [pc1, pc2, ...], # as strings, to preserve precision
                           ...
                       }
            numPCs -- an integer of the number of PCs in the file
        '''
        pcaDict = {}
        with open(eigenvecFile, "r") as fileIn:
            header = fileIn.readline().rstrip("\r\n").lstrip("#").split("\t")
            iidIndex = header.index("IID")
            pcIndices = [ i for _, i in sorted( (int(c[2:]), i) for i, c in enumerate(header)
                                                 if re.fullmatch(r"PC\d+", c) ) ]
            if pcIndices == []:
                raise ValueError(f"No PC columns found in the header of '{eigenvecFile}'")
            for line in fileIn:
                sl = line.rstrip("\r\n").split("\t")
                pcaDict[sl[iidIndex]] = [ sl[i] for i in pcIndices ]
        return pcaDict, len(pcIndices)

    @staticmethod
    def parse_dosages(rawFile):
        '''
        Parses a plink2 --export A .raw file, imputing missing dosages with the
        mean dosage of each variant.

        Parameters:
            rawFile -- a string pointing to a .raw file
        Returns:
            covariates -- a list of covariate names, being 'SNP_' and the variant ID
                          with non-alphanumeric characters replaced by '_'
            dosageDict -- a dictionary with structure like:
                          {
                              iid: [dosage1, dosage2, ...], # as strings
                              ...
                          }
            means -- a list of mean dosages as strings, for samples absent from
                     the .raw file
        '''
        with open(rawFile, "r") as fileIn:
            header = fileIn.readline().rstrip("\r\n").split("\t")
            covariates = [ "SNP_" + re.sub(r"[^A-Za-z0-9]", "_", c.rsplit("_", 1)[0]) for c in header[6:] ]
            rows = { sl[1]: sl[6:] for sl in ( line.rstrip("\r\n").split("\t") for line in fileIn ) }

        means = []
        for i in range(len(covariates)):
            called = [ float(d[i]) for d in rows.values() if d[i] != "NA" ]
            if called == []:
                raise ValueError(f"Covariate variant '{header[6+i]}' has no called genotypes")
            means.append(f"{sum(called) / len(called):.6g}")
        dosageDict = { iid: [ means[i] if v == "NA" else v for i, v in enumerate(d) ]
                       for iid, d in rows.items() }
        return covariates, dosageDict, means

    def write(self, outputFileName, pcList, rawFile=None):
        '''
        Writes the SAIGE phenotype file. Sample counts, and a summary of the
        phenotypes written (see the trait's summarise()), are recorded in the
        counts property.

        Parameters:
            outputFileName -- a string of the location to write the file
            pcList -- a list of integers of the PCs to include as covariates
            rawFile -- (OPTIONAL) a string pointing to a plink2 --export A .raw
                       file of variant dosages to include as covariates
        Returns:
            covariates -- a list of the covariate column names
        '''
        if pcList != [] and max(pcList) > self.numPCs:
            raise ValueError(f"PC{max(pcList)} requested but only {self.numPCs} PCs are in '{self.eigenvecFile}'")
        pcNames = [ f"PC{i}" for i in pcList ]
        snpNames, dosageDict, means = SAIGEPhenotype.parse_dosages(rawFile) if rawFile is not None else ([], {}, [])

        self.counts = { "written": 0, "missingPheno": 0, "missingPCA": 0 }
        phenotypes = []
        with open(outputFileName, "w") as fileOut:
            fileOut.write("\t".join([self.phenoName] + pcNames + ["IID"] + snpNames) + "\n")
            for iid in self.sampleOrder:
                pheno = self.trait.recode(self.sampleDict[iid][self.phenoColumn])
                if pheno is None:
                    self.counts["missingPheno"] += 1
                    continue
                if iid not in self.pcaDict:
                    self.counts["missingPCA"] += 1
                    continue
                row = [pheno] + [ self.pcaDict[iid][i-1] for i in pcList ] + [iid] + dosageDict.get(iid, means)
                fileOut.write("\t".join(row) + "\n")
                self.counts["written"] += 1
                phenotypes.append(pheno)
        if self.counts["written"] == 0:
            raise ValueError("No samples have both a phenotype and PCs; check that IIDs match")
        self.counts["phenotype"] = self.trait.summarise(phenotypes)
        return pcNames + snpNames

class SAIGE(PLINKTool):
    '''
    Runs a SAIGE single-variant association test of a binary or quantitative
    trait, using the R scripts of a SAIGE source checkout (by default,
    variantopia's external/SAIGE submodule) with an Rscript whose R library has
    the SAIGE package installed:
        1. A sparse GRM from LD-pruned variants
        2. The null mixed model (logistic for a binary trait, linear for a
           quantitative one), fitted with the dense GRM of the pruned variants,
           leaving one chromosome out (LOCO)
        3. Per-chromosome tests, run in parallel and merged; binary traits are
           tested with SPA, with Firth correction for P < 0.05
    The trait's SAIGE settings come from its trait handler (see traits.py).

    SAIGE reads its genotypes as a bed fileset, so pgen and bgen inputs are
    converted within a temporary working directory. plink2 settings are
    inherited from PLINKTool. SAIGE fixes its random seed, but its null model
    fit is not bit-for-bit reproducible (e.g., tau can differ in the 7th
    significant figure between runs), so re-runs give equivalent rather than
    identical results.

    Example:
        saige = SAIGE()
        saige.rscriptExe = "/path/to/SAIGE/env/bin/Rscript"
        with saige.work_dir("output/prefix") as workDir:
            grmBed = saige.as_bed(PLINK("prefix.pruned", "pgen"), workDir)
            sparseGRM = saige.sparse_grm(grmBed, "output/prefix")
            saige.fit_null(grmBed, sparseGRM, "output/prefix.pheno", ["PC1", "PC2"], "output/prefix")
            saige.test(grmBed, sparseGRM, "output/prefix", "output/prefix.saige.txt", workDir)
    '''
    SCRIPTS = { "grm": "createSparseGRM.R", "null": "step1_fitNULLGLMM.R", "test": "step2_SPAtests.R" }

    def __init__(self):
        super().__init__()
        self.outputFormat = "bed"
        self.saigeDir = SAIGE_SUBMODULE_DIR
        self._rscriptExe = shutil.which("Rscript") # may be None; validated when SAIGE is first run
        self.jobs = 1
        self.trait = BinaryTrait()
        self.relCutoff = 0.125
        self.grmMarkers = 20000
        self.minMAF = 0.01
        self.minMAC = 20

    @property
    def saigeDir(self):
        return self._saigeDir

    @saigeDir.setter
    def saigeDir(self, value):
        '''
        A SAIGE source checkout, from which the extdata R scripts are run.
        '''
        if not isinstance(value, str):
            raise TypeError("saigeDir must be a string")
        value = os.path.abspath(value)
        missing = [ s for s in SAIGE.SCRIPTS.values() if not os.path.isfile(os.path.join(value, "extdata", s)) ]
        if missing != []:
            raise FileNotFoundError(f"SAIGE directory '{value}' is missing extdata script(s) {missing}")
        self._saigeDir = value

    def source_version(self):
        '''
        Returns the SAIGE version of saigeDir's DESCRIPTION file, and its git
        commit (or None if it is not a git checkout).
        '''
        with open(os.path.join(self.saigeDir, "DESCRIPTION"), "r") as fileIn:
            version = next( (line.split(":", 1)[1].strip() for line in fileIn if line.startswith("Version:")), None )
        result = subprocess.run(["git", "-C", self.saigeDir, "rev-parse", "HEAD"],
                                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        commit = result.stdout.decode("utf-8").strip() if result.returncode == 0 else None
        return version, commit

    @property
    def rscriptExe(self):
        return self._rscriptExe

    @rscriptExe.setter
    def rscriptExe(self, value):
        '''
        An Rscript executable whose R library has the SAIGE package installed.
        '''
        if not isinstance(value, str):
            raise TypeError("rscriptExe must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"Rscript executable not found at '{value}'")
        self._rscriptExe = os.path.abspath(value)

    @property
    def jobs(self):
        return self._jobs

    @jobs.setter
    def jobs(self, value):
        '''
        Number of chromosomes to test in parallel; each job needs its own memory
        (~2.5 GB for ~250 samples).
        '''
        if not isinstance(value, int) or value < 1:
            raise ValueError("jobs must be a positive integer")
        self._jobs = value

    @property
    def relCutoff(self):
        return self._relCutoff

    @relCutoff.setter
    def relCutoff(self, value):
        '''
        Relatedness below which sparse GRM entries are set to 0.
        '''
        if not isinstance(value, (int, float)) or not 0 <= value <= 1:
            raise ValueError("relCutoff must be a number from 0 to 1")
        self._relCutoff = value

    @property
    def grmMarkers(self):
        return self._grmMarkers

    @grmMarkers.setter
    def grmMarkers(self, value):
        '''
        Number of randomly selected markers for the sparse GRM, capped at the
        number of variants in the GRM fileset.
        '''
        if not isinstance(value, int) or value < 1:
            raise ValueError("grmMarkers must be a positive integer")
        self._grmMarkers = value

    @property
    def minMAF(self):
        return self._minMAF

    @minMAF.setter
    def minMAF(self, value):
        if not isinstance(value, (int, float)) or not 0 <= value <= 0.5:
            raise ValueError("minMAF must be a number from 0 to 0.5")
        self._minMAF = value

    @property
    def minMAC(self):
        return self._minMAC

    @minMAC.setter
    def minMAC(self, value):
        '''
        Minimum minor allele count to test; the single variance ratio is
        estimated from markers with MAC >= 20, matching the default.
        '''
        if not isinstance(value, (int, float)) or value < 0:
            raise ValueError("minMAC must be a number >= 0")
        self._minMAC = value

    def plink(self, prefix, fileFormat):
        '''
        Returns a PLINK object for a fileset, with this object's plink2 settings.
        '''
        plink = PLINK(prefix, fileFormat)
        plink.threads, plink.memory = self.threads, self.memory
        if self.plinkExe is not None:
            plink.plinkExe = self.plinkExe
        return plink

    def as_bed(self, plink, workDir):
        '''
        Returns the prefix of a fileset in bed format, converting it into workDir
        if necessary.

        Parameters:
            plink -- a PLINK object of the fileset
            workDir -- a string pointing to the working directory
        Returns:
            bedPrefix -- a string of the bed fileset's prefix
        '''
        if plink.fileFormat == "bed":
            return plink.prefix
        bedPrefix = os.path.join(workDir, os.path.basename(plink.prefix))
        self.write_fileset(plink.inputArgs, bedPrefix, workDir)
        return bedPrefix

    def export_dosages(self, plink, variantIDs, workDir):
        '''
        Exports the ALT allele dosages of variants with plink2 --export A.

        Parameters:
            plink -- a PLINK object of the fileset holding the variants
            variantIDs -- a list of strings of variant IDs
            workDir -- a string pointing to the working directory
        Returns:
            rawFile -- a string of the .raw file's location
        '''
        idFile = os.path.join(workDir, "covar_snps.txt")
        with open(idFile, "w") as fileOut:
            fileOut.write("\n".join(variantIDs) + "\n")

        # Find each variant's ALT allele, since --export A counts REF by default
        dosagePrefix = os.path.join(workDir, "covar_snps")
        plink.run_plink(plink.inputArgs + ["--extract", idFile, "--make-just-pvar", "cols="], dosagePrefix)
        alleleFile = dosagePrefix + ".alleles.txt"
        found = []
        with open(dosagePrefix + ".pvar", "r") as fileIn, open(alleleFile, "w") as fileOut:
            for line in fileIn:
                if line.startswith("#"):
                    continue
                chrom, pos, varId, ref, alt = line.rstrip("\r\n").split("\t")[:5]
                if "," in alt:
                    raise ValueError(f"Covariate variant '{varId}' is multiallelic")
                fileOut.write(f"{varId}\t{alt}\n")
                found.append(varId)
        missing = [ v for v in variantIDs if v not in found ]
        if missing != []:
            raise ValueError(f"Covariate variant(s) {missing} were not found in '{plink.prefix}'")

        plink.run_plink(plink.inputArgs + ["--extract", idFile, "--export", "A",
                                           "--export-allele", alleleFile], dosagePrefix)
        return dosagePrefix + ".raw"

    def run_saige(self, step, args, logFile):
        '''
        Runs a SAIGE extdata R script, writing its stdout and stderr to logFile.

        Parameters:
            step -- a string of the SCRIPTS key of the script to run
            args -- a list of strings of '--option=value' arguments
            logFile -- a string of the location to write the log
        '''
        if self.rscriptExe is None:
            raise FileNotFoundError("Rscript is not discoverable in your system PATH and rscriptExe was not set")
        cmd = [self.rscriptExe, os.path.join(self.saigeDir, "extdata", SAIGE.SCRIPTS[step])] + args
        with open(logFile, "w") as logOut:
            result = subprocess.run(cmd, stdout=logOut, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            with open(logFile, "r") as logIn:
                logTail = "".join(logIn.readlines()[-15:])
            raise Exception(f"SAIGE {SAIGE.SCRIPTS[step]} failed with exit code {result.returncode}; command " +
                            f"was '{' '.join(cmd)}'. Log tail ({logFile}):\n{logTail}")

    def sparse_grm(self, grmBedPrefix, outputPrefix, workDir):
        '''
        Creates a sparse GRM from a bed fileset of LD-pruned variants.

        Parameters:
            grmBedPrefix -- a string of the bed fileset's prefix
            outputPrefix -- a string of the output prefix; files are written to
                            '<outputPrefix>.sparseGRM.{mtx,mtx.sampleIDs.txt,log}'
            workDir -- a string pointing to the working directory
        Returns:
            sparseGRM -- a string of the .sparseGRM.mtx location
        '''
        sparseGRM = outputPrefix + ".sparseGRM.mtx"
        numMarkers = min(self.grmMarkers, PLINK.count_lines(grmBedPrefix + ".bim"))

        # SAIGE names its outputs by its settings, so write them to a fresh directory
        grmDir = os.path.join(workDir, "sparse_grm")
        os.makedirs(grmDir)
        self.run_saige("grm", [f"--plinkFile={grmBedPrefix}", f"--nThreads={self.threads}",
                               f"--outputPrefix={os.path.join(grmDir, 'grm')}",
                               f"--numRandomMarkerforSparseKin={numMarkers}",
                               f"--relatednessCutoff={self.relCutoff}"],
                       outputPrefix + ".sparseGRM.log")
        mtxFiles = glob.glob(os.path.join(grmDir, "*.sparseGRM.mtx"))
        if len(mtxFiles) != 1:
            raise Exception(f"Expected one sparse GRM from SAIGE in '{grmDir}', but found {len(mtxFiles)}")
        shutil.move(mtxFiles[0], sparseGRM)
        shutil.move(mtxFiles[0] + ".sampleIDs.txt", sparseGRM + ".sampleIDs.txt")
        return sparseGRM

    def fit_null(self, grmBedPrefix, sparseGRM, phenoFile, covariates, outputPrefix):
        '''
        Fits the null mixed model with LOCO. A single variance ratio is
        estimated from pruned markers with MAC >= 20; categorical variance ratios
        would force LOCO off in SAIGE, and are intended for rare-variant tests.

        Parameters:
            grmBedPrefix -- a string of the LD-pruned bed fileset's prefix
            sparseGRM -- a string pointing to the .sparseGRM.mtx file
            phenoFile -- a string pointing to the SAIGE phenotype file
            covariates -- a list of strings of covariate column names
            outputPrefix -- a string of the output prefix; files are written to
                            '<outputPrefix>.null.{rda,varianceRatio.txt,log}'
        Returns:
            tau -- a float of the GRM variance component estimate
        '''
        args = [f"--plinkFile={grmBedPrefix}", f"--sparseGRMFile={sparseGRM}",
                f"--sparseGRMSampleIDFile={sparseGRM}.sampleIDs.txt", f"--phenoFile={phenoFile}",
                f"--phenoCol={self.trait.PHENO_NAME}", "--sampleIDColinphenoFile=IID"] + self.trait.null_args() + \
               ["--isCateVarianceRatio=FALSE", "--useSparseGRMtoFitNULL=FALSE",
                "--useSparseGRMforVarRatio=TRUE", "--LOCO=TRUE", f"--nThreads={self.threads}",
                f"--outputPrefix={outputPrefix}.null", "--IsOverwriteVarianceRatioFile=TRUE"]
        if covariates != []:
            args.append(f"--covarColList={','.join(covariates)}")
        self.run_saige("null", args, outputPrefix + ".null.log")
        return SAIGE.parse_tau(outputPrefix + ".null.log")

    @staticmethod
    def parse_tau(nullLog):
        '''
        Parses the final (converged) GRM variance component, tau, from a SAIGE
        step 1 log (see parse_variance_components()).
        '''
        return SAIGE.parse_variance_components(nullLog)[1]

    @staticmethod
    def parse_variance_components(nullLog):
        '''
        Parses the final (converged) variance components from a SAIGE step 1 log;
        SAIGE logs them as 'Final <tau1> <tau2> :', where tau1 is the dispersion
        (fixed at 1 for a binary trait) or the residual variance (for a
        quantitative trait), and tau2 is the GRM variance component.

        Returns:
            tau1, tau2 -- floats of the variance components
        '''
        components = None
        with open(nullLog, "r") as fileIn:
            for line in fileIn:
                if line.startswith("Final "):
                    components = [float(line.split()[1]), float(line.split()[2])]
        if components is None:
            raise ValueError(f"No final tau estimate found in SAIGE null model log '{nullLog}'")
        return components

    def test(self, testBedPrefix, sparseGRM, nullPrefix, outputFileName, workDir):
        '''
        Tests each chromosome of a bed fileset in parallel (jobs at a time) and
        merges the results in fileset chromosome order.

        Parameters:
            testBedPrefix -- a string of the prefix of the bed fileset to test
            sparseGRM -- a string pointing to the .sparseGRM.mtx file
            nullPrefix -- a string of the null model's prefix, as given to fit_null()
            outputFileName -- a string of the location to write merged results
            workDir -- a string pointing to the working directory; a failed
                       chromosome's log is copied alongside outputFileName
        Returns:
            numTested -- an integer of the number of variants tested
        '''
        chroms = []
        with open(testBedPrefix + ".bim", "r") as fileIn:
            for line in fileIn:
                chrom = line.split("\t", 1)[0]
                if chroms == [] or chroms[-1] != chrom:
                    if chrom in chroms:
                        raise ValueError(f"Variants of chromosome '{chrom}' in '{testBedPrefix}.bim' are not contiguous")
                    chroms.append(chrom)

        def test_chrom(chrom):
            chromPrefix = os.path.join(workDir, f"saige.chr{chrom}")
            try:
                self.run_saige("test", [f"--bedFile={testBedPrefix}.bed", f"--bimFile={testBedPrefix}.bim",
                                        f"--famFile={testBedPrefix}.fam", "--AlleleOrder=alt-first",
                                        f"--SAIGEOutputFile={chromPrefix}.txt",
                                        f"--GMMATmodelFile={nullPrefix}.null.rda",
                                        f"--varianceRatioFile={nullPrefix}.null.varianceRatio.txt",
                                        f"--sparseGRMFile={sparseGRM}",
                                        f"--sparseGRMSampleIDFile={sparseGRM}.sampleIDs.txt",
                                        f"--chrom={chrom}", "--LOCO=TRUE", f"--minMAF={self.minMAF}",
                                        f"--minMAC={self.minMAC}"] + self.trait.test_args() +
                                       ["--is_output_moreDetails=TRUE"],
                               chromPrefix + ".log")
            except Exception:
                shutil.copy(chromPrefix + ".log", f"{os.path.splitext(outputFileName)[0]}.chr{chrom}.log")
                raise
            return chromPrefix + ".txt"

        with ThreadPoolExecutor(max_workers=self.jobs) as executor:
            chromFiles = list(executor.map(test_chrom, chroms))

        numTested = 0
        with open(outputFileName, "w") as fileOut:
            for i, chromFile in enumerate(chromFiles):
                with open(chromFile, "r") as fileIn:
                    header = fileIn.readline()
                    if i == 0:
                        fileOut.write(header)
                    for line in fileIn:
                        fileOut.write(line)
                        numTested += 1
        return numTested

    @staticmethod
    def lambda_gc(resultFile):
        '''
        Computes the genomic inflation factor from a SAIGE results file, as the
        median 1-df chi-squared statistic of the P values over its expectation.

        Parameters:
            resultFile -- a string pointing to a SAIGE step 2 results file
        Returns:
            lambdaGC -- a float of the genomic inflation factor
            numTested -- an integer of the number of non-missing P values
            numSuggestive -- an integer of the number of P values < 1e-5
        '''
        normal = NormalDist()
        chisq = []
        numSuggestive = 0
        with open(resultFile, "r") as fileIn:
            pIndex = fileIn.readline().rstrip("\r\n").split("\t").index("p.value")
            for line in fileIn:
                p = line.rstrip("\r\n").split("\t")[pIndex]
                if p == "NA":
                    continue
                p = float(p)
                chisq.append(normal.inv_cdf(p / 2) ** 2 if p > 0 else float("inf"))
                numSuggestive += p < 1e-5
        if chisq == []:
            raise ValueError(f"No P values found in '{resultFile}'")
        return statistics.median(chisq) / normal.inv_cdf(0.25) ** 2, len(chisq), numSuggestive

def parse_pc_set(value, asList=False):
    '''
    Parses a PC set as given on the command line: an integer N for the first N
    PCs (0 for none), or a comma-separated list of PCs (e.g., '1,2,4', or '3,'
    for PC3 alone).

    Parameters:
        value -- a string or integer of the PC set
        asList -- a boolean indicating that value is always a list of PCs, so
                  that '3' means PC3 alone
    Returns:
        pcList -- a list of integers of PCs
        tag -- a string naming the set (e.g., 'pc2' or 'pc1-2-4')
    '''
    value = str(value).strip()
    if asList or "," in value:
        try:
            pcList = [ int(x) for x in value.split(",") if x != "" ]
        except ValueError:
            raise ValueError(f"PC set '{value}' is not a comma-separated list of integers")
        if pcList == [] or any(x < 1 for x in pcList) or len(set(pcList)) != len(pcList):
            raise ValueError(f"PC set '{value}' must list unique PCs >= 1")
        return pcList, "pc" + "-".join(map(str, pcList))
    try:
        numPCs = int(value)
    except ValueError:
        raise ValueError(f"PC set '{value}' is not an integer or a comma-separated list of integers")
    if numPCs < 0:
        raise ValueError(f"PC set '{value}' must be >= 0")
    return list(range(1, numPCs + 1)), f"pc{numPCs}"

def parse_calibration(calibrationFile):
    '''
    Finds the suggested PC set in a 'saige calibrate' table.

    Returns:
        pcList -- a list of integers of the suggested PCs
        tag -- a string of the suggested set's tag
    '''
    with open(calibrationFile, "r") as fileIn:
        header = fileIn.readline().rstrip("\r\n").split("\t")
        for line in fileIn:
            row = dict(zip(header, line.rstrip("\r\n").split("\t")))
            if row.get("SUGGESTED") == "yes":
                if row["COVARIATES"] == "none":
                    return [], row["TAG"]
                return [ int(c[2:]) for c in row["COVARIATES"].split(",") ], row["TAG"]
    raise ValueError(f"No suggested PC set found in calibration table '{calibrationFile}'")

def configure_saige(args):
    '''
    Creates a SAIGE object from the arguments shared by "saige" mode commands.
    '''
    saige = SAIGE()
    saige.plinkExe = args.plink2
    saige.threads = args.threads
    saige.memory = args.memory
    saige.tmpDir = args.tmpDir
    saige.keepTmp = args.keepTmp
    saige.saigeDir = args.saigeDir
    saige.rscriptExe = args.saigeRscript
    saige.jobs = args.jobs
    saige.relCutoff = args.relCutoff
    saige.grmMarkers = args.grmMarkers
    saige.minMAF = args.minMAF
    saige.minMAC = args.minMAC
    saige.trait = get_trait(args.traitType, **({ "invNormalize": not args.noInvNormalize }
                                               if args.traitType == "quantitative" else {}))
    return saige

def saige_calibrate(args):
    '''
    Handles "saige calibrate" mode of variantopia.
    '''
    saige = configure_saige(args)
    plink = saige.plink(args.plinkPrefix, args.inputFormat)
    outputPrefix = os.path.join(args.outputDirectory, os.path.basename(args.plinkPrefix))
    tableFile = outputPrefix + ".pc_calibration.tsv"
    plotFile = outputPrefix + ".pc_calibration.png"
    pcSets = [ parse_pc_set(s) for s in args.pcSets ]
    saige.check_outputs_absent([tableFile, plotFile, outputPrefix + ".sparseGRM.mtx"] +
                               [ f"{outputPrefix}.{tag}.null.log" for _, tag in pcSets ])

    # Fit the null model and test the pruned variants with each PC set; the pruned
    # variants are approximately independent and mostly null, so lambda should
    # approach 1 once confounding structure is adequately modelled
    header = ["TAG", "COVARIATES", "TAU", "LAMBDA_GC", "N_TESTED", "N_P_LT_1E-5"]
    rows = []
    with saige.work_dir(outputPrefix) as workDir:
        bedPrefix = saige.as_bed(plink, workDir)
        sparseGRM = saige.sparse_grm(bedPrefix, outputPrefix, workDir)
        pheno = SAIGEPhenotype(plink.write_psam(os.path.join(workDir, "samples")), args.eigenvecFile, saige.trait)
        pheno.phenoColumn = args.phenoColumn

        print(f"# Trait type: {saige.trait.TYPE}")
        print("# " + "\t".join(header))
        for pcList, tag in pcSets:
            runPrefix = f"{outputPrefix}.{tag}"
            phenoFile = os.path.join(workDir, f"{tag}.pheno")
            covariates = pheno.write(phenoFile, pcList)
            tau = saige.fit_null(bedPrefix, sparseGRM, phenoFile, covariates, os.path.join(workDir, tag))
            shutil.move(os.path.join(workDir, f"{tag}.null.log"), runPrefix + ".null.log")
            resultFile = os.path.join(workDir, f"{tag}.saige.txt")
            saige.test(bedPrefix, sparseGRM, os.path.join(workDir, tag), resultFile, workDir)
            lambdaGC, numTested, numSuggestive = SAIGE.lambda_gc(resultFile)
            rows.append([tag, ",".join(covariates) or "none", f"{tau:.7g}", f"{lambdaGC:.4f}",
                         str(numTested), str(numSuggestive)])
            print("# " + "\t".join(rows[-1]))
            for f in glob.glob(os.path.join(workDir, f"{tag}.*")) + glob.glob(os.path.join(workDir, "saige.chr*")):
                os.remove(f)

    # Suggest the set with lambda closest to 1, preferring fewer covariates on ties
    suggested = min(range(len(rows)), key=lambda i: (abs(float(rows[i][3]) - 1), len(pcSets[i][0])))
    with open(tableFile, "w") as fileOut:
        fileOut.write("\t".join(header + ["SUGGESTED"]) + "\n")
        for i, row in enumerate(rows):
            fileOut.write("\t".join(row + ["yes" if i == suggested else "no"]) + "\n")
    print(f"# Suggested PC set (lambda closest to 1): {rows[suggested][0]} ({rows[suggested][1]})")
    print(f"# Wrote {tableFile}")

    run_rscript("plot_calibration.R", ["-i", tableFile, "-o", outputPrefix + ".pc_calibration",
                                       "--title", args.title],
                args.rscript, logFile=outputPrefix + ".pc_calibration.plot.log")
    print(f"# Wrote {plotFile}")

def saige_gwas(args):
    '''
    Handles "saige gwas" mode of variantopia.
    '''
    saige = configure_saige(args)
    plink = saige.plink(args.plinkPrefix, args.inputFormat)
    grmPlink = saige.plink(args.grmPrefix, args.grmFormat)
    outputPrefix = os.path.join(args.outputDirectory, os.path.basename(args.plinkPrefix))
    resultFile = outputPrefix + ".saige.txt"
    saige.check_outputs_absent([ outputPrefix + suffix for suffix in
                                 [".pheno", ".sparseGRM.mtx", ".null.rda", ".null.varianceRatio.txt", ".saige.txt"] ])

    # Choose PC covariates
    if args.calibrationFile is not None:
        pcList, tag = parse_calibration(args.calibrationFile)
        print(f"# Using the suggested PC set from '{args.calibrationFile}': {tag}")
    elif args.pcs is not None:
        pcList, tag = parse_pc_set(args.pcs, asList=True)
    else:
        pcList, tag = parse_pc_set(args.npcs)

    with saige.work_dir(outputPrefix) as workDir:
        # Build the phenotype file, optionally with lead variant dosages as covariates
        pheno = SAIGEPhenotype(plink.write_psam(os.path.join(workDir, "samples")), args.eigenvecFile, saige.trait)
        pheno.phenoColumn = args.phenoColumn
        rawFile = saige.export_dosages(plink, args.covarSNPs, workDir) if args.covarSNPs is not None else None
        covariates = pheno.write(outputPrefix + ".pheno", pcList, rawFile)
        counts = pheno.counts
        print(f"# Samples: {counts['written']} ({saige.trait.TYPE} trait; " +
              f"{saige.trait.describe(counts['phenotype'])})")
        print(f"# Excluded {counts['missingPheno']} samples with a missing phenotype, and " +
              f"{counts['missingPCA']} phenotyped samples absent from the PCA")
        print(f"# Covariates: {','.join(covariates) or 'none'}")

        # Fit the null model with the pruned variants, then test all variants
        grmBed = saige.as_bed(grmPlink, workDir)
        grmVariants = PLINK.count_lines(grmBed + ".bim")
        sparseGRM = saige.sparse_grm(grmBed, outputPrefix, workDir)
        tau = saige.fit_null(grmBed, sparseGRM, outputPrefix + ".pheno", covariates, outputPrefix)
        tauResidual = SAIGE.parse_variance_components(outputPrefix + ".null.log")[0]
        print(f"# Null model GRM variance (tau): {tau:.7g}")
        if saige.trait.TYPE == "quantitative":
            print(f"# Null model residual variance: {tauResidual:.7g}")
        testBed = saige.as_bed(plink, workDir)
        numTested = saige.test(testBed, sparseGRM, outputPrefix, resultFile, workDir)
    print(f"# Variants tested: {numTested}")
    print(f"# Wrote {resultFile}")

    # Record the inputs and settings for 'saige summarise'
    manifestFile = write_gwas_manifest(outputPrefix, saige, plink, grmPlink, {
        "eigenvec": args.eigenvecFile, "phenoColumn": pheno.phenoColumn, "pcs": pcList, "pcTag": tag,
        "covarSNPs": args.covarSNPs, "samples": counts["written"], "phenotype": counts["phenotype"],
        "tau": tau, "tauResidual": tauResidual, "variantsTested": numTested, "grmVariants": grmVariants })
    print(f"# Wrote {manifestFile}")

def file_record(fileName, relativeTo):
    '''
    Describes a file by its location relative to a directory, and its size, so
    that it can be found after the directories are moved together and changes
    to it can be noticed.
    '''
    return { "path": os.path.relpath(fileName, relativeTo), "size": os.path.getsize(fileName) }

def check_file_record(record, relativeTo, description):
    '''
    Returns the absolute location of a file described by file_record(), raising
    an error if it is missing or its size has changed.
    '''
    fileName = os.path.normpath(os.path.join(relativeTo, record["path"]))
    if not os.path.isfile(fileName):
        raise FileNotFoundError(f"{description} '{fileName}' no longer exists")
    if os.path.getsize(fileName) != record["size"]:
        raise ValueError(f"{description} '{fileName}' has changed since it was recorded")
    return fileName

def write_gwas_manifest(outputPrefix, saige, plink, grmPlink, analysis):
    '''
    Writes a record of a 'saige gwas' run's inputs and settings to
    '<outputPrefix>.gwas.json', for use by 'saige summarise'. Locations are
    recorded relative to the output directory. The trait type and its settings
    are taken from the SAIGE object; manifests without a trait type are of a
    binary trait.

    Parameters:
        outputPrefix -- a string of the run's output prefix
        saige -- the SAIGE object used for the run
        plink -- a PLINK object of the tested fileset
        grmPlink -- a PLINK object of the GRM fileset
        analysis -- a dictionary of run details, structured like:
                    {
                        "eigenvec": eigenvecFile, # string
                        "phenoColumn": column, # string
                        "pcs": [1, 2], # list of integers
                        "pcTag": tag, # string
                        "covarSNPs": [variantID, ...], # list of strings, or None
                        "samples": int, # samples in the null model
                        "phenotype": summary, # dictionary from the trait's summarise()
                        "tau": float, "tauResidual": float, # variance components
                        "variantsTested": int,
                        "grmVariants": int # variants in the GRM fileset
                    }
    Returns:
        manifestFile -- a string of the manifest's location
    '''
    manifestFile = outputPrefix + GWAS_MANIFEST_SUFFIX
    outputDir = os.path.dirname(manifestFile)
    saigeVersion, saigeCommit = saige.source_version()
    manifest = {
        "variantopia": __version__,
        "saige": { "version": saigeVersion, "commit": saigeCommit },
        "basename": os.path.basename(outputPrefix),
        "results": file_record(outputPrefix + ".saige.txt", outputDir),
        "pheno": file_record(outputPrefix + ".pheno", outputDir),
        "test": { "prefix": os.path.relpath(plink.prefix, outputDir), "format": plink.fileFormat,
                  "files": [ file_record(f, outputDir) for f in plink.files ] },
        "grm": { "prefix": os.path.relpath(grmPlink.prefix, outputDir), "format": grmPlink.fileFormat,
                 "files": [ file_record(f, outputDir) for f in grmPlink.files ],
                 "variants": analysis["grmVariants"] },
        "eigenvec": file_record(analysis["eigenvec"], outputDir),
        "traitType": saige.trait.TYPE,
        "traitSettings": saige.trait.settings,
        "phenoColumn": analysis["phenoColumn"],
        "pcs": analysis["pcs"],
        "pcTag": analysis["pcTag"],
        "covarSNPs": analysis["covarSNPs"],
        "settings": { "relCutoff": saige.relCutoff, "grmMarkers": saige.grmMarkers,
                      "minMAF": saige.minMAF, "minMAC": saige.minMAC },
        "samples": analysis["samples"],
        "phenotype": analysis["phenotype"],
        "tau": analysis["tau"],
        "tauResidual": analysis["tauResidual"],
        "variantsTested": analysis["variantsTested"]
    }
    with open(manifestFile, "w") as fileOut:
        json.dump(manifest, fileOut, indent=4)
        fileOut.write("\n")
    return manifestFile

def find_gwas_manifest(gwasDir):
    '''
    Returns the location of the single 'saige gwas' manifest within a directory.
    '''
    manifestFiles = glob.glob(os.path.join(gwasDir, "*" + GWAS_MANIFEST_SUFFIX))
    if len(manifestFiles) != 1:
        raise FileNotFoundError(f"Expected one '*{GWAS_MANIFEST_SUFFIX}' file from 'saige gwas' in '{gwasDir}', " +
                                f"but found {len(manifestFiles)}")
    return manifestFiles[0]

class SAIGESummary:
    '''
    Summarises a 'saige gwas' run, writing into its 'summarise' subdirectory:
        1. plot -- Manhattan and QQ plots, top hits, per-chromosome lambda GC and
                   associated regions (plot_manhattan.R), with Bonferroni and
                   Meff thresholds; Meff is the number of GRM (LD-pruned) variants
        2. clump -- LD clumping of associated variants into independent loci
                    (plink2 --clump), with LD from the samples in the null model
        3. inspect -- the direction of each lead variant's effect, the families its
                      effect allele segregates in, its parental origin and
                      Mendelian errors (where parents are genotyped), local support,
                      and LD with stronger lead variants on other chromosomes,
                      judged against the LD between unlinked (LD-pruned) variants
                      on different chromosomes
        4. annotate -- (OPTIONAL) genes near the lead variants of regions and clumps

    The tested fileset and samples are taken from the run's manifest. Each step's
    settings (and the variantopia version) are recorded in '<basename>.summary.json',
    so that a re-run skips steps whose settings, version and outputs are unchanged,
    and re-uses recorded settings that are not given again (e.g., to add an
    annotation to an earlier summary). The trait type is taken from the manifest.

    Example:
        summary = SAIGESummary("output/saige/pc2")
        options = summary.resolve_options({ "suggestive": None, ... })
    '''
    DEFAULTS = { "title": None, "suggestive": 1e-5, "regionGap": 2000000, "topN": 50,
                 "clumpP2": 0.01, "clumpR2": 0.2, "clumpKb": 2000, "supportKb": 250, "linkR2": "auto",
                 "gff3File": None, "annotableFile": None, "chromMap": None,
                 "radius": 50000, "annotateP": None, "annotariumDir": None }
    ANNOTATION_FILES = ["gff3File", "annotableFile", "chromMap"]
    ANNOTATION_OPTIONS = ANNOTATION_FILES + ["radius", "annotateP"]
    STEP_OPTIONS = {
        "plot": ["title", "suggestive", "regionGap", "topN"],
        "clump": ["suggestive", "clumpP2", "clumpR2", "clumpKb"],
        "inspect": ["suggestive", "supportKb", "linkR2"],
        "annotate": ["radius", "annotateP"]
    }
    NULL_VARIANTS = 2000 # LD-pruned variants sampled for the null r^2 distribution
    LINK_ALPHA = 0.05 # adjusted empirical P below which leads are called linked, for --linkR2 auto
    NULL_SEED = 1

    def __init__(self, gwasDir):
        self.gwasDir = os.path.abspath(gwasDir)
        self.manifestFile = find_gwas_manifest(self.gwasDir)
        with open(self.manifestFile, "r") as fileIn:
            self.manifest = json.load(fileIn)
        self.trait = get_trait(self.manifest.get("traitType", "binary"), **self.manifest.get("traitSettings", {}))
        self.outputDir = os.path.join(self.gwasDir, "summarise")
        self.outputPrefix = os.path.join(self.outputDir, self.manifest["basename"])
        self.stateFile = self.outputPrefix + ".summary.json"
        if os.path.isfile(self.stateFile):
            with open(self.stateFile, "r") as fileIn:
                self.state = json.load(fileIn)
        else:
            self.state = { "options": {}, "steps": {} }

    def check_inputs(self):
        '''
        Checks that the run's results, phenotype file and filesets are unchanged.

        Returns:
            inputs -- a dictionary with structure like:
                      {
                          "results": resultFile, # string
                          "pheno": phenoFile, # string
                          "test": (prefix, format) # strings
                      }
        '''
        inputs = {}
        for key, description in [["results", "SAIGE results file"], ["pheno", "SAIGE phenotype file"]]:
            inputs[key] = check_file_record(self.manifest[key], self.gwasDir, description)
        for key, description in [["test", "Tested fileset file"], ["grm", "GRM fileset file"]]:
            for record in self.manifest[key]["files"]:
                check_file_record(record, self.gwasDir, description + " (re-run 'saige gwas' if it was remade)")
        inputs["test"] = (os.path.normpath(os.path.join(self.gwasDir, self.manifest["test"]["prefix"])),
                          self.manifest["test"]["format"])
        return inputs

    def resolve_options(self, given):
        '''
        Resolves summary options from those given (None if not given), then
        those recorded by an earlier summary, then defaults. Annotation file
        locations are recorded relative to the output directory.

        Parameters:
            given -- a dictionary of option names (as in DEFAULTS) to values or None
        Returns:
            options -- a dictionary of option names to resolved values, with
                       absolute annotation file locations
        '''
        recorded = self.state["options"]
        options = {}
        for key, default in SAIGESummary.DEFAULTS.items():
            if given.get(key) is not None:
                options[key] = given[key]
            elif key in recorded:
                options[key] = recorded[key]
                if key in SAIGESummary.ANNOTATION_FILES and options[key] is not None:
                    options[key] = os.path.normpath(os.path.join(self.outputDir, options[key]))
            else:
                options[key] = default

        # Title the plots with the run's covariates by default
        if options["title"] is None:
            covariates = self.manifest["pcTag"] + \
                         ("".join( f" + {v}" for v in self.manifest["covarSNPs"] ) if self.manifest["covarSNPs"] else "")
            options["title"] = f"{self.manifest['basename']} ({covariates})"

        # Annotation options other than the GFF3 need one
        if options["gff3File"] is None:
            unused = [ key for key in SAIGESummary.ANNOTATION_OPTIONS
                       if key != "gff3File" and given.get(key) is not None ]
            if unused != []:
                raise ValueError(f"Annotation option(s) {unused} need a GFF3 (--gff3), which was not given now " +
                                 "or in an earlier summary")
        elif options["annotariumDir"] is None:
            raise ValueError("Gene annotation needs the annotarium repository location (--annotarium)")
        return options

    def record_options(self, options):
        self.state["options"] = { key: (os.path.relpath(value, self.outputDir)
                                        if key in SAIGESummary.ANNOTATION_FILES and value is not None else value)
                                  for key, value in options.items() }
        self.save_state()

    def save_state(self):
        with open(self.stateFile, "w") as fileOut:
            json.dump(self.state, fileOut, indent=4)
            fileOut.write("\n")

    def is_current(self, step, params, required=[]):
        '''
        Checks whether a step was previously run with these parameters, and its
        outputs still exist, along with any required outputs (e.g., those added
        to a step since an earlier summary was made).
        '''
        previous = self.state["steps"].get(step)
        return previous is not None and previous["params"] == params and \
               all( os.path.isfile(os.path.join(self.outputDir, f)) for f in previous["outputs"] ) and \
               all( os.path.isfile(f) for f in required )

    def start_step(self, step, possibleOutputs):
        '''
        Removes the outputs of a step's previous run before it is re-run.
        '''
        previous = self.state["steps"].pop(step, { "outputs": [] })
        self.save_state()
        for fileName in set(possibleOutputs) | { os.path.join(self.outputDir, f) for f in previous["outputs"] }:
            if os.path.isfile(fileName):
                os.remove(fileName)

    def finish_step(self, step, params, outputs, results=None):
        self.state["steps"][step] = { "params": params, "outputs": [ os.path.basename(f) for f in outputs ] }
        if results is not None:
            self.state["steps"][step]["results"] = results
        self.save_state()

def configure_plink(plink, args):
    '''
    Applies the plink2 runtime arguments of 'saige summarise' to a PLINK object.
    '''
    plink.threads, plink.memory, plink.tmpDir, plink.keepTmp = args.threads, args.memory, args.tmpDir, args.keepTmp
    if args.plink2 is not None:
        plink.plinkExe = args.plink2
    return plink

def inspect_leads(summary, inputs, leads, options, args, leadsFile, familiesFile):
    '''
    Writes the lead variant inspections of 'saige summarise' (see LeadInspector).

    Returns:
        results -- a dictionary of the linkage threshold and counts, to record
    '''
    trait = summary.trait
    leadStats = parse_saige_leads(inputs["results"], leads, options["supportKb"], options["suggestive"],
                                  trait.SAIGE_RESULT_COLUMNS)
    with open(inputs["pheno"], "r") as fileIn:
        header = fileIn.readline().rstrip("\r\n").split("\t")
        phenoIndex, iidIndex = header.index(trait.PHENO_NAME), header.index("IID")
        phenos = { sl[iidIndex]: float(sl[phenoIndex]) for sl in ( line.rstrip("\r\n").split("\t") for line in fileIn ) }

    plink = configure_plink(PLINK(*inputs["test"]), args)
    grmPlink = configure_plink(PLINK(os.path.normpath(os.path.join(summary.gwasDir, summary.manifest["grm"]["prefix"])),
                                     summary.manifest["grm"]["format"]), args)
    with plink.work_dir(summary.outputPrefix) as workDir:
        leadIDs = { lead["id"] for lead in leads }
        nullVariants = [ v for v in grmPlink.thin_variants(min(SAIGESummary.NULL_VARIANTS, summary.manifest["grm"]["variants"]),
                                                           workDir, seed=SAIGESummary.NULL_SEED)
                         if v[1] not in leadIDs ]
        extractFile = os.path.join(workDir, "inspect.txt")
        with open(extractFile, "w") as fileOut:
            fileOut.write("\n".join([ lead["id"] for lead in leads ] + [ v[1] for v in nullVariants ]) + "\n")
        iids, variantIDs, countedAlleles, dosages = plink.export_dosages(extractFile, workDir)
        sampleOrder, sampleDict = parse_sample_info(plink.write_psam(os.path.join(workDir, "samples")))

    # Orient lead dosages to SAIGE's effect allele (Allele2)
    columnIndex = { variantID: i for i, variantID in enumerate(variantIDs) }
    missing = [ lead["id"] for lead in leads if lead["id"] not in columnIndex ]
    if missing != []:
        raise ValueError(f"{len(missing)} lead variant(s) are not in '{plink.prefix}' (e.g., '{missing[0]}')")
    leadDosages = np.empty((len(iids), len(leads)))
    for i, lead in enumerate(leads):
        stats, counted = leadStats[lead["id"]], countedAlleles[columnIndex[lead["id"]]]
        column = dosages[:, columnIndex[lead["id"]]]
        if counted.upper() == stats["Allele2"].upper(): # SAIGE upper-cases soft-masked alleles
            leadDosages[:, i] = column
        elif counted.upper() == stats["Allele1"].upper():
            leadDosages[:, i] = 2 - column
        else:
            raise ValueError(f"plink2 counted allele '{counted}' of '{lead['id']}' is neither SAIGE allele " +
                             f"('{stats['Allele1']}', '{stats['Allele2']}')")

    # Null distribution of r^2 between unlinked variants, among the tested samples
    testedRows = [ i for i, iid in enumerate(iids) if iid in phenos ]
    nullColumns = [ columnIndex[v[1]] for v in nullVariants if v[1] in columnIndex ]
    nullChroms = [ v[0] for v in nullVariants if v[1] in columnIndex ]
    nullR2 = inter_chromosomal_r2(dosages[np.ix_(testedRows, nullColumns)], nullChroms)

    inspector = LeadInspector(leads, leadStats, iids, leadDosages, sampleDict, phenos, trait)
    inspector.set_null(nullR2, alpha=SAIGESummary.LINK_ALPHA,
                       linkR2=None if options["linkR2"] == "auto" else options["linkR2"])
    counts = inspector.write(leadsFile, familiesFile)
    return { "nullVariants": len(nullColumns), "nullPairs": len(nullR2),
             "nullR2Quantiles": { str(q): round(float(np.quantile(nullR2, q)), 4) for q in [0.5, 0.99, 0.9999] },
             "hasParents": inspector.hasParents, **counts }

def saige_summarise(args):
    '''
    Handles "saige summarise" mode of variantopia.
    '''
    summary = SAIGESummary(args.gwasDir)
    inputs = summary.check_inputs()
    options = summary.resolve_options({ key: getattr(args, key) for key in SAIGESummary.DEFAULTS })
    os.makedirs(summary.outputDir, exist_ok=True)
    summary.record_options(options)
    outputPrefix = summary.outputPrefix
    resultsRecord = summary.manifest["results"]
    print(f"# Summarising '{inputs['results']}' ({summary.trait.TYPE} trait; {summary.manifest['pcTag']}; " +
          f"{summary.manifest['samples']} samples) into '{summary.outputDir}'")

    # Step 1: Manhattan and QQ plots, top hits and associated regions
    meff = summary.manifest["grm"]["variants"]
    plotParams = { "version": __version__, "results": resultsRecord, "meff": meff,
                   **{ key: options[key] for key in SAIGESummary.STEP_OPTIONS["plot"] } }
    plotOutputs = [ outputPrefix + suffix for suffix in [".manhattan.png", ".qq.png", ".top_hits.tsv",
                                                          ".chrom_lambda.tsv", ".plot.log"] ]
    regionsFile = outputPrefix + ".regions.tsv"
    plotRan = not summary.is_current("plot", plotParams, plotOutputs)
    if plotRan:
        summary.start_step("plot", plotOutputs + [regionsFile])
        suggestive = f"{options['suggestive']:g}"
        run_rscript("plot_manhattan.R", ["-i", inputs["results"], "-o", outputPrefix, "--sig", "bonferroni",
                                         "--meff", str(meff), "--suggestive", suggestive, "--annotate", suggestive,
                                         "--region_p", suggestive, "--region_gap", f"{options['regionGap']:g}",
                                         "--top_n", str(options["topN"]), "--title", options["title"]],
                    args.rscript, logFile=outputPrefix + ".plot.log")
        with open(outputPrefix + ".plot.log", "r") as fileIn:
            for line in fileIn:
                if not line.startswith("Program completed"):
                    print("# " + line.rstrip("\r\n"))
        summary.finish_step("plot", plotParams, plotOutputs + ([regionsFile] if os.path.isfile(regionsFile) else []))
        print(f"# Wrote plots, top hits, per-chromosome lambda and regions to {outputPrefix}.*")
    else:
        print("# Plots and regions are up to date; skipping")

    # Step 2: Clump associated variants into independent loci, using LD among the
    # samples in the null model; index variants must reach the suggestive threshold,
    # and clumped variants are those within --clumpKb and r^2 >= --clumpR2 of one
    with open(inputs["pheno"], "r") as fileIn:
        iidIndex = fileIn.readline().rstrip("\r\n").split("\t").index("IID")
        keepIIDs = [ line.rstrip("\r\n").split("\t")[iidIndex] for line in fileIn ]
    clumpParams = { "version": __version__, "results": resultsRecord, "samples": len(keepIIDs),
                    **{ key: options[key] for key in SAIGESummary.STEP_OPTIONS["clump"] } }
    clumpsFile = outputPrefix + ".clumps"
    clumpRan = not summary.is_current("clump", clumpParams)
    if clumpRan:
        summary.start_step("clump", [clumpsFile, clumpsFile + ".log"])
        plink = configure_plink(PLINK(*inputs["test"]), args)
        with plink.work_dir(outputPrefix) as workDir:
            clumped = plink.clump(inputs["results"], outputPrefix, workDir, idField="MarkerID", pField="p.value",
                                  p1=options["suggestive"], p2=options["clumpP2"], r2=options["clumpR2"],
                                  kb=options["clumpKb"], keepIIDs=keepIIDs)
        if clumped is None:
            print(f"# No variants reached P < {options['suggestive']:g}; no clumps written")
        else:
            print(f"# Independent loci (index P < {options['suggestive']:g}): {PLINK.count_lines(clumpsFile) - 1}")
        summary.finish_step("clump", clumpParams, [clumpsFile + ".log"] + ([clumpsFile] if clumped else []))
    else:
        print("# Clumps are up to date; skipping")

    # Step 3: Inspect lead variants of regions and clumps
    leads = parse_lead_variants(regionsFile if os.path.isfile(regionsFile) else None,
                                clumpsFile if os.path.isfile(clumpsFile) else None)
    inspectParams = { "version": __version__, "plot": plotParams, "clump": clumpParams, "grm": summary.manifest["grm"]["files"],
                      "nullVariants": SAIGESummary.NULL_VARIANTS, "linkAlpha": SAIGESummary.LINK_ALPHA,
                      "nullSeed": SAIGESummary.NULL_SEED,
                      **{ key: options[key] for key in SAIGESummary.STEP_OPTIONS["inspect"] } }
    leadsFile, familiesFile = outputPrefix + ".leads.tsv", outputPrefix + ".lead_families.tsv"
    if not (plotRan or clumpRan) and summary.is_current("inspect", inspectParams):
        print("# Lead variant inspection is up to date; skipping")
    else:
        summary.start_step("inspect", [leadsFile, familiesFile])
        if leads == []:
            print("# No lead variants to inspect")
            summary.finish_step("inspect", inspectParams, [])
        else:
            results = inspect_leads(summary, inputs, leads, options, args, leadsFile, familiesFile)
            quantiles = results["nullR2Quantiles"]
            print(f"# Null r^2 between {results['nullVariants']} LD-pruned variants on different chromosomes: " +
                  f"median {quantiles['0.5']}, 99% {quantiles['0.99']}, 99.99% {quantiles['0.9999']}")
            threshold = f"r^2 >= {options['linkR2']:g}" if options["linkR2"] != "auto" else \
                        f"adjusted empirical P < {SAIGESummary.LINK_ALPHA:g}"
            print(f"# Lead variants inspected: {results['leads']}; {results['linked']} are linked to a stronger " +
                  f"lead on another chromosome ({threshold})")
            if results["hasParents"]:
                print(f"# Lead variants with Mendelian errors in >= {LeadInspector.MENDEL_RATE_FLAG:.0%} of offspring " +
                      f"with a genotyped parent: {results['withMendelErrors']}")
            else:
                print("# No family has a genotyped parent; parental origin and Mendelian errors are not assessed")
            summary.finish_step("inspect", inspectParams, [leadsFile, familiesFile], results)
            print(f"# Wrote {leadsFile} and {familiesFile}")

    # Step 4: Optionally, relate lead variants of regions and clumps to nearby genes
    if options["gff3File"] is None:
        print("# No GFF3 given; skipping gene annotation")
        return
    annotateParams = { key: (file_record(options[key], summary.outputDir) if options[key] is not None else None)
                       for key in SAIGESummary.ANNOTATION_FILES }
    annotateParams.update({ key: options[key] for key in SAIGESummary.STEP_OPTIONS["annotate"] })
    annotateParams.update({ "version": __version__, "plot": plotParams, "clump": clumpParams })
    leadGenesFile = outputPrefix + ".lead_genes.tsv"
    if not (plotRan or clumpRan) and summary.is_current("annotate", annotateParams):
        print("# Gene annotation is up to date; skipping")
        return
    summary.start_step("annotate", [leadGenesFile])
    GFF3Feature, GFF3Tarium = import_annotarium_gff3(options["annotariumDir"])
    reporter = LeadGeneReporter(GFF3Tarium(options["gff3File"]),
                                import_annotarium_annotable(options["annotariumDir"])(options["annotableFile"])
                                if options["annotableFile"] is not None else None)
    reporter.chromMap = options["chromMap"]
    reporter.radius = options["radius"]
    if options["annotateP"] is not None:
        leads = [ lead for lead in leads if float(lead["p"]) < options["annotateP"] ]
    counts = reporter.write(leads, leadGenesFile)
    if counts["missingSeqids"] != []:
        print(f"# WARNING: lead variant seqid(s) {counts['missingSeqids']} have no genes in the GFF3; " +
              "check --chromMap")
    print(f"# Lead variants annotated: {counts['leads']} ({counts['withGenes']} with genes within " +
          f"{options['radius']} bp; {counts['genes']} lead variant-gene pairs)")
    summary.finish_step("annotate", annotateParams, [leadGenesFile])
    print(f"# Wrote {leadGenesFile}")
