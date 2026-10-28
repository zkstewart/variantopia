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

import os, re, sys, glob, shutil, statistics, subprocess
from concurrent.futures import ThreadPoolExecutor
from statistics import NormalDist

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from plink import PLINK, PLINKTool, parse_sample_info
from subprocess_runners import run_rscript

SAIGE_SUBMODULE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "external", "SAIGE")

class SAIGEPhenotype:
    '''
    Builds a SAIGE phenotype file (for step1_fitNULLGLMM.R --phenoFile) from a
    fileset's sample information and a PLINK2 PCA, with the chosen PCs, and
    optionally the ALT allele dosages of variants, as covariates.

    Samples are kept in sample information order, and are excluded if their
    phenotype is missing or they are absent from the PCA. A binary phenotype
    in PLINK coding (1 = control, 2 = case; 0 or -9 = missing) is recoded to 0/1.

    Example:
        pheno = SAIGEPhenotype("prefix.psam", "prefix.pca.eigenvec")
        covariates = pheno.write("prefix.pheno", [1, 2])
    '''
    MISSING_VALUES = {"-9", "0", "NA", "na", "nan", "NaN", "."}
    ID_COLUMNS = {"FID", "IID", "SID", "PAT", "MAT", "SEX"}

    def __init__(self, sampleFile, eigenvecFile):
        self.sampleFile = sampleFile
        self.eigenvecFile = eigenvecFile
        self.phenoColumn = None
        self.phenoName = "y_binary"
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
        return self._phenoName

    @phenoName.setter
    def phenoName(self, value):
        '''
        The phenotype's column name in the SAIGE phenotype file.
        '''
        if not isinstance(value, str) or value == "":
            raise ValueError("phenoName must be a non-empty string")
        self._phenoName = value

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
    def recode_binary(value):
        '''
        Recodes a PLINK binary phenotype (1/2) to 0/1, or None if missing.
        '''
        if value in SAIGEPhenotype.MISSING_VALUES:
            return None
        if value not in {"1", "2"}:
            raise ValueError(f"Phenotype value '{value}' is not a binary PLINK phenotype (1 = control, 2 = case)")
        return str(int(value) - 1)

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
        Writes the SAIGE phenotype file. Sample counts are recorded in the counts
        property.

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

        self.counts = { "written": 0, "cases": 0, "missingPheno": 0, "missingPCA": 0 }
        with open(outputFileName, "w") as fileOut:
            fileOut.write("\t".join([self.phenoName] + pcNames + ["IID"] + snpNames) + "\n")
            for iid in self.sampleOrder:
                pheno = SAIGEPhenotype.recode_binary(self.sampleDict[iid][self.phenoColumn])
                if pheno is None:
                    self.counts["missingPheno"] += 1
                    continue
                if iid not in self.pcaDict:
                    self.counts["missingPCA"] += 1
                    continue
                row = [pheno] + [ self.pcaDict[iid][i-1] for i in pcList ] + [iid] + dosageDict.get(iid, means)
                fileOut.write("\t".join(row) + "\n")
                self.counts["written"] += 1
                self.counts["cases"] += int(pheno)
        if self.counts["written"] == 0:
            raise ValueError("No samples have both a phenotype and PCs; check that IIDs match")
        return pcNames + snpNames

class SAIGE(PLINKTool):
    '''
    Runs a SAIGE single-variant association test of a binary trait, using the
    R scripts of a SAIGE source checkout (by default, variantopia's external/SAIGE
    submodule) with an Rscript whose R library has the SAIGE package installed:
        1. A sparse GRM from LD-pruned variants
        2. The null logistic mixed model, fitted with the dense GRM of the pruned
           variants, leaving one chromosome out (LOCO)
        3. Per-chromosome tests (SPA, with Firth correction for P < 0.05), run in
           parallel and merged

    SAIGE reads its genotypes as a bed fileset, so pgen and bgen inputs are
    converted within a temporary working directory. plink2 settings are
    inherited from PLINKTool. SAIGE fixes its random seed, so results are
    reproducible.

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
        self.phenoName = "y_binary"
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
        Fits the null logistic mixed model with LOCO. A single variance ratio is
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
                f"--phenoCol={self.phenoName}", "--sampleIDColinphenoFile=IID", "--traitType=binary",
                "--isCateVarianceRatio=FALSE", "--useSparseGRMtoFitNULL=FALSE",
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
        step 1 log; SAIGE logs it as 'Final <dispersion> <tau> :'.
        '''
        tau = None
        with open(nullLog, "r") as fileIn:
            for line in fileIn:
                if line.startswith("Final "):
                    tau = float(line.split()[2])
        if tau is None:
            raise ValueError(f"No final tau estimate found in SAIGE null model log '{nullLog}'")
        return tau

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
                                        f"--minMAC={self.minMAC}", "--is_Firth_beta=TRUE",
                                        "--pCutoffforFirth=0.05", "--is_output_moreDetails=TRUE"],
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
        pheno = SAIGEPhenotype(plink.write_psam(os.path.join(workDir, "samples")), args.eigenvecFile)
        pheno.phenoColumn = args.phenoColumn

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
        pheno = SAIGEPhenotype(plink.write_psam(os.path.join(workDir, "samples")), args.eigenvecFile)
        pheno.phenoColumn = args.phenoColumn
        rawFile = saige.export_dosages(plink, args.covarSNPs, workDir) if args.covarSNPs is not None else None
        covariates = pheno.write(outputPrefix + ".pheno", pcList, rawFile)
        counts = pheno.counts
        print(f"# Samples: {counts['written']} ({counts['cases']} cases, {counts['written'] - counts['cases']} controls)")
        print(f"# Excluded {counts['missingPheno']} samples with a missing phenotype, and " +
              f"{counts['missingPCA']} phenotyped samples absent from the PCA")
        print(f"# Covariates: {','.join(covariates) or 'none'}")

        # Fit the null model with the pruned variants, then test all variants
        grmBed = saige.as_bed(grmPlink, workDir)
        sparseGRM = saige.sparse_grm(grmBed, outputPrefix, workDir)
        tau = saige.fit_null(grmBed, sparseGRM, outputPrefix + ".pheno", covariates, outputPrefix)
        print(f"# Null model GRM variance (tau): {tau:.7g}")
        testBed = saige.as_bed(plink, workDir)
        numTested = saige.test(testBed, sparseGRM, outputPrefix, resultFile, workDir)
    print(f"# Variants tested: {numTested}")
    print(f"# Wrote {resultFile}")
