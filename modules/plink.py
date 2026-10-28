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

import os, re, sys, shutil, statistics, subprocess, tempfile
from collections import Counter
from contextlib import contextmanager

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, parse_2col_tsv_as_dict
from subprocess_runners import run_rscript

FAM_COLUMNS = ["FID", "IID", "PAT", "MAT", "SEX", "PHENO1"]

def parse_sample_info(sampleFile):
    '''
    Parses a .psam (with a '#'-prefixed header) or .fam file.

    Parameters:
        sampleFile -- a string pointing to a .psam or .fam file
    Returns:
        sampleOrder -- a list of IIDs in file order
        sampleDict -- a dictionary with structure like:
                      {
                          iid: {column: value, ...},
                          ...
                      }
                      where .fam columns are named as in a variantopia .psam
                      (FID, IID, PAT, MAT, SEX, PHENO1)
    '''
    sampleOrder, sampleDict = [], {}
    with open(sampleFile, "r") as fileIn:
        columns = FAM_COLUMNS
        for line in fileIn:
            if line.startswith("#"):
                columns = line.rstrip("\r\n").lstrip("#").split("\t")
                continue
            sl = line.split()
            if sl == []:
                continue
            row = dict(zip(columns, sl))
            if row["IID"] in sampleDict:
                raise ValueError(f"Duplicate IID '{row['IID']}' in '{sampleFile}'")
            sampleOrder.append(row["IID"])
            sampleDict[row["IID"]] = row
    return sampleOrder, sampleDict

class PLINKTool:
    '''
    Settings and behaviour shared by all plink2-running classes: the plink2
    executable and its resources, a temporary working directory for intermediate
    files, and the format of any PLINK2 fileset that is written.

    Properties are given sensible defaults and can be adjusted prior to calling
    the methods of child classes.
    '''
    OUTPUT_FORMATS = ["pgen", "bgen", "bed"]
    COMPRESSION_TYPES = { # first value is the default for that format
        "pgen": ["zstd", None], # zstd compresses the .pvar ('vzs')
        "bgen": ["zlib", "zstd"], # zlib == bgen v1.2, zstd == bgen v1.3
        "bed": [None]
    }
    BGEN_VERSIONS = { "zlib": "1.2", "zstd": "1.3" }

    def __init__(self):
        self.outputFormat = "pgen" # also sets the default compression
        self.bgenBits = 8
        self.threads = 1
        self.memory = 16000
        self.tmpDir = None
        self.keepTmp = False
        self._plinkExe = shutil.which("plink2") # may be None; validated when plink2 is first run

    @property
    def outputFormat(self):
        return self._outputFormat

    @outputFormat.setter
    def outputFormat(self, value):
        if value not in PLINKTool.OUTPUT_FORMATS:
            raise ValueError(f"outputFormat must be one of {PLINKTool.OUTPUT_FORMATS}")
        self._outputFormat = value

        # Reset compression to the format default if it is now incompatible
        if getattr(self, "_compression", "unset") not in PLINKTool.COMPRESSION_TYPES[value]:
            self._compression = PLINKTool.COMPRESSION_TYPES[value][0]

    @property
    def compression(self):
        return self._compression

    @compression.setter
    def compression(self, value):
        '''
        For pgen, 'zstd' (Zstd-compressed .pvar) or None; for bgen, 'zlib'
        (bgen v1.2) or 'zstd' (bgen v1.3); for bed, None.
        '''
        allowed = PLINKTool.COMPRESSION_TYPES[self.outputFormat]
        if value not in allowed:
            raise ValueError(f"compression for '{self.outputFormat}' output must be one of {allowed}")
        self._compression = value

    @property
    def bgenBits(self):
        return self._bgenBits

    @bgenBits.setter
    def bgenBits(self, value):
        '''
        Bits of precision per probability in bgen output; 8 is lossless for
        hardcalls and is widely supported (e.g., by SAIGE and REGENIE).
        '''
        if not isinstance(value, int) or not 1 <= value <= 32:
            raise ValueError("bgenBits must be an integer from 1 to 32")
        self._bgenBits = value

    @property
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, value):
        if not isinstance(value, int) or value < 1:
            raise ValueError("threads must be a positive integer")
        self._threads = value

    @property
    def memory(self):
        return self._memory

    @memory.setter
    def memory(self, value):
        '''
        Memory for plink2 in megabytes, as an integer (or integer-like string).
        '''
        try:
            value = int(value)
        except (TypeError, ValueError):
            raise ValueError("memory must be an integer number of megabytes")
        if value < 1:
            raise ValueError("memory must be a positive integer")
        self._memory = value

    @property
    def tmpDir(self):
        return self._tmpDir

    @tmpDir.setter
    def tmpDir(self, value):
        '''
        Parent directory for the temporary working directory; None (default)
        uses the output directory, since intermediate files can be large.
        '''
        if value is not None:
            if not isinstance(value, str):
                raise TypeError("tmpDir must be None or a string")
            if not os.path.isdir(value):
                raise FileNotFoundError(f"Temporary directory location '{value}' does not exist")
            value = os.path.abspath(value)
        self._tmpDir = value

    @property
    def keepTmp(self):
        return self._keepTmp

    @keepTmp.setter
    def keepTmp(self, value):
        if not isinstance(value, bool):
            raise TypeError("keepTmp must be a boolean")
        self._keepTmp = value

    @property
    def plinkExe(self):
        return self._plinkExe

    @plinkExe.setter
    def plinkExe(self, value):
        if not isinstance(value, str):
            raise TypeError("plinkExe must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"plink2 executable not found at '{value}'")
        self._plinkExe = os.path.abspath(value)

    def output_files(self, outputPrefix):
        '''
        The files of a fileset written by write_fileset() to outputPrefix for the
        current outputFormat and compression.
        '''
        if self.outputFormat == "pgen":
            pvarSuffix = ".pvar.zst" if self.compression == "zstd" else ".pvar"
            suffixes = [".pgen", pvarSuffix, ".psam"]
        elif self.outputFormat == "bgen":
            suffixes = [".bgen", ".sample"]
        else:
            suffixes = [".bed", ".bim", ".fam"]
        return [ outputPrefix + suffix for suffix in suffixes ]

    def check_outputs_absent(self, outputFiles):
        existing = [ f for f in outputFiles if os.path.exists(f) ]
        if existing != []:
            raise FileExistsError(f"Output file(s) already exist: {existing}")

    def run_plink(self, args, outPrefix):
        '''
        Runs plink2 with the given arguments, adding --threads, --memory,
        --allow-extra-chr and --out.

        Parameters:
            args -- a list of strings of plink2 arguments
            outPrefix -- a string of the --out prefix for this run
        '''
        if self.plinkExe is None:
            raise FileNotFoundError("plink2 is not discoverable in your system PATH and plinkExe was not set")
        cmd = [self.plinkExe, "--threads", str(self.threads), "--memory", str(self.memory),
               "--allow-extra-chr"] + args + ["--out", outPrefix]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            logTail = ""
            if os.path.isfile(outPrefix + ".log"):
                with open(outPrefix + ".log", "r") as logIn:
                    logTail = "".join(logIn.readlines()[-15:])
            raise Exception(f"plink2 failed with exit code {result.returncode}; command was '{' '.join(cmd)}'. " +
                            f"Log tail:\n{logTail}\nstderr:\n{result.stderr.decode('utf-8')}")

    @contextmanager
    def work_dir(self, outputPrefix):
        '''
        Context manager yielding a temporary working directory, created within
        tmpDir (or alongside outputPrefix) and removed afterwards unless keepTmp.
        '''
        workDir = tempfile.mkdtemp(prefix=".variantopia_plink_",
                                   dir=self.tmpDir or os.path.dirname(outputPrefix))
        try:
            yield workDir
        finally:
            if not self.keepTmp:
                shutil.rmtree(workDir, ignore_errors=True)

    def write_fileset(self, args, outputPrefix, workDir, sortVars=False):
        '''
        Writes a fileset in the configured outputFormat and compression.

        Parameters:
            args -- a list of strings of plink2 arguments specifying the input
                    and any filters to apply
            outputPrefix -- a string of the output fileset's prefix
            workDir -- a string pointing to the working directory; only used
                       when sorting variants for bgen or bed output
            sortVars -- a boolean indicating whether variants should be sorted
        Returns:
            outputFiles -- a list of the output files written
        '''
        if self.outputFormat == "pgen":
            makeArgs = (["--sort-vars"] if sortVars else []) + ["--make-pgen"] + \
                       (["vzs"] if self.compression == "zstd" else [])
            self.run_plink(args + makeArgs, outputPrefix)
        else:
            # --sort-vars requires --make-pgen, so sort into the working directory first
            if sortVars:
                sortedPrefix = os.path.join(workDir, "sorted")
                self.run_plink(args + ["--sort-vars", "--make-pgen"], sortedPrefix)
                args = ["--pfile", sortedPrefix]
            if self.outputFormat == "bgen":
                self.run_plink(args + ["--export", f"bgen-{PLINKTool.BGEN_VERSIONS[self.compression]}",
                                       f"bits={self.bgenBits}", "ref-first"], outputPrefix)
            else:
                self.run_plink(args + ["--make-bed"], outputPrefix)
        return self.output_files(outputPrefix)

    @staticmethod
    def parse_log_counts(logFile):
        '''
        Parses the numbers of samples and variants remaining after filtering from
        a plink2 log file.

        Parameters:
            logFile -- a string pointing to a plink2 .log file
        Returns:
            counts -- a dictionary with structure like:
                      {
                          "samples": int,
                          "variants": int
                      }
        '''
        with open(logFile, "r") as fileIn:
            log = " ".join(fileIn.read().split()) # plink2 wraps long lines
        counts = {}
        for key, pattern in [["samples", r"(\d+) samples \([^)]*\) remaining after main filters"],
                             ["variants", r"(\d+) variants remaining after main filters"]]:
            match = re.findall(pattern, log)
            if match != []:
                counts[key] = int(match[-1])
        return counts

class PLINK(PLINKTool):
    '''
    Encapsulates an existing PLINK fileset (pgen, bgen or bed) for QC reporting,
    sample and variant filtering, LD pruning, PCA, and KING kinship estimation.

    Sample and variant selection properties (removeFile, requirePheno,
    nonfounders, snpsOnly, geno, maf) apply to both filter() and prune(); they
    default to no filtering and can be adjusted prior to calling either method.
    Pruning properties are given sensible defaults. Intermediate files are written to a temporary directory which is
    removed afterwards (unless keepTmp is True), so only the final outputs are
    written to the output location.

    Example:
        plink = PLINK("input/prefix", "pgen")
        plink.sample_qc("output/prefix")
        plink.nonfounders = True
        plink.requirePheno = True
        plink.geno, plink.maf = 0.1, 0.01
        plink.filter("output/prefix.qc")
        plink.snpsOnly = True
        plink.geno, plink.maf = 0.05, 0.05
        plink.densityCap = 3
        plink.prune("output/prefix.pruned")
        PLINK("output/prefix.pruned", "pgen").pca("output/pca/prefix.pruned.pca")
        plink.king("output/king/prefix.king", extractFile="output/prefix.pruned.pvar.zst")
    '''
    INPUT_FORMATS = ["pgen", "bgen", "bed"]

    def __init__(self, prefix, fileFormat="pgen"):
        super().__init__()
        self._set_fileset(prefix, fileFormat)

        # Set default values for adjustable sample and variant selection properties
        self.removeFile = None
        self.requirePheno = False
        self.nonfounders = False
        self.snpsOnly = False
        self.geno = None
        self.maf = None

        # Set default values for adjustable pruning properties
        self.r2 = 0.2
        self.windowVariants = 50
        self.stepVariants = 5
        self.windowKb = 1000
        self.densityCap = None
        self.densityBinSize = 1000000

        # Set default values for adjustable PCA properties
        self.npcs = 20
        self.filterCounts = {}
        self.pruneCounts = {}

    def _set_fileset(self, prefix, fileFormat):
        '''
        Validates that the fileset exists, noting whether a pgen fileset has a
        Zstd-compressed .pvar.
        '''
        if fileFormat not in PLINK.INPUT_FORMATS:
            raise ValueError(f"fileFormat must be one of {PLINK.INPUT_FORMATS}")
        if not isinstance(prefix, str):
            raise TypeError("prefix must be a string")
        prefix = os.path.abspath(prefix)

        inputCompression = None
        if fileFormat == "pgen":
            if os.path.isfile(prefix + ".pvar.zst"):
                inputCompression = "zstd"
                suffixes = [".pgen", ".pvar.zst", ".psam"]
            else:
                suffixes = [".pgen", ".pvar", ".psam"]
        elif fileFormat == "bgen":
            suffixes = [".bgen", ".sample"]
        else:
            suffixes = [".bed", ".bim", ".fam"]
        missing = [ prefix + s for s in suffixes if not os.path.isfile(prefix + s) ]
        if missing != []:
            raise FileNotFoundError(f"{fileFormat} fileset '{prefix}' is missing file(s): {missing}")

        self._prefix = prefix
        self._fileFormat = fileFormat
        self._inputCompression = inputCompression

    @property
    def prefix(self):
        return self._prefix

    @property
    def fileFormat(self):
        return self._fileFormat

    @property
    def inputArgs(self):
        '''
        plink2 arguments for reading this fileset.
        '''
        if self.fileFormat == "pgen":
            args = ["--pfile", self.prefix]
            if self._inputCompression == "zstd":
                args.append("vzs")
        elif self.fileFormat == "bgen":
            args = ["--bgen", self.prefix + ".bgen", "ref-first",
                    "--sample", self.prefix + ".sample"]
        else:
            args = ["--bfile", self.prefix]
        return args

    @property
    def removeFile(self):
        return self._removeFile

    @removeFile.setter
    def removeFile(self, value):
        '''
        None, or a file of samples to exclude (one IID per line, or FID IID).
        '''
        if value is not None:
            if not isinstance(value, str):
                raise TypeError("removeFile must be None or a string")
            if not os.path.isfile(value):
                raise FileNotFoundError(f"Sample removal file '{value}' is not a file")
            value = os.path.abspath(value)
        self._removeFile = value

    @property
    def requirePheno(self):
        return self._requirePheno

    @requirePheno.setter
    def requirePheno(self, value):
        '''
        If True, samples with a missing phenotype are excluded.
        '''
        if not isinstance(value, bool):
            raise TypeError("requirePheno must be a boolean")
        self._requirePheno = value

    @property
    def nonfounders(self):
        return self._nonfounders

    @nonfounders.setter
    def nonfounders(self, value):
        '''
        If True, all samples are treated as founders for allele frequencies (and
        hence MAF filtering) and LD calculation. This suits family-based designs
        where most samples have parents recorded in the fileset.
        '''
        if not isinstance(value, bool):
            raise TypeError("nonfounders must be a boolean")
        self._nonfounders = value

    @property
    def snpsOnly(self):
        return self._snpsOnly

    @snpsOnly.setter
    def snpsOnly(self, value):
        '''
        If True, only SNPs with A/C/G/T alleles are retained; this suits pruning,
        since indel calls are noisier.
        '''
        if not isinstance(value, bool):
            raise TypeError("snpsOnly must be a boolean")
        self._snpsOnly = value

    @property
    def geno(self):
        return self._geno

    @geno.setter
    def geno(self, value):
        '''
        Maximum variant missingness for a variant to be retained, or None for
        no filter.
        '''
        if value is not None and (not isinstance(value, (int, float)) or not 0 <= value <= 1):
            raise ValueError("geno must be None or a number from 0 to 1")
        self._geno = value

    @property
    def maf(self):
        return self._maf

    @maf.setter
    def maf(self, value):
        '''
        Minimum minor allele frequency for a variant to be retained, or None for
        no filter.
        '''
        if value is not None and (not isinstance(value, (int, float)) or not 0 <= value <= 0.5):
            raise ValueError("maf must be None or a number from 0 to 0.5")
        self._maf = value

    @property
    def r2(self):
        return self._r2

    @r2.setter
    def r2(self, value):
        '''
        r^2 threshold for both pruning stages.
        '''
        if not isinstance(value, (int, float)) or not 0 < value <= 1:
            raise ValueError("r2 must be a number greater than 0 and at most 1")
        self._r2 = value

    @property
    def windowVariants(self):
        return self._windowVariants

    @windowVariants.setter
    def windowVariants(self, value):
        '''
        Window size (in variants) for the first, local pruning stage.
        '''
        if not isinstance(value, int) or value < 2:
            raise ValueError("windowVariants must be an integer >= 2")
        self._windowVariants = value

    @property
    def stepVariants(self):
        return self._stepVariants

    @stepVariants.setter
    def stepVariants(self, value):
        '''
        Window step size (in variants) for the first, local pruning stage.
        '''
        if not isinstance(value, int) or value < 1:
            raise ValueError("stepVariants must be a positive integer")
        self._stepVariants = value

    @property
    def windowKb(self):
        return self._windowKb

    @windowKb.setter
    def windowKb(self, value):
        '''
        Window size (in kb) for the second, wide physical pruning stage, or None
        to skip it. LD in family-based populations extends over long distances
        that the first stage's variant-count windows do not span.
        '''
        if value is not None and (not isinstance(value, int) or value < 1):
            raise ValueError("windowKb must be None or a positive integer")
        self._windowKb = value

    @property
    def densityCap(self):
        return self._densityCap

    @densityCap.setter
    def densityCap(self, value):
        '''
        None to skip density capping, or a multiple of the genome-wide median
        number of pruned SNPs per bin; bins holding more than this are thinned
        to it by keeping evenly spaced SNPs. Regions with noisy genotypes can
        retain many pruned SNPs whose pairwise r^2 is low but which jointly
        dominate a PC.
        '''
        if value is not None and (not isinstance(value, (int, float)) or value <= 0):
            raise ValueError("densityCap must be None or a number greater than 0")
        self._densityCap = value

    @property
    def densityBinSize(self):
        return self._densityBinSize

    @densityBinSize.setter
    def densityBinSize(self, value):
        '''
        Bin size (in bp) used for density capping.
        '''
        if not isinstance(value, int) or value < 1:
            raise ValueError("densityBinSize must be a positive integer")
        self._densityBinSize = value

    @property
    def npcs(self):
        return self._npcs

    @npcs.setter
    def npcs(self, value):
        '''
        Number of principal components to compute.
        '''
        if not isinstance(value, int) or value < 1:
            raise ValueError("npcs must be a positive integer")
        self._npcs = value

    @staticmethod
    def count_lines(fileName, skipHeader=False):
        with open(fileName, "r") as fileIn:
            return sum( 1 for line in fileIn if not (skipHeader and line.startswith("#")) )

    @staticmethod
    def clear_parents(psamFile, outputFileName):
        '''
        Writes a copy of a .psam with paternal and maternal IDs set to '0', so
        that plink2 treats every sample as a founder. LD pruning has no
        --nonfounders override, and requires at least 50 founders.

        Parameters:
            psamFile -- a string pointing to the .psam file
            outputFileName -- a string of the location to write the copy
        '''
        with open(psamFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
            header = fileIn.readline()
            fileOut.write(header)
            columns = header.rstrip("\r\n").lstrip("#").split("\t")
            parentIndices = [ columns.index(c) for c in ["PAT", "MAT"] if c in columns ]
            for line in fileIn:
                sl = line.rstrip("\r\n").split("\t")
                for i in parentIndices:
                    sl[i] = "0"
                fileOut.write("\t".join(sl) + "\n")

    def density_cap(self, pvarFile, pruneInFile, outputFileName):
        '''
        Thins any bin holding more than densityCap x the median bin count of
        pruned variants down to that cap, keeping evenly spaced variants so the
        result is deterministic.

        Parameters:
            pvarFile -- a string pointing to an uncompressed .pvar file which
                        contains (at least) the pruned variants
            pruneInFile -- a string pointing to a file of pruned variant IDs
            outputFileName -- a string of the location to write the IDs of
                              variants retained after capping
        Returns:
            medianBin -- an integer of the median number of variants per bin
            binCap -- an integer of the maximum number of variants per bin
        '''
        with open(pruneInFile, "r") as fileIn:
            pruned = set( line.strip() for line in fileIn if line.strip() != "" )

        # Assign pruned variants to bins in .pvar order
        variantBins = []
        with open(pvarFile, "r") as fileIn:
            for line in fileIn:
                if line.startswith("#"):
                    continue
                chrom, pos, varId = line.split("\t", 3)[:3]
                if varId in pruned:
                    variantBins.append([varId, f"{chrom}:{int(pos) // self.densityBinSize}"])

        binCounts = Counter( b for _, b in variantBins )
        sortedCounts = sorted(binCounts.values())
        medianBin = sortedCounts[(len(sortedCounts) + 1) // 2 - 1] # lower median if even
        binCap = int(medianBin * self.densityCap)

        # Keep a variant when it crosses an integer step of cap/n along its bin
        seen = Counter()
        with open(outputFileName, "w") as fileOut:
            for varId, b in variantBins:
                seen[b] += 1
                n, i = binCounts[b], seen[b]
                if n <= binCap or (i * binCap) // n != ((i - 1) * binCap) // n:
                    fileOut.write(varId + "\n")
        return medianBin, binCap

    def selection_args(self):
        '''
        plink2 arguments for reading this fileset with the configured sample and
        variant selection applied.
        '''
        args = self.inputArgs + (["--nonfounders"] if self.nonfounders else [])
        if self.removeFile is not None:
            args += ["--remove", self.removeFile]
        if self.requirePheno:
            args.append("--require-pheno")
        if self.snpsOnly:
            args += ["--snps-only", "just-acgt"]
        if self.geno is not None:
            args += ["--geno", str(self.geno)]
        if self.maf is not None:
            args += ["--maf", str(self.maf)]
        return args

    def check_output_prefix(self, outputPrefix):
        outputPrefix = os.path.abspath(outputPrefix)
        if not os.path.isdir(os.path.dirname(outputPrefix)):
            raise FileNotFoundError(f"Output directory '{os.path.dirname(outputPrefix)}' does not exist")
        self.check_outputs_absent(self.output_files(outputPrefix))
        return outputPrefix

    def filter(self, outputPrefix):
        '''
        Writes a fileset of the selected samples (removeFile, requirePheno) and
        variants (snpsOnly, geno, maf), e.g., as the variant set for association
        testing. Retained sample and variant counts are recorded in the
        filterCounts property.

        Parameters:
            outputPrefix -- a string of the output fileset's prefix
        Returns:
            outputFiles -- a list of the output files written
        '''
        outputPrefix = self.check_output_prefix(outputPrefix)
        outputFiles = self.write_fileset(self.selection_args(), outputPrefix, workDir=None)
        self.filterCounts = PLINKTool.parse_log_counts(outputPrefix + ".log")
        return outputFiles

    def prune(self, outputPrefix):
        '''
        Writes a fileset of LD-pruned, high-quality variants (e.g., for a GRM,
        PCA, or kinship estimation). Samples and variants are selected as for
        filter() (typically with snpsOnly and stricter geno and maf) then pruned in
        two stages (a local window of windowVariants, then a wide window of
        windowKb), and optionally density capped. Per-stage variant counts are
        recorded in the pruneCounts property.

        Parameters:
            outputPrefix -- a string of the output fileset's prefix
        Returns:
            outputFiles -- a list of the output files written
        '''
        outputPrefix = self.check_output_prefix(outputPrefix)
        nonfounderArgs = ["--nonfounders"] if self.nonfounders else []

        with self.work_dir(outputPrefix) as workDir:
            # Select samples and pruning candidate variants
            filteredPrefix = os.path.join(workDir, "filtered")
            self.run_plink(self.selection_args() + ["--make-pgen"], filteredPrefix)
            self.pruneCounts = { "samples": PLINK.count_lines(filteredPrefix + ".psam", skipHeader=True),
                                 "filtered": PLINK.count_lines(filteredPrefix + ".pvar", skipHeader=True) }

            # Prune using a .psam without parents if all samples are to be founders;
            # the output retains the original pedigree information
            psamFile = filteredPrefix + ".psam"
            if self.nonfounders:
                psamFile = os.path.join(workDir, "founders.psam")
                PLINK.clear_parents(filteredPrefix + ".psam", psamFile)
            ldArgs = ["--pgen", filteredPrefix + ".pgen", "--pvar", filteredPrefix + ".pvar",
                      "--psam", psamFile] + nonfounderArgs

            stage1Prefix = os.path.join(workDir, "stage1")
            self.run_plink(ldArgs + ["--indep-pairwise", str(self.windowVariants),
                                     str(self.stepVariants), str(self.r2)], stage1Prefix)
            pruneIn = stage1Prefix + ".prune.in"
            self.pruneCounts["stage1"] = PLINK.count_lines(pruneIn)

            if self.windowKb is not None:
                stage2Prefix = os.path.join(workDir, "stage2")
                self.run_plink(ldArgs + ["--extract", pruneIn, "--indep-pairwise",
                                         f"{self.windowKb}kb", "1", str(self.r2)], stage2Prefix)
                pruneIn = stage2Prefix + ".prune.in"
                self.pruneCounts["stage2"] = PLINK.count_lines(pruneIn)

            if self.densityCap is not None:
                cappedFile = os.path.join(workDir, "capped.prune.in")
                medianBin, binCap = self.density_cap(filteredPrefix + ".pvar", pruneIn, cappedFile)
                pruneIn = cappedFile
                self.pruneCounts["densityMedianBin"] = medianBin
                self.pruneCounts["densityBinCap"] = binCap
                self.pruneCounts["densityCapped"] = PLINK.count_lines(pruneIn)

            return self.write_fileset(["--pfile", filteredPrefix, "--extract", pruneIn],
                                      outputPrefix, workDir)

    def write_psam(self, outputPrefix):
        '''
        Writes this fileset's sample information as a .psam, whatever its format.

        Parameters:
            outputPrefix -- a string of the output prefix; '<outputPrefix>.psam'
                            is written
        Returns:
            psamFile -- a string of the .psam location
        '''
        self.run_plink(self.inputArgs + ["--make-just-psam"], outputPrefix)
        return outputPrefix + ".psam"

    def pca(self, outputPrefix):
        '''
        Runs a PCA of all samples, with SNP loadings. The fileset should already
        be LD pruned (see prune()). If nonfounders is True, all samples are used
        for the allele frequencies that standardise genotypes.

        Parameters:
            outputPrefix -- a string of the output prefix; files are written to
                            '<outputPrefix>.{eigenvec,eigenval,eigenvec.var,log}'
        Returns:
            outputFiles -- a list of the output files written
        '''
        outputPrefix = os.path.abspath(outputPrefix)
        outputFiles = [ outputPrefix + suffix for suffix in [".eigenvec", ".eigenval", ".eigenvec.var"] ]
        self.check_outputs_absent(outputFiles)

        args = self.inputArgs + (["--nonfounders"] if self.nonfounders else [])
        self.run_plink(args + ["--pca", str(self.npcs), "biallelic-var-wts"], outputPrefix)
        return outputFiles

    def variant_ids(self, extractFile, workDir):
        '''
        Writes the variant IDs listed in, or described by, a file into workDir
        for use with plink2 --extract.

        Parameters:
            extractFile -- a string pointing to a .pvar, .pvar.zst or .bim file,
                           or a plain list of variant IDs (one per line)
            workDir -- a string pointing to the working directory
        Returns:
            idFile -- a string of the variant ID list's location
        '''
        # Decompress .zst files with plink2 itself
        if extractFile.endswith(".zst"):
            decompressedFile = os.path.join(workDir, os.path.basename(extractFile)[:-4])
            result = subprocess.run([self.plinkExe, "--zst-decompress", extractFile, decompressedFile],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode != 0:
                raise Exception(f"plink2 failed to decompress '{extractFile}'; stderr was:\n{result.stderr.decode('utf-8')}")
            extractFile = decompressedFile

        if not (extractFile.endswith(".pvar") or extractFile.endswith(".bim")):
            return extractFile # assume a plain ID list

        idFile = os.path.join(workDir, "extract.txt")
        with open(extractFile, "r") as fileIn, open(idFile, "w") as fileOut:
            idIndex = 1 # .bim column
            for line in fileIn:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    idIndex = line.rstrip("\r\n").split("\t").index("ID")
                    continue
                fileOut.write(line.split()[idIndex] + "\n")
        return idFile

    def king(self, outputPrefix, extractFile=None):
        '''
        Estimates KING-robust kinship between all pairs of samples. LD-pruned
        variants should be used, either by running this on a pruned fileset or
        by giving the pruned variants as extractFile; the latter allows samples
        excluded from pruning (e.g., unphenotyped parents) to be included.

        Parameters:
            outputPrefix -- a string of the output prefix; files are written to
                            '<outputPrefix>.{kin0,log}'
            extractFile -- (OPTIONAL) a string pointing to a .pvar, .pvar.zst or
                           .bim file, or a plain list of variant IDs, to restrict
                           the variants used
        Returns:
            kin0File -- a string of the kinship table's location
        '''
        outputPrefix = os.path.abspath(outputPrefix)
        kin0File = outputPrefix + ".kin0"
        self.check_outputs_absent([kin0File])

        with self.work_dir(outputPrefix) as workDir:
            args = list(self.inputArgs)
            if extractFile is not None:
                args += ["--extract", self.variant_ids(os.path.abspath(extractFile), workDir)]
            self.run_plink(args + ["--make-king-table"], outputPrefix)
        return kin0File

    def sample_qc(self, outputPrefix):
        '''
        Writes a per-sample QC report of missingness (all variants) and
        heterozygosity rate (SNPs only, since indel calls are noisier), alongside
        the plink2 reports it was derived from.

        Parameters:
            outputPrefix -- a string of the output prefix; files are written to
                            '<outputPrefix>.sample_qc.{tsv,smiss,scount}'
        Returns:
            outputFiles -- a list of the output files written
        '''
        outputPrefix = os.path.abspath(outputPrefix) + ".sample_qc"
        outputFiles = [ outputPrefix + suffix for suffix in [".tsv", ".smiss", ".scount"] ]
        self.check_outputs_absent(outputFiles)

        with self.work_dir(outputPrefix) as workDir:
            qcPrefix = os.path.join(workDir, "sample_qc")
            self.run_plink(self.inputArgs + ["--missing", "sample-only"], qcPrefix)
            self.run_plink(self.inputArgs + ["--snps-only", "just-acgt", "--sample-counts"], qcPrefix)

            # Combine the plink2 reports
            def parse_report(fileName):
                with open(fileName, "r") as fileIn:
                    header = fileIn.readline().rstrip("\r\n").lstrip("#").split("\t")
                    return [ dict(zip(header, line.rstrip("\r\n").split("\t"))) for line in fileIn ]
            smiss = parse_report(qcPrefix + ".smiss")
            scount = { row["IID"]: row for row in parse_report(qcPrefix + ".scount") }

            with open(outputFiles[0], "w") as fileOut:
                fileOut.write("FID\tIID\tF_MISS\tHET_RATE\n")
                for row in smiss:
                    counts = scount[row["IID"]]
                    called = int(counts["HOM_REF_CT"]) + int(counts["HOM_ALT_SNP_CT"]) + int(counts["HET_SNP_CT"])
                    hetRate = f"{int(counts['HET_SNP_CT']) / called:.6g}" if called > 0 else "NA"
                    fileOut.write(f"{row.get('FID', row['IID'])}\t{row['IID']}\t{row['F_MISS']}\t{hetRate}\n")
            shutil.move(qcPrefix + ".smiss", outputFiles[1])
            shutil.move(qcPrefix + ".scount", outputFiles[2])

        return outputFiles

class KING:
    '''
    Interprets a PLINK2 KING-robust kinship table (.kin0) alongside the recorded
    pedigree of its samples, to verify recorded parent-offspring relationships.

    Expected values for a true parent-offspring pair are kinship ~0.25 (first-degree
    range 0.177-0.354) and IBS0 ~0 (an offspring always shares >= 1 allele with its
    parent). Genotyping error from low-depth calls inflates IBS0 and deflates
    kinship, so the observed values are also summarised as the median across all
    recorded pairs.

    Example:
        king = KING("prefix.king.kin0", "prefix.psam")
        king.minKinship = 0.177
        stats = king.check_pedigree("prefix.king.pedigree_check.tsv")
    '''
    PEDIGREE_CHECK_HEADER = ["IID", "FID", "ROLE", "RECORDED_PARENT", "KINSHIP", "IBS0",
                             "STATUS", "BEST_CANDIDATE", "BEST_KINSHIP", "BEST_IBS0"]

    def __init__(self, kin0File, sampleFile):
        self.kin0File = kin0File
        self.sampleFile = sampleFile
        self.minKinship = 0.177

    @property
    def kin0File(self):
        return self._kin0File

    @kin0File.setter
    def kin0File(self, value):
        if not isinstance(value, str):
            raise TypeError("kin0File must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"KING kinship table '{value}' is not a file")
        self._kin0File = os.path.abspath(value)
        self.kinDict = KING.parse_kin0(self._kin0File)

    @property
    def sampleFile(self):
        return self._sampleFile

    @sampleFile.setter
    def sampleFile(self, value):
        '''
        A .psam or .fam file recording each sample's family and parents.
        '''
        if not isinstance(value, str):
            raise TypeError("sampleFile must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"Sample file '{value}' is not a file")
        self._sampleFile = os.path.abspath(value)
        self.sampleDict = KING.parse_samples(self._sampleFile)

    @property
    def minKinship(self):
        return self._minKinship

    @minKinship.setter
    def minKinship(self, value):
        '''
        Minimum kinship for a recorded parent-offspring pair to pass; 0.177 is
        the lower bound of the first-degree range.
        '''
        if not isinstance(value, (int, float)) or not 0 < value <= 0.5:
            raise ValueError("minKinship must be a number greater than 0 and at most 0.5")
        self._minKinship = value

    @property
    def genotyped(self):
        return { iid for pair in self.kinDict for iid in pair }

    @staticmethod
    def parse_kin0(kin0File):
        '''
        Parameters:
            kin0File -- a string pointing to a PLINK2 .kin0 file
        Returns:
            kinDict -- a dictionary with structure like:
                       {
                           frozenset({iid1, iid2}): (kinship, ibs0),
                           ...
                       }
        '''
        kinDict = {}
        with open(kin0File, "r") as fileIn:
            header = fileIn.readline().rstrip("\r\n").lstrip("#").split("\t")
            id1, id2 = header.index("IID1"), header.index("IID2")
            kinIndex, ibs0Index = header.index("KINSHIP"), header.index("IBS0")
            for line in fileIn:
                sl = line.rstrip("\r\n").split("\t")
                kinDict[frozenset((sl[id1], sl[id2]))] = (float(sl[kinIndex]), float(sl[ibs0Index]))
        return kinDict

    @staticmethod
    def parse_samples(sampleFile):
        '''
        Parameters:
            sampleFile -- a string pointing to a .psam (with a '#'-prefixed header)
                          or .fam file
        Returns:
            sampleDict -- a dictionary with structure like:
                          {
                              iid: {"fid": str, "pat": str, "mat": str},
                              ...
                          }
                          where FID defaults to the IID, and parents to '0', when
                          a .psam lacks those columns
        '''
        _, rows = parse_sample_info(sampleFile)
        return { iid: { "fid": row.get("FID", iid), "pat": row.get("PAT", "0"), "mat": row.get("MAT", "0") }
                 for iid, row in rows.items() }

    def parent_offspring_pairs(self):
        '''
        Returns:
            pairs -- a list of [iid, fid, role, parent, kinship, ibs0] for every
                     recorded parent-offspring pair where both are genotyped
        '''
        genotyped = self.genotyped
        pairs = []
        for iid, info in self.sampleDict.items():
            if iid not in genotyped:
                continue
            for role in ["pat", "mat"]:
                parent = info[role]
                if parent == "0" or parent not in genotyped:
                    continue
                kinship, ibs0 = self.kinDict[frozenset((iid, parent))]
                pairs.append([iid, info["fid"], role.upper(), parent, kinship, ibs0])
        return pairs

    def best_candidate(self, iid, fid):
        '''
        Finds the genotyped sample with the highest kinship to iid which is not
        one of its siblings (i.e., not in its recorded family), as a candidate
        true parent.

        Returns:
            candidate -- a list of [candidateIID, kinship, ibs0], or None if
                         there are no candidates
        '''
        siblings = { s for s, info in self.sampleDict.items() if info["fid"] == fid }
        candidates = [ (self.kinDict[frozenset((iid, other))], other)
                       for other in self.genotyped
                       if other != iid and other not in siblings and frozenset((iid, other)) in self.kinDict ]
        if candidates == []:
            return None
        (kinship, ibs0), candidate = max(candidates)
        return [candidate, kinship, ibs0]

    def check_pedigree(self, outputFileName):
        '''
        Writes a per-pair report of recorded parent-offspring relationships, with
        the best candidate parent for pairs below minKinship.

        Parameters:
            outputFileName -- a string of the TSV location to write
        Returns:
            stats -- a dictionary with keys "pairs", "failed", "medianKinship"
                     and "medianIBS0" (the medians are None if there are no pairs)
        '''
        pairs = self.parent_offspring_pairs()
        stats = { "pairs": len(pairs), "failed": 0, "medianKinship": None, "medianIBS0": None }
        if pairs != []:
            stats["medianKinship"] = statistics.median([ p[4] for p in pairs ])
            stats["medianIBS0"] = statistics.median([ p[5] for p in pairs ])

        with open(outputFileName, "w") as fileOut:
            fileOut.write("\t".join(KING.PEDIGREE_CHECK_HEADER) + "\n")
            for iid, fid, role, parent, kinship, ibs0 in pairs:
                status = "PASS" if kinship >= self.minKinship else "FAIL"
                best = ["NA", "NA", "NA"]
                if status == "FAIL":
                    stats["failed"] += 1
                    candidate = self.best_candidate(iid, fid)
                    if candidate is not None:
                        best = [candidate[0], f"{candidate[1]:.4f}", f"{candidate[2]:.4f}"]
                fileOut.write("\t".join([iid, fid, role, parent, f"{kinship:.4f}", f"{ibs0:.4f}", status] + best) + "\n")
        return stats

class VCFToPLINK(PLINKTool):
    '''
    Encapsulates the conversion of a VCF into a PLINK fileset (pgen, bgen or bed) with
    harmonised chromosome names, unique variant IDs, and sample information taken
    from a .fam file matched by sample IID (rather than by row order, which is how
    plink2 --vcf --fam would otherwise assign it).

    Properties are given sensible defaults and can be adjusted prior to calling
    convert(). Intermediate files are written to a temporary directory which is
    removed afterwards (unless keepTmp is True).

    Example:
        converter = VCFToPLINK("input.vcf.gz", "samples.fam", "output/prefix")
        converter.chromMap = "chrom_map.tsv"
        converter.threads = 8
        converter.convert()
    '''
    HALF_CALL_MODES = ["error", "haploid", "missing", "reference"]
    PSAM_HEADER = ["#FID", "IID", "PAT", "MAT", "SEX", "PHENO1"]

    def __init__(self, vcfFile, famFile, outputPrefix):
        super().__init__()
        self.vcfFile = vcfFile
        self.famFile = famFile
        self.outputPrefix = outputPrefix

        # Set default values for adjustable properties
        self.chromMap = None
        self.sampleMap = None
        self.allowMissing = False
        self.maxAlleles = 2
        self.varIdFormat = "@:#:$r:$a"
        self.halfCall = "missing"

    @property
    def vcfFile(self):
        return self._vcfFile

    @vcfFile.setter
    def vcfFile(self, value):
        if not isinstance(value, str):
            raise TypeError("vcfFile must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"VCF file '{value}' is not a file")
        self._vcfFile = os.path.abspath(value)

    @property
    def famFile(self):
        return self._famFile

    @famFile.setter
    def famFile(self, value):
        if not isinstance(value, str):
            raise TypeError("famFile must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"Fam file '{value}' is not a file")
        self._famFile = os.path.abspath(value)

    @property
    def outputPrefix(self):
        return self._outputPrefix

    @outputPrefix.setter
    def outputPrefix(self, value):
        if not isinstance(value, str):
            raise TypeError("outputPrefix must be a string")
        value = os.path.abspath(value)
        if not os.path.isdir(os.path.dirname(value)):
            raise FileNotFoundError(f"Output directory '{os.path.dirname(value)}' does not exist")
        self._outputPrefix = value

    @property
    def chromMap(self):
        return self._chromMap

    @chromMap.setter
    def chromMap(self, value):
        '''
        Accepts None (keep all contigs with their VCF names), a dictionary of
        VCF contig: output chromosome, or a two-column TSV file of the same.
        Contigs not listed in a map are excluded from the output.
        '''
        if value is None:
            self._chromMap = None
        elif isinstance(value, dict):
            self._chromMap = { str(k): str(v) for k, v in value.items() }
        elif isinstance(value, str):
            if not os.path.isfile(value):
                raise FileNotFoundError(f"Chromosome map file '{value}' is not a file")
            self._chromMap = parse_2col_tsv_as_dict(value)
        else:
            raise TypeError("chromMap must be None, a dictionary, or a file path string")
        if self._chromMap == {}:
            raise ValueError("chromMap must contain at least one contig")

    @property
    def sampleMap(self):
        return self._sampleMap

    @sampleMap.setter
    def sampleMap(self, value):
        '''
        Accepts None (VCF sample names are .fam IIDs), a dictionary of VCF sample:
        .fam IID, or a two-column TSV file of the same. VCF samples not listed in
        a map are excluded from the output.
        '''
        if value is None:
            self._sampleMap = None
            return
        elif isinstance(value, dict):
            sampleMap = { str(k): str(v) for k, v in value.items() }
        elif isinstance(value, str):
            if not os.path.isfile(value):
                raise FileNotFoundError(f"Sample map file '{value}' is not a file")
            sampleMap = parse_2col_tsv_as_dict(value)
        else:
            raise TypeError("sampleMap must be None, a dictionary, or a file path string")
        if len(set(sampleMap.values())) != len(sampleMap):
            raise ValueError("sampleMap renames more than one VCF sample to the same IID")
        self._sampleMap = sampleMap

    @property
    def allowMissing(self):
        return self._allowMissing

    @allowMissing.setter
    def allowMissing(self, value):
        '''
        If True, samples absent from the .fam are kept with FID=IID, unknown
        parents/sex and missing phenotype; if False, they raise an error.
        '''
        if not isinstance(value, bool):
            raise TypeError("allowMissing must be a boolean")
        self._allowMissing = value

    @property
    def maxAlleles(self):
        return self._maxAlleles

    @maxAlleles.setter
    def maxAlleles(self, value):
        '''
        Maximum alleles (REF + ALTs) for a variant to be retained; 2 keeps only
        biallelic variants, and None imposes no limit.
        '''
        if value is not None and (not isinstance(value, int) or value < 2):
            raise ValueError("maxAlleles must be None or an integer >= 2")
        self._maxAlleles = value

    @property
    def varIdFormat(self):
        return self._varIdFormat

    @varIdFormat.setter
    def varIdFormat(self, value):
        '''
        Template for plink2 --set-all-var-ids, or None to keep VCF IDs.
        '''
        if value is not None and not isinstance(value, str):
            raise TypeError("varIdFormat must be None or a string")
        self._varIdFormat = value

    @property
    def halfCall(self):
        return self._halfCall

    @halfCall.setter
    def halfCall(self, value):
        '''
        How plink2 handles half-calls (e.g., './1') via --vcf-half-call.
        '''
        if value not in VCFToPLINK.HALF_CALL_MODES:
            raise ValueError(f"halfCall must be one of {VCFToPLINK.HALF_CALL_MODES}")
        self._halfCall = value

    @property
    def outputFiles(self):
        '''
        The files that convert() writes for the current outputFormat and compression.
        '''
        return self.output_files(self.outputPrefix)

    @staticmethod
    def get_max_allele_length(pvarFile):
        '''
        Finds the longest REF or ALT allele in a .pvar file, which is the value
        needed for plink2 --new-id-max-allele-len so that no variant ID is
        truncated or left missing.

        Parameters:
            pvarFile -- a string pointing to a (possibly gzipped) .pvar file
        Returns:
            maxLength -- an integer of the longest allele length
        '''
        maxLength = 0
        with read_gz_file(pvarFile) as fileIn:
            for line in fileIn:
                if line.startswith("#"):
                    continue
                sl = line.split("\t", 5)
                ref, alts = sl[3], sl[4]
                maxLength = max(maxLength, len(ref), *[ len(alt) for alt in alts.split(",") ])
        return maxLength

    @staticmethod
    def parse_fam(famFile):
        '''
        Parses a .fam file into a dictionary keyed by IID.

        Parameters:
            famFile -- a string pointing to a 6-column .fam file
        Returns:
            famDict -- a dictionary with structure like:
                       {
                           iid: [fid, iid, pat, mat, sex, pheno],
                           ...
                       }
        '''
        famDict = {}
        with read_gz_file(famFile) as fileIn:
            for line in fileIn:
                sl = line.split()
                if sl == []:
                    continue
                if len(sl) != 6:
                    raise ValueError(f"Fam file line should have six columns, but has {len(sl)}; offending line is '{line.rstrip()}'")
                if sl[1] in famDict:
                    raise ValueError(f"Duplicate IID '{sl[1]}' in fam file '{famFile}'")
                famDict[sl[1]] = sl
        return famDict

    def import_vcf(self, workDir):
        '''
        Imports the VCF into an uncompressed pgen fileset within workDir, keeping
        variants on mapped contigs with at most maxAlleles alleles, and samples
        listed in the sample map (if one is given).

        Parameters:
            workDir -- a string pointing to the working directory
        Returns:
            rawPrefix -- a string of the imported fileset's prefix
        '''
        rawPrefix = os.path.join(workDir, "raw")
        args = ["--vcf", self.vcfFile, "--vcf-half-call", self.halfCall]
        if self.chromMap is not None:
            args += ["--chr", ",".join(self.chromMap.keys())]
        if self.maxAlleles is not None:
            args += ["--max-alleles", str(self.maxAlleles)]
        if self.sampleMap is not None:
            keepFile = os.path.join(workDir, "keep.txt")
            with open(keepFile, "w") as fileOut:
                fileOut.write("#IID\n" + "\n".join(self.sampleMap.keys()) + "\n")
            args += ["--keep", keepFile]
        self.run_plink(args + ["--make-pgen"], rawPrefix)
        return rawPrefix

    def rename_chromosomes(self, pvarFile):
        '''
        Renames chromosomes in an uncompressed .pvar file according to chromMap.
        plink2 normalises standard names on import (e.g., 'chr1' -> '1'), so map
        keys are matched with and without a leading 'chr'.

        Parameters:
            pvarFile -- a string pointing to an uncompressed .pvar file; it is
                        rewritten in place
        '''
        lookup = {}
        for key, value in self.chromMap.items():
            lookup[key] = value
            lookup[key[3:] if key.startswith("chr") else key] = value

        tmpFile = pvarFile + ".tmp"
        with open(pvarFile, "r") as fileIn, open(tmpFile, "w") as fileOut:
            for line in fileIn:
                if not line.startswith("#"):
                    chrom, rest = line.split("\t", 1)
                    if chrom not in lookup:
                        raise ValueError(f"Chromosome '{chrom}' in '{pvarFile}' is not in the chromosome map")
                    line = f"{lookup[chrom]}\t{rest}"
                fileOut.write(line)
        os.replace(tmpFile, pvarFile)

    def update_psam(self, psamFile):
        '''
        Rewrites a .psam with sample information from the .fam file, matched by
        IID after applying the sample map (if any), preserving the .psam sample
        order so that it remains valid alongside the existing .pgen.

        Parameters:
            psamFile -- a string pointing to the .psam file; it is rewritten in place
        '''
        famDict = VCFToPLINK.parse_fam(self.famFile)

        # Read the IIDs in .psam order
        with open(psamFile, "r") as fileIn:
            header = fileIn.readline().rstrip("\r\n").lstrip("#").split("\t")
            iidIndex = header.index("IID")
            iids = [ line.rstrip("\r\n").split("\t")[iidIndex] for line in fileIn if line.strip() != "" ]

        # Rename and check against the .fam
        if self.sampleMap is not None:
            iids = [ self.sampleMap[iid] for iid in iids ] # import_vcf() kept only mapped samples
        missing = [ iid for iid in iids if iid not in famDict ]
        if missing != [] and not self.allowMissing:
            raise ValueError(f"{len(missing)} samples are absent from the fam file '{self.famFile}', e.g., {missing[:5]}; " +
                             "set allowMissing to keep them as unphenotyped")

        # Write the updated .psam
        tmpFile = psamFile + ".tmp"
        with open(tmpFile, "w") as fileOut:
            fileOut.write("\t".join(VCFToPLINK.PSAM_HEADER) + "\n")
            for iid in iids:
                fid, _, pat, mat, sex, pheno = famDict.get(iid, [iid, iid, "0", "0", "0", "-9"])
                sex = sex if sex in ["1", "2"] else "NA"
                fileOut.write("\t".join([fid, iid, pat, mat, sex, pheno]) + "\n")
        os.replace(tmpFile, psamFile)

    def convert(self):
        '''
        Converts the VCF into the configured PLINK2 output fileset, via a
        temporary working directory, setting variant IDs and sorting variants.

        Returns:
            outputFiles -- a list of the output files written
        '''
        self.check_outputs_absent(self.outputFiles)

        with self.work_dir(self.outputPrefix) as workDir:
            rawPrefix = self.import_vcf(workDir)
            if self.chromMap is not None:
                self.rename_chromosomes(rawPrefix + ".pvar")
            self.update_psam(rawPrefix + ".psam")

            args = ["--pfile", rawPrefix]
            if self.varIdFormat is not None:
                maxAlleleLength = VCFToPLINK.get_max_allele_length(rawPrefix + ".pvar")
                args += ["--set-all-var-ids", self.varIdFormat,
                         "--new-id-max-allele-len", str(max(maxAlleleLength, 1))]
            return self.write_fileset(args, self.outputPrefix, workDir, sortVars=True)

def plink_filter(args):
    '''
    Handles "plink filter" mode of variantopia.
    '''
    plink = PLINK(args.plinkPrefix, args.inputFormat)
    plink.plinkExe = args.plink2
    plink.threads = args.threads
    plink.memory = args.memory
    plink.outputFormat = args.outputFormat
    plink.removeFile = args.removeFile
    plink.requirePheno = args.requirePheno
    plink.nonfounders = args.nonfounders
    plink.snpsOnly = args.snpsOnly
    plink.geno = args.geno
    plink.maf = args.maf

    # Filter and report on what remains
    outputFiles = plink.filter(args.outputPrefix)
    counts = plink.filterCounts
    print(f"# Samples: {counts.get('samples', 'NA')}")
    print(f"# Variants passing filters: {counts.get('variants', 'NA')}")
    print(f"# Wrote {', '.join(outputFiles)}")

def plink_king(args):
    '''
    Handles "plink king" mode of variantopia.
    '''
    plink = PLINK(args.plinkPrefix, args.inputFormat)
    plink.plinkExe = args.plink2
    plink.threads = args.threads
    plink.memory = args.memory

    # Estimate kinship
    kingPrefix = os.path.join(args.outputDirectory, os.path.basename(args.plinkPrefix) + ".king")
    checkFile = kingPrefix + ".pedigree_check.tsv"
    heatmapFile = kingPrefix + ".heatmap.png"
    plink.check_outputs_absent([checkFile, heatmapFile])
    kin0File = plink.king(kingPrefix, args.extractFile)
    print(f"# Wrote {kin0File}")

    with plink.work_dir(kingPrefix) as workDir:
        psamFile = plink.write_psam(os.path.join(workDir, "samples"))

        # Check recorded parent-offspring pairs
        king = KING(kin0File, psamFile)
        king.minKinship = args.minKinship
        stats = king.check_pedigree(checkFile)
        if stats["pairs"] == 0:
            print("# No genotyped parent-offspring pairs are recorded in the sample information")
        else:
            print(f"# {stats['pairs']} genotyped parent-offspring pairs; median kinship = " +
                  f"{stats['medianKinship']:.3f}, median IBS0 = {stats['medianIBS0']:.4f}")
            print(f"# {stats['failed']} parent-offspring pairs fall below kinship {args.minKinship}")
        print(f"# Wrote {checkFile}")

        # Plot clustered kinship with recorded parents
        run_rscript("plot_king.R", ["-k", kin0File, "-f", psamFile, "-c", checkFile,
                                    "-o", kingPrefix, "--title", args.title],
                    args.rscript, logFile=kingPrefix + ".plot.log")
    print(f"# Wrote {heatmapFile}")

def plink_pca(args):
    '''
    Handles "plink pca" mode of variantopia.
    '''
    plink = PLINK(args.plinkPrefix, args.inputFormat)
    plink.plinkExe = args.plink2
    plink.threads = args.threads
    plink.memory = args.memory
    plink.nonfounders = args.nonfounders
    plink.npcs = args.npcs

    # Run PCA
    pcaPrefix = os.path.join(args.outputDirectory, os.path.basename(args.plinkPrefix) + ".pca")
    plotFiles = [ pcaPrefix + suffix for suffix in [".summary.tsv", ".scree.png", ".pcs.png", ".loadings.png"] ]
    plink.check_outputs_absent(plotFiles)
    eigenvecFile, eigenvalFile, loadingsFile = plink.pca(pcaPrefix)
    print(f"# Wrote {', '.join([eigenvecFile, eigenvalFile, loadingsFile])}")

    # Plot PCA diagnostics, using a .psam to obtain family IDs and phenotypes
    with plink.work_dir(pcaPrefix) as workDir:
        psamFile = plink.write_psam(os.path.join(workDir, "samples"))
        run_rscript("plot_pca.R", ["-e", eigenvecFile, "-v", eigenvalFile, "-f", psamFile,
                                   "-l", loadingsFile, "-o", pcaPrefix, "--title", args.title,
                                   "--plot_pcs", str(args.plotPCs)],
                    args.rscript, logFile=pcaPrefix + ".plot.log")
    print(f"# Wrote {', '.join(plotFiles)}")

def plink_prune(args):
    '''
    Handles "plink prune" mode of variantopia.
    '''
    plink = PLINK(args.plinkPrefix, args.inputFormat)
    plink.plinkExe = args.plink2
    plink.threads = args.threads
    plink.memory = args.memory
    plink.tmpDir = args.tmpDir
    plink.keepTmp = args.keepTmp
    plink.outputFormat = args.outputFormat
    plink.removeFile = args.removeFile
    plink.requirePheno = args.requirePheno
    plink.nonfounders = args.nonfounders
    plink.snpsOnly = True # indel calls are noisier
    plink.geno = args.geno
    plink.maf = args.maf
    plink.r2 = args.r2
    plink.windowKb = args.windowKb
    plink.densityCap = args.densityCap

    # Prune and report on each stage
    outputFiles = plink.prune(args.outputPrefix)
    counts = plink.pruneCounts
    print(f"# Samples: {counts['samples']}")
    print(f"# Variants passing filters: {counts['filtered']}")
    print(f"# Stage 1 pruned variants: {counts['stage1']}")
    if "stage2" in counts:
        print(f"# Stage 2 pruned variants: {counts['stage2']}")
    if "densityCapped" in counts:
        print(f"# Density cap: median {counts['densityMedianBin']} variants per bin; cap {counts['densityBinCap']} per bin")
        print(f"# Density-capped variants: {counts['densityCapped']}")
    print(f"# Wrote {', '.join(outputFiles)}")

def plink_qc(args):
    '''
    Handles "plink qc" mode of variantopia.
    '''
    plink = PLINK(args.plinkPrefix, args.inputFormat)
    plink.plinkExe = args.plink2
    plink.threads = args.threads
    plink.memory = args.memory

    # Report on samples
    outputFiles = plink.sample_qc(os.path.join(args.outputDirectory, os.path.basename(args.plinkPrefix)))
    print(f"# Wrote {', '.join(outputFiles)}")
