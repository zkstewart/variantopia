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

import os, sys, shutil, subprocess, tempfile

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, parse_2col_tsv_as_dict

class PLINK:
    '''
    Encapsulates the conversion of a VCF into a PLINK2 fileset (pgen or bgen) with
    harmonised chromosome names, unique variant IDs, and sample information taken
    from a .fam file matched by sample IID (rather than by row order, which is how
    plink2 --vcf --fam would otherwise assign it).

    Properties are given sensible defaults and can be adjusted prior to calling
    convert(). Intermediate files are written to a temporary directory which is
    removed afterwards (unless keepTmp is True).

    Example:
        plink = PLINK("input.vcf.gz", "samples.fam", "output/prefix")
        plink.chromMap = "chrom_map.tsv"
        plink.threads = 8
        plink.convert()
        plink.sample_qc()
    '''
    OUTPUT_FORMATS = ["pgen", "bgen"]
    COMPRESSION_TYPES = { # first value is the default for that format
        "pgen": [None, "zstd"], # zstd compresses the .pvar ('vzs')
        "bgen": ["zlib", "zstd"] # zlib == bgen v1.2, zstd == bgen v1.3
    }
    BGEN_VERSIONS = { "zlib": "1.2", "zstd": "1.3" }
    HALF_CALL_MODES = ["error", "haploid", "missing", "reference"]
    PSAM_HEADER = ["#FID", "IID", "PAT", "MAT", "SEX", "PHENO1"]

    def __init__(self, vcfFile, famFile, outputPrefix):
        self.vcfFile = vcfFile
        self.famFile = famFile
        self.outputPrefix = outputPrefix

        # Set default values for adjustable properties
        self.chromMap = None
        self.sampleMap = None
        self.allowMissing = False
        self.outputFormat = "pgen"
        self.compression = None
        self.bgenBits = 8
        self.maxAlleles = 2
        self.varIdFormat = "@:#:$r:$a"
        self.halfCall = "missing"
        self.threads = 1
        self.memory = 16000
        self.tmpDir = None
        self.keepTmp = False
        self._plinkExe = shutil.which("plink2") # may be None; validated when plink2 is first run

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
    def outputFormat(self):
        return self._outputFormat

    @outputFormat.setter
    def outputFormat(self, value):
        if value not in PLINK.OUTPUT_FORMATS:
            raise ValueError(f"outputFormat must be one of {PLINK.OUTPUT_FORMATS}")
        self._outputFormat = value

        # Reset compression to the format default if it is now incompatible
        if getattr(self, "_compression", "unset") not in PLINK.COMPRESSION_TYPES[value]:
            self._compression = PLINK.COMPRESSION_TYPES[value][0]

    @property
    def compression(self):
        return self._compression

    @compression.setter
    def compression(self, value):
        '''
        For pgen, None or 'zstd' (Zstd-compressed .pvar); for bgen, 'zlib'
        (bgen v1.2) or 'zstd' (bgen v1.3).
        '''
        allowed = PLINK.COMPRESSION_TYPES[self.outputFormat]
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
        if value not in PLINK.HALF_CALL_MODES:
            raise ValueError(f"halfCall must be one of {PLINK.HALF_CALL_MODES}")
        self._halfCall = value

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

    @property
    def outputFiles(self):
        '''
        The files that convert() writes for the current outputFormat and compression.
        '''
        if self.outputFormat == "pgen":
            pvarSuffix = ".pvar.zst" if self.compression == "zstd" else ".pvar"
            suffixes = [".pgen", pvarSuffix, ".psam"]
        else:
            suffixes = [".bgen", ".sample"]
        return [ self.outputPrefix + suffix for suffix in suffixes ]

    @property
    def inputArgs(self):
        '''
        plink2 arguments for reading the output fileset back in, e.g., for sample_qc().
        '''
        if self.outputFormat == "pgen":
            args = ["--pfile", self.outputPrefix]
            if self.compression == "zstd":
                args.append("vzs")
        else:
            args = ["--bgen", self.outputPrefix + ".bgen", "ref-first",
                    "--sample", self.outputPrefix + ".sample"]
        return args

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
        famDict = PLINK.parse_fam(self.famFile)

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
            fileOut.write("\t".join(PLINK.PSAM_HEADER) + "\n")
            for iid in iids:
                fid, _, pat, mat, sex, pheno = famDict.get(iid, [iid, iid, "0", "0", "0", "-9"])
                sex = sex if sex in ["1", "2"] else "NA"
                fileOut.write("\t".join([fid, iid, pat, mat, sex, pheno]) + "\n")
        os.replace(tmpFile, psamFile)

    def write_output(self, rawPrefix, workDir):
        '''
        Writes the final output fileset from the adjusted intermediate fileset,
        setting variant IDs and sorting variants.

        Parameters:
            rawPrefix -- a string of the adjusted intermediate fileset's prefix
            workDir -- a string pointing to the working directory
        '''
        args = ["--pfile", rawPrefix]
        if self.varIdFormat is not None:
            maxAlleleLength = PLINK.get_max_allele_length(rawPrefix + ".pvar")
            args += ["--set-all-var-ids", self.varIdFormat,
                     "--new-id-max-allele-len", str(max(maxAlleleLength, 1))]

        if self.outputFormat == "pgen":
            makeArgs = ["--sort-vars", "--make-pgen"] + (["vzs"] if self.compression == "zstd" else [])
            self.run_plink(args + makeArgs, self.outputPrefix)
        else:
            # --sort-vars requires --make-pgen, so sort into the working directory first
            sortedPrefix = os.path.join(workDir, "sorted")
            self.run_plink(args + ["--sort-vars", "--make-pgen"], sortedPrefix)
            self.run_plink(["--pfile", sortedPrefix,
                            "--export", f"bgen-{PLINK.BGEN_VERSIONS[self.compression]}",
                            f"bits={self.bgenBits}", "ref-first"], self.outputPrefix)

    def convert(self):
        '''
        Converts the VCF into the configured PLINK2 output fileset, via a
        temporary working directory.

        Returns:
            outputFiles -- a list of the output files written
        '''
        existing = [ f for f in self.outputFiles if os.path.exists(f) ]
        if existing != []:
            raise FileExistsError(f"Output file(s) already exist: {existing}")

        workDir = tempfile.mkdtemp(prefix=".variantopia_plink_",
                                   dir=self.tmpDir or os.path.dirname(self.outputPrefix))
        try:
            rawPrefix = self.import_vcf(workDir)
            if self.chromMap is not None:
                self.rename_chromosomes(rawPrefix + ".pvar")
            self.update_psam(rawPrefix + ".psam")
            self.write_output(rawPrefix, workDir)
        finally:
            if not self.keepTmp:
                shutil.rmtree(workDir, ignore_errors=True)

        return self.outputFiles

    def sample_qc(self, outputFileName=None):
        '''
        Writes a per-sample QC report of missingness (all variants) and
        heterozygosity rate (SNPs only, since indel calls are noisier) for the
        output fileset produced by convert().

        Parameters:
            outputFileName -- (OPTIONAL) a string of the report location; default
                              is '<outputPrefix>.sample_qc.tsv'
        Returns:
            outputFileName -- a string of the report location
        '''
        if outputFileName is None:
            outputFileName = self.outputPrefix + ".sample_qc.tsv"
        missing = [ f for f in self.outputFiles if not os.path.isfile(f) ]
        if missing != []:
            raise FileNotFoundError(f"Output file(s) {missing} not found; run convert() first")

        qcPrefix = self.outputPrefix + ".sample_qc"
        self.run_plink(self.inputArgs + ["--missing", "sample-only"], qcPrefix)
        self.run_plink(self.inputArgs + ["--snps-only", "just-acgt", "--sample-counts"], qcPrefix)

        # Combine the plink2 reports
        def parse_report(fileName):
            with open(fileName, "r") as fileIn:
                header = fileIn.readline().rstrip("\r\n").lstrip("#").split("\t")
                return [ dict(zip(header, line.rstrip("\r\n").split("\t"))) for line in fileIn ]
        smiss = parse_report(qcPrefix + ".smiss")
        scount = { row["IID"]: row for row in parse_report(qcPrefix + ".scount") }

        with open(outputFileName, "w") as fileOut:
            fileOut.write("FID\tIID\tF_MISS\tHET_RATE\n")
            for row in smiss:
                counts = scount[row["IID"]]
                called = int(counts["HOM_REF_CT"]) + int(counts["HOM_ALT_SNP_CT"]) + int(counts["HET_SNP_CT"])
                hetRate = f"{int(counts['HET_SNP_CT']) / called:.6g}" if called > 0 else "NA"
                fileOut.write(f"{row.get('FID', row['IID'])}\t{row['IID']}\t{row['F_MISS']}\t{hetRate}\n")

        return outputFileName
