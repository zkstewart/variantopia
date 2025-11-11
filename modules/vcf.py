#! python3
# vcf.py
# Classes and functions for parsing VCF files into suitable
# data structures for variantopia.

import os, shutil, subprocess, warnings
import numpy as np
import pandas as pd
from cyvcf2 import VCF
from collections import Counter
from copy import deepcopy

class VCFTopia:
    def __init__(self, vcfFile):
        self.vcfFile = vcfFile
        self.vcf = VCF(vcfFile)
        self.isVCFTopia = True
    
    @staticmethod
    def count_alleles(variant):
        '''
        Takes a cyvcf2.Variant object and returns a sorted list of tuples
        containing each allele (as an integer code) and its count.
        
        Parameters:
            variant -- a cyvcf2.Variant object
        Returns:
            alleleCounts -- a list of tuples where each tuple contains an allele
                            and its count, sorted by count in descending order.
        '''
        alleles = [ allele for genotype in variant.genotypes for allele in genotype[:-1]]
        counts = Counter(alleles)
        alleleCounts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        
        return alleleCounts
    
    @staticmethod
    def count_homhet(variant, countAbsence=False):
        '''
        Takes a cyvcf2.Variant object and returns a Counter object with
        two keys: "hom" and "het", where "hom" contains the count of samples with
        homozygous genotype, and "het" contains the count of samples with
        heterozygous genotype.
        
        Parameters:
            variant -- a cyvcf2.Variant object
            countAbsence -- (OPTIONAL) a boolean controlling whether we should
                            count missing data as being homozygous (True) or
                            skip the datapoint (False); default is False
        Returns:
            homHet -- a collections.Counter object with structure like:
                      {
                          "hom": count_of_homozygous_samples,
                          "het": count_of_heterozygous_samples
                      }
        '''
        genotypeCategories = [
            "hom" if len(set(genotype[:-1])) == 1 else "het"
            for genotype in variant.genotypes
            if countAbsence or (not -1 in genotype[:-1])
        ]
        homHet = Counter(genotypeCategories)
        homHet.setdefault("hom", 0)
        homHet.setdefault("het", 0)
        
        return homHet
    
    @staticmethod
    def calculate_mac(variant):
        '''
        Takes a cyvcf2.Variant object and returns the minor allele count (MAC),
        which is the count of the second most common allele.
        
        Parameters:
            variant -- a cyvcf2.Variant object
        Returns:
            mac -- the minor allele count (MAC) for the variant.
        '''
        alleleCounts = VCFTopia.count_alleles(variant)
        #majorAlleles = [ allele for allele, count in orderedAlleles if count == orderedAlleles[0][1] ]
        #return sum(count for allele, count in alleleCounts if not allele in majorAlleles)
        
        for allele, count in alleleCounts[1:]:
            if allele != -1:
                return count
        return 0
    
    @staticmethod
    def calculate_maf(variant):
        '''
        Takes a cyvcf2.Variant object and returns the minor allele frequency (MAF),
        which is the frequency of the second most common allele.
        
        Parameters:
            variant -- a cyvcf2.Variant object
        Returns:
            maf -- the minor allele frequency (MAF) for the variant.
        '''
        alleleCounts = VCFTopia.count_alleles(variant)
        #majorAlleles = [ allele for allele, count in orderedAlleles if count == orderedAlleles[0][1] ]
        #return sum(count for allele, count in alleleCounts if not allele in majorAlleles) / sum(count for allele, count in alleleCounts)
        
        for allele, count in alleleCounts[1:]:
            if allele != -1:
                return count / sum(count for allele, count in alleleCounts if allele != -1)
        return 0.0
    
    @staticmethod
    def calculate_callrate(variant):
        '''
        Takes a cyvcf2.Variant object and returns the callrate, which is the
        proportion of samples that have a non-missing genotype for the variant.
        
        Parameters:
            variant -- a cyvcf2.Variant object
        Returns:
            callrate -- the callrate for the variant, as a float between 0 and 1.
        '''
        alleleCounts = VCFTopia.count_alleles(variant)
        
        for allele, count in alleleCounts[1:]:
            if allele == -1:
                return 1 - (count / sum(count for allele, count in alleleCounts))
        return 1.0
    
    @staticmethod
    def calculate_heterozygosity(variant):
        '''
        Takes a cyvcf2.Variant object and returns the heterozygosity, which is the
        proportion of samples that have a heterozygous genotype call for the variant.
        
        Parameters:
            variant -- a cyvcf2.Variant object
        Returns:
            heterozygosity -- the sample heterozygosity for the variant, as a float
                              between 0 and 1.
        '''
        homHet = VCFTopia.count_homhet(variant)
        return homHet["het"] / (homHet["hom"] + homHet["het"]) if (homHet["hom"] + homHet["het"]) > 0 else 0.0
    
    @staticmethod
    def refAlt(variant):
        '''
        Takes a cyvcf2.Variant object and returns a list
        representing the reference and alternate alleles
        in order of their VCF integer value.
        
        Parameters:
            variant -- a cyvcf2.Variant object
        
        Returns:
            refAlt -- a list of strings representing the reference (index 0)
                      and alternate alleles (index 1 and onwards) e.g., ["A", "C", ...]
        '''
        return [variant.REF, *variant.ALT]
    
    @staticmethod
    def genotype_codes(variant):
        '''
        Takes a cyvcf2.Variant object and returns a string
        representing the genotype in VCF-encoded format.
        
        Parameters:
            variant -- a cyvcf2.Variant object
        
        Returns:
            genotypeCodes -- a list of strings representing the genotype codes
                             e.g., ["0/0", "0/1", "1/1", ...]
        '''
        return [
            "/".join(str(code if code != -1 else ".") for code in gt[:-1]) if gt[-1] == False else
            "|".join(str(code if code != -1 else ".") for code in gt[:-1])
            for gt in variant.genotypes
        ]
    
    @staticmethod
    def sample_codes(samples, variant):
        '''
        Takes a list of sample names (from a VCFTopia.samples property) alongside
        a cyvcf2.Variant object and returns a dictionary where keys are sample names
        and values are strings representing the genotype codes as formatted by
        VCFTopia.genotype_codes.
        '''
        return { sample: code for sample, code in zip(samples, VCFTopia.genotype_codes(variant)) }
    
    @staticmethod
    def genotype_alleles(variant):
        '''
        Takes a cyvcf2.Variant object and returns a string
        representing the genotype by substituting VCF-encoded
        genotype codes (i.e., integers) with the actual allele
        (i.e., 'A' / 'T' / 'C' / 'G').
        
        Parameters:
            variant -- a cyvcf2.Variant object
        
        Returns:
            genotypeAlleles -- a list of strings representing the genotypes as
                               alleles e.g., ["A/A", "A/T", "T/T", ...]
        '''
        refAlt = VCFTopia.refAlt(variant)
        return [
            "/".join(refAlt[code] if code != -1 else "." for code in gt[:-1]) if gt[-1] == False else
            "|".join(refAlt[code] if code != -1 else "." for code in gt[:-1])
            for gt in variant.genotypes
        ]
    
    @staticmethod
    def sample_alleles(samples, variant):
        '''
        Takes a list of sample names (from a VCFTopia.samples property) alongside
        a cyvcf2.Variant object and returns a dictionary where keys are sample names
        and values are strings representing the genotype alleles as formatted by
        VCFTopia.genotype_alleles.
        '''
        return { sample: alleles for sample, alleles in zip(samples, VCFTopia.genotype_alleles(variant)) }
    
    @staticmethod
    def sample_integers(samples, variant):
        '''
        Takes a list of sample names (from a VCFTopia.samples property) alongside
        a cyvcf2.Variant object and returns a dictionary where keys are sample names
        and values are the integer genotype lists from variant.genotypes, omitting
        the final phase boolean
        '''
        return { sample: integerList[:-1] for sample, integerList in zip(samples, variant.genotypes) }
    
    @staticmethod
    def sample_format_field(samples, variant, field):
        '''
        Takes a list of sample names (from a VCFTopia.samples property) alongside
        a cyvcf2.Variant object and a value that is found in variant.FORMAT
        and returns a dictionary where keys are sample names
        and values are the matching FORMAT value.
        '''
        return { sample: alleles for sample, alleles in zip(samples, variant.format(field)) }
    
    @property
    def vcfFile(self):
        return self._vcfFile
    
    @vcfFile.setter
    def vcfFile(self, value):
        if not os.path.exists(value):
            raise FileNotFoundError(f"VCF file not found at '{value}'")
        self._vcfFile = value
    
    @property
    def samples(self):
        return self.vcf.samples
    
    @property
    def contigs(self):
        '''
        Accessor for the chromosome names in the VCF file.
        '''
        return self.vcf.seqnames
    
    @property
    def chroms(self):
        '''
        Accessor for the chromosome names in the VCF file.
        '''
        return self.vcf.seqnames
    
    @property
    def lengths(self):
        if hasattr(self.vcf, "seqlens"):
            return { chrom:seqLength for chrom, seqLength in zip(self.vcf.seqnames, self.vcf.seqlens) }
        return None
    
    @property
    def is_bgzipped(self):
        '''
        Check if the VCF file is bgzipped.
        
        Returns:
            True if the VCF file is bgzipped, False otherwise.
        '''
        return self.vcfFile.endswith(".vcf.gz") or self.vcfFile.endswith(".vcf.bgz")
    
    @property
    def is_indexed(self):
        '''
        Check if the VCF file is indexed.
        
        Returns:
            True if the VCF file is indexed, False otherwise.
        '''
        indexFile = self.vcfFile + ".tbi" if self.is_bgzipped else self.vcfFile + ".csi"
        return os.path.exists(indexFile)
    
    def index(self, tabixExe=None, bcftoolsExe=None):
        '''
        Create an index for the VCF to allow for fast querying. If one of either the tabix
        or bcftools executables is provided, it will be used to create the index. If neither
        is provided, VCFTopia will attempt to find tabix in the system PATH and use it, with
        bcftools as a fallback. If this fails, an error will be raised.
        
        Parameters:
            tabixExe -- path to the tabix executable; if None, uses the default in PATH
            bcftoolsExe -- path to the bcftools executable; if None, uses the default in PATH
        '''
        if not self.is_bgzipped:
            raise ValueError(f"Cannot index VCF file '{self.vcfFile}' since it is not bgzipped. " +
                             "Please bgzip the file before indexing it.")
        
        if tabixExe != None:
            self.tabix(tabixExe)
        elif bcftoolsExe != None:
            self.bcftools_index(bcftoolsExe)
        else:
            # Try to find tabix in the system PATH
            tabixPath = shutil.which("tabix")
            if tabixPath != None:
                self.tabix(tabixPath)
            else:
                # If tabix is not found, try to find bcftools
                bcftoolsPath = shutil.which("bcftools")
                if bcftoolsPath != None:
                    self.bcftools_index(bcftoolsPath)
                else:
                    raise FileNotFoundError(f"Neither tabix nor bcftools found in PATH thus cannot index '{self.vcfFile}'; " +
                                            "please index this file manually before using it in variantopia.")
    
    def tabix(self, tabixExe):
        '''
        Create an index for the VCF file using the tabix executable.
        
        Parameters:
            tabixExe -- path to the tabix executable
        '''
        if not self.is_bgzipped:
            raise ValueError(f"Cannot index VCF file '{self.vcfFile}' since it is not bgzipped. " +
                             "Please bgzip the file before using it with variantopia.")
        if not os.path.exists(tabixExe):
            raise FileNotFoundError(f"tabix executable not found at '{tabixExe}'")
        
        # Run tabix to index the VCF file
        cmd = [tabixExe, "-f", "-p", "vcf", self.vcfFile]
        run_tabix = subprocess.Popen(" ".join(cmd), shell = True,
                                    stdout = subprocess.PIPE,
                                    stderr = subprocess.PIPE)
        tabixout, tabixerr = run_tabix.communicate()
        
        # Check that tabix ran successfully
        if tabixout.decode("utf-8") != "":
            raise Exception(("tabix stdout is not empty as expected, indicating a probable error. " +
                            f'Please check the stdout ({tabixout.decode("utf-8")}) and stderr ' + 
                            f'({tabixerr.decode("utf-8")}) to make sense of this.'))
        elif tabixerr.decode("utf-8") != "":
            raise Exception(("tabix encountered an error; have a look " +
                            f'at the stdout ({tabixout.decode("utf-8")}) and stderr ' + 
                            f'({tabixerr.decode("utf-8")}) to make sense of this.'))
        print(f"# Indexed '{self.vcfFile}' with tabix")
    
    def bcftools_index(self, bcftoolsExe):
        '''
        Create an index for the VCF file using the bcftools executable.
        
        Parameters:
            bcftoolsExe -- path to the tabix executable
        '''
        #BAD_WORDS = ["failed", "error", "warning", "abort", "exception", "fatal", "fail", "unrecogni"]
        if not self.is_bgzipped:
            raise ValueError(f"Cannot index VCF file '{self.vcfFile}' since it is not bgzipped. " +
                             "Please bgzip the file before using it with variantopia.")
        if not os.path.exists(bcftoolsExe):
            raise FileNotFoundError(f"bcftools executable not found at '{bcftoolsExe}'")
        
        # Run tabix to index the VCF file
        cmd = [bcftoolsExe, "index", self.vcfFile]
        run_index = subprocess.Popen(" ".join(cmd), shell = True,
                                    stdout = subprocess.PIPE,
                                    stderr = subprocess.PIPE)
        indexout, indexerr = run_index.communicate()
        
        # Check that bcftools index ran successfully
        if indexout.decode("utf-8") != "":
            raise Exception(("bcftools index stdout is not empty as expected, indicating a probable error. " +
                            f'Please check the stdout ({indexout.decode("utf-8")}) and stderr ' + 
                            f'({indexerr.decode("utf-8")}) to make sense of this.'))
        elif indexerr.decode("utf-8") != "":
            raise Exception(("bcftools index encountered an error; have a look " +
                            f'at the stdout ({indexout.decode("utf-8")}) and stderr ' + 
                            f'({indexerr.decode("utf-8")}) to make sense of this.'))
        print(f"# Indexed '{self.vcfFile}' with bcftools index")
    
    def query(self, chrom, startEnd=None):
        '''
        Query the VCF file for variants in a specific region.
        
        Parameters:
            chrom -- a string representing the chromosome or scaffold name
            startEnd -- a tuple (start, end) representing the range to query.
                        If None, queries the entire chromosome.
        '''
        queryRange = f"{chrom}:{startEnd[0]}-{startEnd[1]}" if startEnd else chrom
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning) # ignore warnings from cyvcf2 about 'no intervals found'
            try:
                yield from self.vcf(queryRange)
            except AssertionError as e:
                if "error loading tabix index" in str(e):
                    print("Attempting to index the VCF file due to an indexing error...")
                    self.index()  # Attempt to index the VCF file if the error is related to indexing
                    yield from self.vcf(queryRange)
                else:
                    raise e # re-raise the exception if it's not related to indexing
    
    def as_genotype_dict(self, multialleles=True, indels=True, mnps=True):
        '''
        Returns a data structure commonly utilised in pre-variantopia methods
        which is still of use in some cases.
        
        Note that this method does not account for variants at the same chrom+position
        location which are presented on different lines in the VCF file.
        
        Parameters:
            multialleles -- if True, include multiallelic variants; otherwise, filter them out
            indels -- if True, include indels; otherwise, filter them out
            mnps -- if True, include multi-nucleotide polymorphisms (MNPs); otherwise, filter
                    to only retain single-nucleotide polymorphisms (SNPs)
        Returns:
            genotypeDict -- a dictionary where keys are chromosome names and values
                            are dictionaries with positions as keys and lists of
                            [refAlt, genotypeCodes] as values.
                refAlt -- a list of strings representing the reference (index 0) and alternate
                          alleles (index 1 and onwards) e.g., ["A", "C"] for a biallelic variant,
                          with additional alleles for multiallelic variants ordered according
                          to the VCF integer value.
                genotypeCodes -- a list of strings representing the genotype codes
                                 in VCF-encoded format, e.g., ["./.", "0/0", "0/1", "1|1", ...]
        '''
        genotypeDict = {}
        for variant in self:
            refAlt = VCFTopia.refAlt(variant)
            
            # Skip variants if applicable
            if (not multialleles) and len(variant.ALT) > 1: # skip multiallelic
                continue
            if (not indels) and variant.is_indel: # skip indels
                continue
            if (not mnps) and variant.is_mnp: # skip non-SNPs
                continue
            
            # Store genotype information
            genotypeDict.setdefault(variant.CHROM, {})
            genotypeDict[variant.CHROM][variant.POS] = [refAlt, VCFTopia.genotype_codes(variant)]
        return genotypeDict
    
    def as_dict(self, multialleles=True, indels=True, mnps=True):
        '''
        Returns a dictionary-based data structure which encapsulates the most frequently
        used information from a VCF file in a format that is easy to iterate through
        and work with.
        
        Parameters:
            multialleles -- if True, include multiallelic variants; otherwise, filter them out
            indels -- if True, include indels; otherwise, filter them out
            mnps -- if True, include multi-nucleotide polymorphisms (MNPs); otherwise, filter
                    to only retain single-nucleotide polymorphisms (SNPs)
        Returns:
            vcfDict -- a dictionary where keys are variant identifiers (chrom:pos:ref:alt)
                       and values are dictionaries with keys "chrom", "pos", "ref", "alt",
                       and "sampleData" (a dictionary of sample names to their genotypes)
        '''
        vcfDict = {}
        for variant in self:
            # Skip variants if applicable
            if (not multialleles) and len(variant.ALT) > 1: # skip multiallelic
                continue
            if (not indels) and variant.is_indel: # skip indels
                continue
            if (not mnps) and variant.is_mnp: # skip non-SNPs
                continue
            
            # Store the sample data in the dictionary
            ref = variant.REF
            alt = ",".join(variant.ALT)
            
            key = f"{chrom}:{pos}:{ref}:{alt}" # ensures unique key
            vcfDict[key] = {
                "chrom": variant.CHROM,
                "pos": variant.POS,
                "ref": ref,
                "alt": alt,
                "sampleData": VCFTopia.sample_alleles(variant)
            }
        return vcfDict
    
    def comprehensive_statistics(self):
        '''
        Calculates and tallies all potentially relevant statistics on the
        composition of this VCF file.
        
        Note:
        1) insertion and deletion are not mutually incompatible for a single
        variant site in cases where multiallelic variants exist
        2) multiallelic and biallelic variants are mutuallly incompatible;
        a variant is either biallelic or multiallelic
        3) the final step of adding a row called "Total" assumes the genome
        does NOT have an existing contig called "Total"; it would be silly for this
        to ever be true, but if it is, maybe name your genome contigs something better?
        
        Returns:
            genomeStatsDF -- a pd.Dataframe where rows correspond to individual contigs
                             (with a final Total row at the bottom) and columns correspond
                             to integer counts of different variant types
            sampleStatsDF -- a pd.Dataframe where rows correspond to individual samples
                             and columns correspond to integer or float values of different
                             variant statistics
        '''
        contigDefault = {"num_variants": 0, "num_multiallelic": 0, "num_biallelic": 0,
                         "num_insertions": 0, "num_deletions": 0,
                         "num_snps": 0, "num_mnps": 0}
        sampleDefault = {"num_called": 0, "num_uncalled": 0,
                         "num_hom": 0, "num_het": 0,
                         "DP": [], "AD": []} # lists require dict to be deepcopied later
        
        contigStats = {}
        sampleStats = {}
        for variant in self:
            # Store contig-level stats
            contigStats.setdefault(variant.CHROM, deepcopy(contigDefault))
            contigStats[variant.CHROM]["num_variants"] += 1
            
            if len(variant.ALT) == 1:
                contigStats[variant.CHROM]["num_biallelic"] += 1
            else:
                contigStats[variant.CHROM]["num_multiallelic"] += 1
            
            if any([ len(alt) > len(variant.REF) for alt in variant.ALT ]):
                contigStats[variant.CHROM]["num_insertions"] += 1
            if any([ len(alt) < len(variant.REF) for alt in variant.ALT ]):
                contigStats[variant.CHROM]["num_deletions"] += 1
            
            if len(variant.REF) == 1:
                contigStats[variant.CHROM]["num_snps"] += 1
            else:
                contigStats[variant.CHROM]["num_mnps"] += 1
            
            # Store sample-level stats
            hasDP = False
            if "DP" in variant.FORMAT:
                hasDP = True
                sampleDP = VCFTopia.sample_format_field(self.samples, variant, "DP")
            hasAD = False
            if "AD" in variant.FORMAT:
                hasAD = True
                sampleAD = VCFTopia.sample_format_field(self.samples, variant, "AD")
            
            sampleIntegers = VCFTopia.sample_integers(self.samples, variant)
            for sample, integerList in sampleIntegers.items():
                sampleStats.setdefault(sample, deepcopy(sampleDefault))
                
                if not -1 in integerList:
                    sampleStats[sample]["num_called"] += 1
                    if len(set(integerList)) == 1:
                        sampleStats[sample]["num_hom"] += 1
                    else:
                        sampleStats[sample]["num_het"] += 1
                    
                    if hasDP:
                        sampleStats[sample]["DP"].append(sum(sampleDP[sample])) # needs testing
                    if hasAD:
                        sampleStats[sample]["AD"].append(sum(sampleAD[sample]))
                else:
                    sampleStats[sample]["num_uncalled"] += 1
        
        # Average sample-level statistics arrays
        for sample, statsDict in sampleStats.items():
            statsDict["DP"] = np.array(statsDict["DP"])
            statsDict["AD"] = np.array(statsDict["AD"])
            
            # Median and mean for DP
            if len(statsDict["DP"]) != 0:
                statsDict["DP_median"] = np.median(statsDict["DP"])
                statsDict["DP_mean"] = np.mean(statsDict["DP"])
            else:
                statsDict["DP_median"] = None
                statsDict["DP_mean"] = None
            
            # Median and mean for AD
            if len(statsDict["AD"]) != 0:
                statsDict["AD_median"] = np.median(statsDict["AD"])
                statsDict["AD_mean"] = np.mean(statsDict["AD"])
            else:
                statsDict["AD_median"] = None
                statsDict["AD_mean"] = None
            
            # Calculate DP-AD value (useful for ploidy change inference)
            if len(statsDict["DP"]) != 0 and len(statsDict["AD"]) != 0:
                if len(statsDict["DP"]) == len(statsDict["AD"]):
                    dpMinusAd = statsDict["DP"] - statsDict["AD"]
                    statsDict["DP-AD_median"] = np.median(dpMinusAd)
                    statsDict["DP-AD_mean"] = np.mean(dpMinusAd)
                else:
                    print("WARNING: Inconsistency between count of 'DP' and 'AD' values; cannot calculate DP-AD statistic")
                    statsDict["DP-AD_median"] = None
                    statsDict["DP-AD_mean"] = None
            else:
                statsDict["DP-AD_median"] = None
                statsDict["DP-AD_mean"] = None
            
            # Wipe AP and DP before later pandas dataframe conversion
            del statsDict["DP"]
            del statsDict["AD"]
        
        # Format pandas dataframes for statistics
        genomeStatsDF = pd.DataFrame.from_dict(contigStats, orient="index")
        sampleStatDF = pd.DataFrame.from_dict(sampleStats, orient="index")
        
        # Add totals row for genome statistics
        genomeStatsDF.loc["Total"] = genomeStatsDF.sum()
        
        return genomeStatsDF, sampleStatDF
    
    def __iter__(self):
        """
        Iterate over all variants in the VCF file.
        """
        try:
            yield from self.vcf
            self.vcf = VCF(self.vcfFile) # allow iteration after successful completion
        except Exception as e:
            if str(e) == "attempt to iterate over closed/invalid VCF":
                self.vcf = VCF(self.vcfFile) # prevent "attempt to iterate over closed/invalid VCF"
                yield from self.vcf
            else:
                raise e
    
    def __repr__(self):
        return "<VCFTopia object;file='{0}';num_samples={1};num_contigs={2}>".format(
            self.vcfFile,
            len(self.vcf.samples),
            len(self.chroms)
        )

def vcf_stats(args):
    '''
    Handles "reformat geno" mode of variantopia.
    '''
    vcf = VCFTopia(args.vcfFile)
    genomeStatsDF, sampleStatDF = vcf.comprehensive_statistics()
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("# genome-level statistics\n")
        genomeStatsDF.to_csv(fileOut, sep="\t", na_rep=".")
        fileOut.write("##\n# sample-level statistics\n")
        sampleStatDF.to_csv(fileOut, sep="\t", na_rep=".")
