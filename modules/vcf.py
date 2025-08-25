#! python3
# vcf.py
# Classes and functions for parsing VCF files into suitable
# data structures for variantopia.

import os, shutil, subprocess, warnings
from cyvcf2 import VCF
from collections import Counter

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
    def count_homhet(variant):
        '''
        Takes a cyvcf2.Variant object and returns a Counter object with
        two keys: "hom" and "het", where "hom" contains the count of samples with
        homozygous genotype, and "het" contains the count of samples with
        heterozygous genotype.
        
        Parameters:
            variant -- a cyvcf2.Variant object
        Returns:
            homHet -- a collections.Counter object with structure like:
                      {
                          "hom": count_of_homozygous_samples,
                          "het": count_of_heterozygous_samples
                      }
        '''
        genotypeCategories = [
            "hom" if len(set(genotype[:-1])) == 1 else "het"
            for genotype in variant.genotypes if not -1 in genotype[:-1]
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
    def format_refAlt(variant):
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
    def format_genotype_codes(variant):
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
    def format_sample_codes(samples, variant):
        '''
        Takes a list of sample names (from a VCFTopia.samples property) alongside
        a cyvcf2.Variant object and returns a dictionary where keys are sample names
        and values are strings representing the genotype codes as formatted by
        VCFTopia.format_genotype_codes.
        '''
        return { sample: code for sample, code in zip(samples, VCFTopia.format_genotype_codes(variant)) }
    
    @staticmethod
    def format_genotype_alleles(variant):
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
        refAlt = VCFTopia.format_refAlt(variant)
        return [
            "/".join(refAlt[code] if code != -1 else "." for code in gt[:-1]) if gt[-1] == False else
            "|".join(refAlt[code] if code != -1 else "." for code in gt[:-1])
            for gt in variant.genotypes
        ]
    
    @staticmethod
    def format_sample_alleles(samples, variant):
        '''
        Takes a list of sample names (from a VCFTopia.samples property) alongside
        a cyvcf2.Variant object and returns a dictionary where keys are sample names
        and values are strings representing the genotype alleles as formatted by
        VCFTopia.format_genotype_alleles.
        '''
        return { sample: alleles for sample, alleles in zip(samples, VCFTopia.format_genotype_alleles(variant)) }
    
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
            refAlt = VCFTopia.format_refAlt(variant)
            
            # Skip variants if applicable
            if (not multialleles) and len(variant.ALT) > 1: # skip multiallelic
                continue
            if (not indels) and variant.is_indel: # skip indels
                continue
            if (not mnps) and variant.is_mnp: # skip non-SNPs
                continue
            
            # Store genotype information
            genotypeDict.setdefault(variant.CHROM, {})
            genotypeDict[variant.CHROM][variant.POS] = [refAlt, VCFTopia.format_genotype_codes(variant)]
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
                "sampleData": VCFTopia.format_sample_alleles(variant)
            }
        return vcfDict
    
    def __iter__(self):
        """
        Iterate over all variants in the VCF file.
        """
        yield from self.vcf
    
    def __repr__(self):
        return "<VCFTopia object;file='{0}';num_samples={1};num_contigs={2}>".format(
            self.vcfFile,
            len(vcf.vcf.samples),
            len(vcf.chroms)
        )

if __name__ == "__main__":
    pass
