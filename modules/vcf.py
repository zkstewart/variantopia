#! python3
# vcf.py
# Classes and functions for parsing VCF files into suitable
# data structures for variantopia.

import os
from cyvcf2 import VCF

class VCFTopia:
    def __init__(self, vcfFile):
        self.vcfFile = vcfFile
        self.vcf = VCF(vcfFile)
    
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
    
    def query(self, chrom, startEnd=None):
        '''
        Query the VCF file for variants in a specific region.
        
        Parameters:
            chrom -- a string representing the chromosome or scaffold name
            startEnd -- a tuple (start, end) representing the range to query.
                        If None, queries the entire chromosome.
        '''
        queryRange = f"{chrom}:{startEnd[0]}-{startEnd[1]}" if startEnd else chrom
        yield from self.vcf(queryRange)
    
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

if __name__ == "__main__":
    pass
