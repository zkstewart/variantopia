import os, sys
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from collections import Counter

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from vcf import VCFTopia

def WMA(s, period):
    """
    See https://stackoverflow.com/questions/74518386/improving-weighted-moving-average-performance
    Code copy-pasted out of psQTL.
    
    Parameters:
        s -- a numpy array of values to smooth
        period -- an integer value indicating the number of previous values to consider
                  during weighted moving average calculation
    Returns:
        sw -- a pandas Series of the smoothed values
    """
    w = np.arange(period)+1
    w_s = w.sum()
    
    try:
        swv = np.lib.stride_tricks.sliding_window_view(s.flatten(), window_shape=period)
    except ValueError:
        "Less data points than period size causes this error"
        return None
    sw = (swv * w).sum(axis=1) / w_s
    
    # Need to now return it as a normal series
    sw = np.concatenate((np.full(period - 1, np.nan), sw))
    try:
        sw[0:period] = sw[period] # set first n=period values to be same as first smoothed value
    except:
        "len(sw)==1 causes this error"
        return None
    return sw

def possible_genotypes(gt1, gt2):
    '''
    Calculate all possible inheritable genotypes from two parent genotypes.
    Code copy-pasted out of psQTL.
    
    Parameters:
        gt1 / gt2 -- a list of integers representing the genotype values
                     with format like:
                     [0, 1, 2, ...]
    Returns:
        possibleGTs -- a set of frozensets containing all possible genotypes
                       that can be formed from the two genotypes; sets
                       mask duplicated alleles but allow for testing
                       of genotype without regard to order
    '''
    numGt1 = int(len(gt1) / 2)
    numGt2 = int(len(gt2) / 2)
    
    possibleGTs = []
    for g1 in combinations(gt1, numGt1):
        for g2 in combinations(gt2, numGt2):
            possibleGTs.append(frozenset(g1 + g2))
    
    return set(possibleGTs)

def as_pangenome_variant(variant, referenceIndex):
    '''
    Receives a standard cyvcf2.cyvcf2.Variant object and modifies its genotypes
    field to correctly indicate the reference sample's genotype.
    
    Note that this is only capable of handling diploids.
    
    Parameters:
        variant -- a cyvcf2.cyvcf2.Variant object
        referenceIndex -- the index of the reference sample as found within variant.genotypes
    '''
    refGT = variant.genotypes[referenceIndex]
    if set(refGT[:-1]) != {-1}: # if this site was called
        variant.genotypes[referenceIndex] = [0, *refGT] # 0 == REF code

def assign_sampleGT_to_parent_haplotypes(sampleGT, p1GT, p2GT, PARENT_FULLY_ASSIGNED):
    '''
    Code copy-pasted out of psQTL.
    
    Parameters:
        sampleGT -- a list containing two integers representing the ref (0) or alt (1) VCF codes
        p1GT / p2GT -- same format as sampleGT, but for the genotypes of the two parents
        PARENT_FULLY_ASSIGNED -- a dict with structure like:
                                 {
                                     "p1": 1,
                                     "p2": 2
                                 }
                                 Where the integer value indicates how many times each haplotype
                                 can be assigned before it is "fully assigned" and could not have been
                                 inherited more than that many times.
    Returns:
        sampleColumns -- a dictionary with structure like:
                         {
                             "p1": {gtIndex1: float, gtIndex2: float}
                             "p2": {...}
                         }
                         Where the gtIndex values are the ref (0) or alt (1) codes
                         and the floats sum to 1 to indicate the proportion likelihood
                         of having inherited that genotype.
    '''
    # Set up the data structure for holding this sample's assigned alleles
    sampleColumns = {
        "p1": {
            gt: 0 for gt in p1GT
        },
        "p2": {
            gt: 0 for gt in p2GT
        }
    }
    
    # Order assignable parent allele by count; ensures that least common alleles are assigned first
    Palleles = Counter([ allele for parent in [p1GT, p2GT] for allele in parent ])
    Palleles = Palleles.most_common()
    Palleles.sort(key=lambda x: x[1]) # sort by count for rare alleles first
    Palleles = { allele: count for allele, count in Palleles }
    
    # Order sample alleles by the parent allele counts
    Salleles = sorted(sampleGT, key=lambda x: Palleles[x])
    
    # Assign sample alleles to most likely parent alleles
    isValidSample = True
    parentAssigned = { parent: 0 for parent in sampleColumns } # how many alleles have been assigned to each parent
    for Sallele in Salleles: # for each sample allele
        # Detect progeny mismatch with parent genotype
        count = Palleles[Sallele]
        
        # Derive the likelihood of this allele being assigned to a specific parental haplotype
        "Fractional values indicate partial assignment of the allele to a haplotype"
        likelihood = 1 / count
        
        # Partially assign the allele to all applicable haplotypes
        for parent, assignedDict in sampleColumns.items():
            # Skip assigned parents
            if sum([ assignedDict[_allele] for _allele in assignedDict ]) == PARENT_FULLY_ASSIGNED[parent]:
                continue
            
            if Sallele in assignedDict:
                assignedDict[Sallele] += likelihood
                
                # Update the parent allele counts to reflect that one allele has been confidently assigned
                if assignedDict[Sallele] == PARENT_FULLY_ASSIGNED[parent]:
                    parentGT = list(sampleColumns[parent].keys())
                    for pGT in parentGT:
                        Palleles[pGT] -= 1
    
    return sampleColumns

def pan_plot_het(args):
    vcf = VCFTopia(args.vcfFile)
    
    try:
        p1Index = vcf.samples.index(args.parents[0])
        p2Index = vcf.samples.index(args.parents[1])
    except ValueError:
        raise ValueError(f"One or both arent IDs '{args.parents}' do not occur within the VCF sample IDs '{vcf.samples}'")
    
    try:
        referenceIndex = vcf.samples.index(args.reference)
    except ValueError:
        raise ValueError(f"Reference ID '{args.reference}' does not occur within the VCF sample IDs '{vcf.samples}'")
    
    assignments = {}
    samplesToPlot = set()
    for variant in vcf:
        as_pangenome_variant(variant, referenceIndex)
        sampleIntegers = VCFTopia.sample_integers(vcf.samples, variant)
        
        # Obtain the parent genotypes
        p1GT = variant.genotypes[p1Index][:-1] # drop the phased boolean
        p2GT = variant.genotypes[p2Index][:-1]
        
        # Skip if we can't use this variant site for het checking
        if -1 in p1GT or -1 in p2GT:
            continue
        
        # Skip if we couldn't have inherited a heterozygous genotype at this site
        possibleGTs = possible_genotypes(p1GT, p2GT)
        if all([len(gt) == 1 for gt in possibleGTs ]):
            continue
        
        # Check each sample for hom/het inheritance
        for sampleID, sampleGT in sampleIntegers.items():
            if sampleID in args.parents:
                continue
            
            sampleGTset = set(sampleGT)
            if not sampleGTset in possibleGTs:
                continue
            
            # Note genotype inheritance
            assignments.setdefault(variant.CHROM, {})
            assignments[variant.CHROM].setdefault(sampleID, [])
            if len(sampleGTset) == 1:
                assignments[variant.CHROM][sampleID].append([variant.POS, 1]) # 1 allele (hom)
            else:
                assignments[variant.CHROM][sampleID].append([variant.POS, 2]) # 2 alleles (het)
            samplesToPlot.add(sampleID)
    
    # Set up the plot object
    nrow = len(samplesToPlot)
    ncol = len(assignments)
    
    fig, ax = plt.subplots(nrows=nrow, # samples
                           ncols=ncol, # contigs
                           figsize=(10*ncol, 5*nrow))
    fig.tight_layout()
    ax = np.reshape(ax, (nrow, ncol)) # ensure shape is as expected
    
    colLabels = list(assignments.keys())
    for topColAx, label in zip(ax[0], colLabels):
        topColAx.set_title(label, fontweight="bold")
    
    # Populate plot with data
    for colIndex, (contig, speciesDict) in enumerate(assignments.items()):
        for rowIndex, (species, hapValues) in enumerate(speciesDict.items()):
            # Extract values out of the list of lists
            x, y = [], []
            for pos, _y in hapValues:
                x.append(pos)
                y.append(_y)
            
            # Get the smoothed value
            y = WMA(np.array(y), 5) # smoothing window of 5 variants
            
            # Plot as overlapping lines
            ax[rowIndex, colIndex].plot(x, y,)
            
            # Set row axis label
            if colIndex == 0:
                ax[rowIndex, colIndex].set_ylabel(species)
    
    # Show the colour scale legend and save
    fig.savefig(args.outputFileName, dpi=300, bbox_inches="tight")
    
def pan_plot_inherit(args):
    vcf = VCFTopia(args.vcfFile)
    
    try:
        p1Index = vcf.samples.index(args.parents[0])
        p2Index = vcf.samples.index(args.parents[1])
    except ValueError:
        raise ValueError(f"One or both arent IDs '{args.parents}' do not occur within the VCF sample IDs '{vcf.samples}'")
    
    try:
        referenceIndex = vcf.samples.index(args.reference)
    except ValueError:
        raise ValueError(f"Reference ID '{args.reference}' does not occur within the VCF sample IDs '{vcf.samples}'")
    
    assignments = {}
    samplesToPlot = set()
    for variant in vcf:
        as_pangenome_variant(variant, referenceIndex)
        sampleIntegers = VCFTopia.sample_integers(vcf.samples, variant)
        
        # Obtain the parent genotypes
        p1GT = variant.genotypes[p1Index][:-1] # drop the phased boolean
        p2GT = variant.genotypes[p2Index][:-1]
        
        # Skip if we can't use this variant site for inheritance inference
        if -1 in p1GT or -1 in p2GT:
            continue
        
        # For each other sample, assign alleles to their possibly inherited parent haplotype
        PARENT_FULLY_ASSIGNED = { f"p{i+1}" : int(len(gt) / 2) for i, gt in enumerate([p1GT, p2GT]) }
        possibleGTs = possible_genotypes(p1GT, p2GT)
        for sampleID, sampleGT in sampleIntegers.items():
            if sampleID in args.parents:
                continue
            
            sampleGTset = set(sampleGT)
            if not sampleGTset in possibleGTs:
                continue
            
            sampleColumns = assign_sampleGT_to_parent_haplotypes(sampleGT, p1GT, p2GT, PARENT_FULLY_ASSIGNED)
            
            # Store data of haplotype assignments for this sample
            assignments.setdefault(variant.CHROM, {})
            assignments[variant.CHROM].setdefault(sampleID, [])
            assignments[variant.CHROM][sampleID].append(
                [
                    variant.POS,
                    sampleColumns["p1"][p1GT[0]],
                    sampleColumns["p1"][p1GT[1]],
                    sampleColumns["p2"][p2GT[0]],
                    sampleColumns["p2"][p2GT[1]]
                ]
            )
            samplesToPlot.add(sampleID)
    
    # Set up the plot object
    nrow = len(samplesToPlot)
    ncol = len(assignments)
    
    fig, ax = plt.subplots(nrows=nrow, # samples
                           ncols=ncol, # contigs
                           figsize=(10*ncol, 5*nrow))
    fig.tight_layout()
    ax = np.reshape(ax, (nrow, ncol)) # ensure shape is as expected
    
    colLabels = list(assignments.keys())
    for topColAx, label in zip(ax[0], colLabels):
        topColAx.set_title(label, fontweight="bold")
    
    # Populate plot with data
    for colIndex, (contig, speciesDict) in enumerate(assignments.items()):
        for rowIndex, (species, hapValues) in enumerate(speciesDict.items()):
            # Extract values out of the list of lists
            x, h1, h2, h3, h4 = [], [], [], [], []
            for pos, _h1, _h2, _h3, _h4 in hapValues:
                x.append(pos)
                h1.append(_h1)
                h2.append(_h2)
                h3.append(_h3)
                h4.append(_h4)
            
            # Get the smoothed value for each haplotype
            h1 = WMA(np.array(h1), 5) # smoothing window of 5 variants
            h2 = WMA(np.array(h2), 5)
            h3 = WMA(np.array(h3), 5)
            h4 = WMA(np.array(h4), 5)
            
            # Plot as overlapping lines
            ax[rowIndex, colIndex].plot(x, h1,
                                        color="#D81B60", alpha=0.3,
                                        label=f"{args.parents[0]}_hap1" if colIndex==0 and rowIndex==0 else None)
            ax[rowIndex, colIndex].plot(x, h2,
                                        color="#1E88E5", alpha=0.3,
                                        label=f"{args.parents[0]}_hap2" if colIndex==0 and rowIndex==0 else None)
            ax[rowIndex, colIndex].plot(x, h3,
                                        color="#FFC107", alpha=0.3,
                                        label=f"{args.parents[1]}_hap1" if colIndex==0 and rowIndex==0 else None)
            ax[rowIndex, colIndex].plot(x, h4,
                                        color="#004D40", alpha=0.3,
                                        label=f"{args.parents[1]}_hap2" if colIndex==0 and rowIndex==0 else None)
            
            # Set row axis label
            if colIndex == 0:
                ax[rowIndex, colIndex].set_ylabel(species)
    
    # Show the colour scale legend and save
    fig.legend(loc="center left", bbox_to_anchor=(1, 0.5), ncol=1)
    fig.savefig(args.outputFileName, dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    pass
