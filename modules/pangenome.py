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

import os, sys, math
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from collections import Counter

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from vcftopia import VCFTopia

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

# Painting states are indexed as 2*i + j where i is the parent 1 haplotype
# and j is the parent 2 haplotype inherited by the offspring
PAINT_P1_HAP = np.array([0, 0, 1, 1])
PAINT_P2_HAP = np.array([0, 1, 0, 1])

def paint_parse_contig(vcf, contig, p1Index, p2Index, referenceIndex, offspringIndices, minGQ):
    '''
    Parses the variants of a contig into arrays suitable for haplotype painting. Only sites
    where both parents are fully genotyped and at least one parent is heterozygous are retained.
    
    Parameters:
        vcf -- a VCFTopia object
        contig -- a string indicating the contig to parse
        p1Index / p2Index -- integer indices of the parent samples
        referenceIndex -- integer index of the reference (haploid column) parent sample
        offspringIndices -- a list of integer indices of the offspring samples
        minGQ -- an integer indicating the minimum GQ for an offspring genotype to be used
    Returns:
        positions -- a numpy array of shape (nSites,) with variant positions
        parentGTs -- a numpy array of shape (nSites, 4) with the alleles of parent 1
                     haplotypes 1 and 2, followed by parent 2 haplotypes 1 and 2
        offspringGTs -- a numpy array of shape (nSites, nOffspring, 2) with each offspring's
                        alleles sorted in ascending order, and -1 indicating a missing
                        or filtered genotype
    '''
    positions, parentGTs, offspringGTs = [], [], []
    for variant in vcf.query(contig):
        as_pangenome_variant(variant, referenceIndex)
        genotypes = variant.genotypes
        
        # Skip sites where parents are not fully genotyped
        p1GT = genotypes[p1Index][:-1] # drop the phased boolean
        p2GT = genotypes[p2Index][:-1]
        if len(p1GT) != 2 or len(p2GT) != 2 or -1 in p1GT or -1 in p2GT:
            continue
        
        # Skip sites which cannot inform haplotype inheritance
        if p1GT[0] == p1GT[1] and p2GT[0] == p2GT[1]:
            continue
        
        # Obtain offspring genotypes that pass filtering
        GQs = variant.format("GQ") if minGQ > 0 else None
        siteGTs = []
        for offspringIndex in offspringIndices:
            sampleGT = genotypes[offspringIndex][:-1]
            if len(sampleGT) != 2 or -1 in sampleGT or (GQs is not None and GQs[offspringIndex][0] < minGQ):
                siteGTs.append([-1, -1])
            else:
                siteGTs.append(sorted(sampleGT))
        
        positions.append(variant.POS)
        parentGTs.append([*p1GT, *p2GT])
        offspringGTs.append(siteGTs)
    
    return (np.array(positions, dtype=np.int64),
            np.array(parentGTs, dtype=np.int32).reshape(-1, 4),
            np.array(offspringGTs, dtype=np.int32).reshape(-1, len(offspringIndices), 2))

def paint_state_matches(parentGTs, offspringGTs):
    '''
    Determines which painting states are consistent with each offspring genotype.
    
    Parameters:
        parentGTs / offspringGTs -- numpy arrays as returned by paint_parse_contig()
    Returns:
        matches -- a boolean numpy array of shape (nSites, nOffspring, 4) indicating whether
                   the offspring genotype is consistent with inheriting each state's haplotypes
    '''
    stateGTs = np.stack([ parentGTs[:, PAINT_P1_HAP], parentGTs[:, 2 + PAINT_P2_HAP] ], axis=2)
    stateGTs.sort(axis=2) # (nSites, 4, 2)
    return (offspringGTs[:, :, None, :] == stateGTs[:, None, :, :]).all(axis=3)

def paint_viterbi(positions, matches, called, recombinationRate, errorRate):
    '''
    Infers the most likely sequence of inherited parental haplotypes for each offspring
    using a hidden Markov model with the four painting states. Transition probabilities
    derive from the physical distance between sites via the Haldane map function, and
    emission probabilities allow for genotyping errors.
    
    Parameters:
        positions -- a numpy array of shape (nSites,) with variant positions
        matches -- a numpy array as returned by paint_state_matches()
        called -- a boolean numpy array of shape (nSites, nOffspring) indicating which
                  offspring genotypes are present
        recombinationRate -- a float indicating the expected recombination rate in cM/Mb
        errorRate -- a float indicating the probability of an offspring genotype being erroneous
    Returns:
        path -- a numpy array of shape (nSites, nOffspring) with the most likely state at each site
    '''
    numSites, numOffspring = called.shape
    
    # Derive log emission probabilities; missing genotypes are uninformative
    logE = np.where(matches, np.log(1 - errorRate), np.log(errorRate))
    logE[~called] = 0
    
    # Derive log transition probabilities between each site and the one before it
    morgans = np.diff(positions, prepend=positions[0]) * recombinationRate / 1e8
    switchProb = np.clip(0.5 * (1 - np.exp(-2 * morgans)), 1e-300, 0.5)
    logSwitch, logStay = np.log(switchProb)[:, None, None], np.log1p(-switchProb)[:, None, None]
    p1Switch = PAINT_P1_HAP[:, None] != PAINT_P1_HAP[None, :]
    p2Switch = PAINT_P2_HAP[:, None] != PAINT_P2_HAP[None, :]
    logT = np.where(p1Switch, logSwitch, logStay) + np.where(p2Switch, logSwitch, logStay) # (nSites, 4, 4)
    
    # Run the forward pass, vectorised across offspring
    backpointers = np.zeros((numSites, numOffspring, 4), dtype=np.int8)
    delta = np.log(0.25) + logE[0]
    for t in range(1, numSites):
        candidates = delta[:, :, None] + logT[t][None, :, :]
        backpointers[t] = candidates.argmax(axis=1)
        delta = candidates.max(axis=1) + logE[t]
    
    # Trace back the most likely path
    path = np.zeros((numSites, numOffspring), dtype=np.int8)
    path[-1] = delta.argmax(axis=1)
    offspringRange = np.arange(numOffspring)
    for t in range(numSites - 1, 0, -1):
        path[t-1] = backpointers[t][offspringRange, path[t]]
    
    return path

def paint_block_counts(start, end, hap, evidencePositions, evidenceHaps):
    '''
    Counts the evidence markers within a block that support or conflict with its haplotype.
    
    Parameters:
        start / end -- integers indicating the block's first and last marker positions
        hap -- an integer indicating the block's haplotype (0 or 1)
        evidencePositions -- a sorted numpy array of evidence marker positions
        evidenceHaps -- a numpy array of the haplotype each evidence marker supports
    Returns:
        supporting / conflicting -- integer counts of evidence markers in the block
    '''
    left = np.searchsorted(evidencePositions, start, side="left")
    right = np.searchsorted(evidencePositions, end, side="right")
    supporting = int((evidenceHaps[left:right] == hap).sum())
    return supporting, int(right - left) - supporting

def paint_smooth_blocks(blocks, evidencePositions, evidenceHaps, minBlockSize, minMarkers):
    '''
    Iteratively removes the most weakly supported block that is either flanked by other blocks
    and smaller than minBlockSize, or has a net support (supporting minus conflicting markers)
    less than minMarkers. Short flanked blocks imply biologically unrealistic double crossovers and
    likely reflect locus-specific artefacts, whereas poorly supported blocks occur where the parental
    haplotypes are too similar to distinguish. A removed block is absorbed by its neighbour(s) when
    they share a haplotype or it is at a contig edge; otherwise, it is dropped to leave a gap where
    the crossover could not be precisely localised.
    
    Parameters:
        blocks -- a list of lists with format like [[start, end, hap], ...] sorted by position
        evidencePositions / evidenceHaps -- numpy arrays as described in paint_block_counts()
        minBlockSize -- an integer indicating the minimum size in base pairs of a block
                        flanked by other blocks
        minMarkers -- an integer indicating the minimum net support for a block
    Returns:
        blocks -- the smoothed list of blocks with format like
                  [[start, end, hap, supporting, conflicting], ...]
    '''
    blocks = [ [start, end, hap, *paint_block_counts(start, end, hap, evidencePositions, evidenceHaps)]
              for start, end, hap in blocks ]
    while len(blocks) > 1:
        failing = [ i for i, (start, end, hap, supporting, conflicting) in enumerate(blocks)
                   if (0 < i < len(blocks) - 1 and end - start < minBlockSize) or supporting - conflicting < minMarkers ]
        if len(failing) == 0:
            break
        weakestIndex = min(failing, key=lambda x: (blocks[x][3] - blocks[x][4], blocks[x][1] - blocks[x][0]))
        start, end = blocks[weakestIndex][:2]
        
        left = blocks[weakestIndex-1] if weakestIndex > 0 else None
        right = blocks[weakestIndex+1] if weakestIndex < len(blocks) - 1 else None
        if left is None:
            right[0] = start
            updated = right
        elif right is None:
            left[1] = end
            updated = left
        elif left[2] == right[2]:
            left[1] = right[1]
            del blocks[weakestIndex+1]
            updated = left
        else:
            updated = None
        del blocks[weakestIndex]
        
        if updated is not None:
            updated[3:] = paint_block_counts(*updated[:3], evidencePositions, evidenceHaps)
    
    return blocks

def paint_offspring_blocks(positions, evidence, pathHaps, minBlockSize, minMarkers):
    '''
    Derives smoothed haplotype blocks for each offspring's inherited haplotype from one parent.
    
    Parameters:
        positions -- a numpy array of shape (nSites,) with variant positions
        evidence -- a numpy array of shape (nSites, nOffspring) indicating the parental haplotype
                    (0 or 1) that each offspring genotype supports, or -1 if it is uninformative
        pathHaps -- a numpy array of shape (nSites, nOffspring) indicating the parental haplotype
                    (0 or 1) inferred to be inherited by each offspring
        minBlockSize / minMarkers -- integers as described in paint_smooth_blocks()
    Returns:
        offspringBlocks -- a list containing, for each offspring, a list of blocks as returned
                           by paint_smooth_blocks()
    '''
    offspringBlocks = []
    for offspringNum in range(evidence.shape[1]):
        isEvidence = evidence[:, offspringNum] != -1
        if not isEvidence.any():
            offspringBlocks.append([])
            continue
        evidencePositions = positions[isEvidence]
        evidenceHaps = evidence[isEvidence, offspringNum]
        evidencePathHaps = pathHaps[isEvidence, offspringNum]
        
        # Form blocks from runs of the same inferred haplotype, then smooth them
        runStarts = np.concatenate(([0], np.flatnonzero(np.diff(evidencePathHaps)) + 1))
        runEnds = np.concatenate((runStarts[1:], [len(evidencePathHaps)])) - 1
        blocks = [ [evidencePositions[s], evidencePositions[e], evidencePathHaps[s]] for s, e in zip(runStarts, runEnds) ]
        offspringBlocks.append(paint_smooth_blocks(blocks, evidencePositions, evidenceHaps, minBlockSize, minMarkers))
    
    return offspringBlocks

def paint_find_phase_switches(offspringBlocks, minShared, window):
    '''
    Identifies breakpoints shared by many offspring at the same location. Meiotic crossovers
    occur independently in each offspring, so shared breakpoints instead indicate switch errors
    in the parent's haplotype phasing.
    
    Parameters:
        offspringBlocks -- a list as returned by paint_offspring_blocks()
        minShared -- an integer indicating the minimum number of offspring sharing a breakpoint
                     for it to be considered a phase switch
        window -- an integer indicating the maximum distance between breakpoint midpoints
                  for them to be considered as occurring at the same location
    Returns:
        switches -- a sorted list of lists with format like [[position, numShared], ...]
    '''
    breakpoints = sorted([ ((left[1] + right[0]) // 2, offspringNum)
                           for offspringNum, blocks in enumerate(offspringBlocks)
                           for left, right in zip(blocks[:-1], blocks[1:]) ])
    
    switches = []
    clusterStart = 0
    while clusterStart < len(breakpoints):
        clusterEnd = clusterStart
        while clusterEnd + 1 < len(breakpoints) and breakpoints[clusterEnd+1][0] - breakpoints[clusterStart][0] <= window:
            clusterEnd += 1
        cluster = breakpoints[clusterStart:clusterEnd+1]
        numShared = len(set([ offspringNum for _, offspringNum in cluster ]))
        if numShared >= minShared:
            switches.append([int(np.median([ midpoint for midpoint, _ in cluster ])), numShared])
        clusterStart = clusterEnd + 1
    
    return switches

def pan_paint(args):
    vcf = VCFTopia(args.vcfFile)
    
    try:
        p1Index = vcf.samples.index(args.parents[0])
        p2Index = vcf.samples.index(args.parents[1])
    except ValueError:
        raise ValueError(f"One or both parent IDs '{args.parents}' do not occur within the VCF sample IDs '{vcf.samples}'")
    referenceIndex = vcf.samples.index(args.reference)
    
    offspringIndices = [ i for i in range(len(vcf.samples)) if not i in [p1Index, p2Index] ]
    if len(offspringIndices) == 0:
        raise ValueError(f"No offspring samples found in the VCF beyond the parents '{args.parents}'")
    offspringIDs = [ vcf.samples[i] for i in offspringIndices ]
    if args.minGQ > 0 and not "##FORMAT=<ID=GQ," in vcf.header:
        print("WARNING: VCF does not contain a GQ FORMAT field; --minGQ filtering will not be applied")
        args.minGQ = 0
    minShared = max(3, math.ceil(args.phaseShared * len(offspringIDs)))
    
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("sample\tchrom\tstart\tend\tparent\thaplotype\tsupporting\tconflicting\n")
        for contig in vcf.contigs:
            positions, parentGTs, offspringGTs = paint_parse_contig(vcf, contig, p1Index, p2Index, referenceIndex,
                                                                    offspringIndices, args.minGQ)
            if len(positions) == 0:
                continue
            
            # Drop sites where too many offspring are inconsistent with any inheritance state
            matches = paint_state_matches(parentGTs, offspringGTs)
            called = offspringGTs[:, :, 0] != -1
            numCalled = called.sum(axis=1)
            numErrors = (called & ~matches.any(axis=2)).sum(axis=1)
            keep = (numCalled > 0) & (numErrors <= args.maxSiteError * numCalled)
            positions, parentGTs, matches, called = positions[keep], parentGTs[keep], matches[keep], called[keep]
            print(f"# {contig}: retained {keep.sum()} of {len(keep)} informative sites")
            if len(positions) == 0:
                continue
            
            path = paint_viterbi(positions, matches, called, args.recombinationRate, args.errorRate)
            
            # Derive haplotype blocks for each offspring's inherited haplotype from each parent
            for parentID, stateHaps in zip(args.parents, [PAINT_P1_HAP, PAINT_P2_HAP]):
                # Haplotype is evidenced where the genotype fits one parental haplotype but not the other
                hapMatches = [ matches[:, :, stateHaps == hap].any(axis=2) for hap in [0, 1] ]
                evidence = np.where(hapMatches[0] & ~hapMatches[1], 0,
                                    np.where(hapMatches[1] & ~hapMatches[0], 1, -1))
                pathHaps = stateHaps[path]
                offspringBlocks = paint_offspring_blocks(positions, evidence, pathHaps,
                                                         args.minBlockSize, args.minMarkers)
                
                # Correct parental phase switches by swapping haplotype labels after each switch
                if args.phaseShared > 0:
                    switches = paint_find_phase_switches(offspringBlocks, minShared, args.minBlockSize)
                    if len(switches) > 0:
                        for switchPosition, numShared in switches:
                            print(f"# {contig}: corrected {parentID} phase switch at ~{switchPosition} " +
                                  f"(shared by {numShared} of {len(offspringIDs)} offspring)")
                        flip = (np.searchsorted([ switchPosition for switchPosition, _ in switches ], positions) % 2)[:, None]
                        evidence = np.where(evidence == -1, -1, evidence ^ flip)
                        pathHaps = pathHaps ^ flip
                        offspringBlocks = paint_offspring_blocks(positions, evidence, pathHaps,
                                                                 args.minBlockSize, args.minMarkers)
                
                for offspringID, blocks in zip(offspringIDs, offspringBlocks):
                    for start, end, hap, supporting, conflicting in blocks:
                        fileOut.write(f"{offspringID}\t{contig}\t{start}\t{end}\t{parentID}\t{hap+1}\t{supporting}\t{conflicting}\n")

if __name__ == "__main__":
    pass
