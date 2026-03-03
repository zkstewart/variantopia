#! python3
# panel.py
# Functions to aid in the production of a variant panel
# from an original VCF file. Code is separated out into a separate module
# to provide tidiness as the multithreaded code takes up some space.

import os, sys, math
import numpy as np
import concurrent.futures
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from vcf import VCFTopia
from parsing import read_gz_file, get_chunking_points, BgzCapableWriter

def get_chunk_index_for_this_position(points, position):
    for chunkIndex, chunkStart in enumerate(points[:-1]):
        if position >= chunkStart and position < points[chunkIndex+1]:
            return chunkIndex
    raise Exception(f"Unhandled scenario in get_chunk_index_for_this_position(); position={position}; points={points}")

def score_variant(pos, callrate, mac, linkage):
    return callrate * (mac * (1-linkage))

def calculate_simple_linkage(variant1, variant2):
    '''
    This method was invented (?) by Z.K.S as a simple measure of how often
    a particular genotype change occurs. For example, we can count the number
    of samples which make a transition between 0/0 -> 1/1 compared to any
    other transition (e.g., 0/0 -> 0/1) and obtain the proportion of samples
    which make the most common transition as the measure of "linkage". Although
    we should always expect the value to be at least ~0.33 (for a diploid, the least
    linkage is an event split of same homo/opposite homo/hetero) it still provides
    a relevant value for comparing linkage between variants where the "less linked"
    variant will have a lower value here.
    
    Parameters:
        variant1 / variant2 -- a cyvcf2.cyvcf2.Variant object as obtained
                               when parsing a VCF file.
    '''  
    gt1 = VCFTopia.genotype_codes(variant1)
    gt2 = VCFTopia.genotype_codes(variant2)
    
    changes = {}
    for _gt1, _gt2 in zip(gt1, gt2):
        if _gt1 != "./." and _gt2 != "./.":
            changes.setdefault(_gt1, {})
            changes[_gt1].setdefault(_gt2, 0)
            changes[_gt1][_gt2] += 1
    
    linkage = 0
    for linkageDict in changes.values():
        linkedProportion = max(linkageDict.values()) / sum(linkageDict.values())
        linkage += linkedProportion
    linkage /= len(changes)
    
    return linkage

def chunk_statistics_task(vcfFile, contig, points, variantTypes, biallelicOnly):
    chunkDict = {
        contig: {
            i:[]
            for i in range(0, len(points)-1)
        }
    }
    
    vcf = VCFTopia(vcfFile)
    for variant in vcf.query(contig):
        # Skip variant if it fails variant type or biallelic filters
        filterSNP = variant.is_snp and (not "snp" in variantTypes)
        filterMNP = variant.is_mnp and (not "mnp" in variantTypes)
        filterINDEL = variant.is_indel and (not "indel" in variantTypes)
        filterMulti = len(variant.ALT) > 1 and (biallelicOnly)
        
        if filterSNP or filterMNP or filterINDEL or filterMulti:
            continue
        
        # Calculate statistics that'll help with picking variants
        callrate = VCFTopia.calculate_callrate(variant)
        mac = VCFTopia.calculate_mac(variant)
        
        # Store statistics within this variant's genome chunk
        chunkIndex = get_chunk_index_for_this_position(points, variant.POS)
        chunkDict[contig][chunkIndex].append([variant.POS, callrate, mac])
    return chunkDict

def run_vcf_chunk_statistics(vcfFile, contigs, chunkPoints, variantTypes, biallelicOnly, threads=1):
    '''
    Will run the chunk_statistics_task() on contigs from a VCF file in parallel
    
    Parameters:
        vcf -- a string pointing to a VCF file location.
        contigs -- the contigs to iterate over; they must be represented in
                   the vcf object.
        chunkPoints -- a dictionary where keys are contigs and values are a list
                       of integers indicating the boundaries of each genomic 'chunk'
        variantTypes -- a list containing strings which indicate allowed variants;
                        handled values include ["snp", "mnp", "indel"]
        biallelicOnly -- a boolean indicating whether variants should only be accepted if they
                         are biallelic (==True); otherwise, allow multiallelic variants (==False)
        threads -- an integer indicating how many threads to run.
    '''
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for contig in contigs:
            futures.append(
                executor.submit(chunk_statistics_task, vcfFile, contig, chunkPoints[contig], variantTypes, biallelicOnly)
            )
    
    chunkVariants = {}
    for f in futures:
        try:
            chunkDict = f.result()
            chunkVariants.update(chunkDict)
        except Exception as e:
            raise Exception(f"run_vcf_chunk_statistics encountered the following error:\n{e}")
    return chunkVariants

def vcf_panel(args):
    '''
    Handles "vcf panel" mode of variantopia.
    '''
    # Parse genome for contig lengths
    with read_gz_file(args.fastaFile) as fileIn:
        genomeRecords = SeqIO.parse(fileIn, "fasta")
        lengthsDict = { record.id:len(record) for record in genomeRecords }
        contigs = list(lengthsDict.keys())
    
    # Figure out how many variants to pick per contig
    perContig = math.ceil(args.numVariants / len(lengthsDict))
    print(f"# Panel will aim for {args.numVariants} variants to be picked from " +
          f"a genome with {len(contigs)} contigs")
    print(f"# With rounding, this is approximately {perContig} variant(s) per contig")
    
    # Get the chunking points for each contig
    chunkPoints = {} # this will have a list of integers where the start is max(0, lastStart) and end is the integer within the list
    for contig, length in lengthsDict.items():
        chunkPoints[contig] = [0, *get_chunking_points(length, perContig, isNumOfChunks=True), length+1] # include start and end for chunk boundaries
    
    # Parallel parsing of VCF with storage of statistics
    chunkVariants = run_vcf_chunk_statistics(args.vcfFile, contigs, chunkPoints, args.variantTypes,
                                             args.biallelicOnly, args.threads)
    
    # Select variants from within each chunk
    vcf = VCFTopia(args.vcfFile)
    chosenStats = {}
    with BgzCapableWriter(args.outputFileName) as fileOut:
        # Write VCF header
        fileOut.write(vcf.vcf.raw_header)
        
        # Iterate over and select variants
        lastVariant = None
        for contig, chunkDict in chunkVariants.items():
            lastVariant = None # new contig, new last variant
            for chunkIndex, variantsList in chunkDict.items():
                # Skip genomic chunks that were empty
                if variantsList == []:
                    continue
                
                # Calculate linkage between this variant and last chunk's selected variant
                queryVariants = vcf.query(contig, (variantsList[0][0], variantsList[-1][0]))
                thisVariant = next(queryVariants)
                for i, (pos, callrate, mac) in enumerate(variantsList):
                    while pos != thisVariant.POS:
                        thisVariant = next(queryVariants)
                    assert pos == thisVariant.POS # sanity check
                    
                    if lastVariant == None:
                        linkage = 0
                    else:
                        linkage = calculate_simple_linkage(thisVariant, lastVariant)
                    
                    variantsList[i].append(linkage)
                
                # Order the variants in this chunk by the scored combination of measured metrics
                scoredVariants = [ [*x, score_variant(*x)] for x in variantsList ]
                scoredVariants.sort(key = lambda x: -x[-1])
                
                # Select the first variant if it has a positive score
                if scoredVariants[0][1] > 0:
                    pass # selection is done after
                
                # Otherwise, select the most central variant
                else:
                    """A score of 0 means it is perfectly linked with the previous variant and strictly speaking 
                    has no benefit to this panel. However, there is a chance that future sampling will reveal a
                    difference, and hence a central variant choice gives us spread over the chromosome and greater
                    adherence to our anticipated number of variants in the final panel (as opposed to just skipping
                    the variant entirely)."""
                    chunkCentre = (chunkPoints[contig][chunkIndex] + chunkPoints[contig][chunkIndex+1]) / 2
                    scoredVariants.sort(key = lambda x: abs(chunkCentre - x[0]))
                
                chosenPos = scoredVariants[0][0]
                chosenScore = scoredVariants[0][-1]
                
                # Write the selected variant to file
                chosenVariant = [ x for x in vcf.query(contig, (chosenPos, chosenPos)) if x.POS == chosenPos ][0] # search through possible overlaps
                fileOut.write(str(chosenVariant))
                
                # Setup for next iteration and store statistics for later reporting
                lastVariant = chosenVariant
                chosenStats.setdefault(contig, [0, 0]) # [numVariants, summedScore]
                chosenStats[contig][0] += 1
                chosenStats[contig][1] += chosenScore
    
    # Print out some basic statistics
    numChosen = sum([ x[0] for x in chosenStats.values() ])
    print(f"# A total of {numChosen} variants were selected. Quick summary of each contig is as follows:")
    for contig, (numVariants, summedScore) in chosenStats.items():
        avgScore = round(summedScore / numVariants, 3) # 3 decimals, why not?
        print(f"# {contig}: {numVariants} variants with avg. score of {avgScore}")
