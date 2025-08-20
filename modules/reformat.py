import os

from .vcf import VCFTopia

def vcf_to_cf(args):
    '''
    Handles "reformat cf" mode of variantopia.
    '''
    # Parse VCF into necessary data structure
    vcf = VCFTopia(args.vcfFile)
    genotypeDict = vcf.as_genotype_dict(multialleles=False, indels=False, mnps=False)
    
    # Calculate nsites for use in header format
    if not args.onlySNPs:
        nsites = sum( seqLen for chrom, seqLen in vcf.lengths.items() )
    else:
        nsites = sum([ 1 for chrom in genotypeDict for pos in genotypeDict[chrom] ])
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("COUNTSFILE NPOP {0} NSITES {1}\n".format(len(vcf.samples), nsites))
        fileOut.write("CHROM POS {0}\n".format(" ".join(vcf.samples)))
        
        for chrom in vcf.chroms:
            ongoingCount = 1 # This is the genomic position
            for nucleotide in record.seq:
                # Speed up onlySNPs operation
                if onlySNPs and ((chrom not in genotypeDict) or (ongoingCount not in genotypeDict[chrom])):
                    ongoingCount += 1
                    continue
                # Generate reference genotype
                if nucleotide.lower() == "a":
                    genotype = "2,0,0,0"
                elif nucleotide.lower() == "c":
                    genotype = "0,2,0,0"
                elif nucleotide.lower() == "g":
                    genotype = "0,0,2,0"
                elif nucleotide.lower() == "t":
                    genotype = "0,0,0,2"
                elif nucleotide.lower() == "n":
                    ongoingCount += 1
                    continue # skip gaps
                else:
                    print("Unknown nucleotide encountered '{0}'".format(nucleotide))
                    print("Program must exit to prevent erroneous behaviour")
                    quit()
                
                isSNP = False
                # Produce genotypes for non-variant position
                if (chrom not in genotypeDict) or (ongoingCount not in genotypeDict[chrom]):
                    cfGenotypes = [genotype for i in range(len(samplesList))]
                # Produce genotypes for variant position
                else:
                    print("Found SNP in chrom {0}, pos {1}".format(chrom, ongoingCount))
                    ref, alt, genotypes = genotypeDict[chrom][ongoingCount]
                    if len(ref) > 1: # counts file format does not support block substitution; we must ignore this SNP
                        print("Skipping block substitution at above position")
                        cfGenotypes = [genotype for i in range(len(samplesList))] # do the same as for non-variant positions
                    else:
                        isSNP = True
                        cfGenotypes = []
                        for g in genotypes:
                            assert "|" not in g # ensure formatting compliance
                            
                            newGenotype = [0,0,0,0]
                            g = g.replace("0", ref).replace("1", alt)
                            for allele in g.split("/"):
                                if allele.lower() == "a":
                                    newGenotype[0] += 1
                                elif allele.lower() == "c":
                                    newGenotype[1] += 1
                                elif allele.lower() == "g":
                                    newGenotype[2] += 1
                                elif allele.lower() == "t":
                                    newGenotype[3] += 1
                                elif allele.lower() == ".":
                                    newGenotype[0] += 1 # impute gaps as reference allele
                                else:
                                    print("Unknown allele encountered '{0}'".format(allele))
                                    print("Program must exit to prevent erroneous behaviour")
                                    quit()
                            cfGenotypes.append(str(newGenotype).replace("[", "").replace("]", "").replace(" ", ""))
                
                # Write line to file
                if onlySNPs == False or isSNP == True:
                    outLine = "{0} {1} {2}\n".format(chrom, ongoingCount, " ".join(cfGenotypes))
                    fileOut.write(outLine)
                
                # Iterate position counter
                ongoingCount += 1
