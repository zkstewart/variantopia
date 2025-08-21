import os, gzip
from Bio import SeqIO

from .vcf import VCFTopia

def vcf_to_cf(args):
    '''
    Handles "reformat cf" mode of variantopia.
    '''
    def onlysnps_mode(args, vcf, genotypeDict):
        # Calculate nsites for use in header format
        nsites = sum([ 1 for chrom in genotypeDict for pos in genotypeDict[chrom] ])
        
        # Write output
        with open(args.outputFileName, "w") as fileOut:
            # Write header
            fileOut.write("COUNTSFILE NPOP {0} NSITES {1}\n".format(len(vcf.samples), nsites))
            fileOut.write("CHROM POS {0}\n".format(" ".join(vcf.samples)))
            
            # Format variant genotypes
            for chrom, posDict in genotypeDict.items():
                ongoingCount = 1
                for pos, refAltAndGenotypes in posDict.items():
                    refAlt, genotypes = refAltAndGenotypes
                    
                    cfGenotypes = []
                    for g in genotypes:
                        newGenotype = [0,0,0,0]
                        g = g.replace("|", "/").replace("0", refAlt[0]).replace("1", refAlt[1]) # refAlt is always len 2 since multialleles=False
                        
                        # Count alleles in genotype
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
                                raise ValueError(f"Unknown allele encountered ('{allele}')")
                        cfGenotypes.append(str(newGenotype).replace("[", "").replace("]", "").replace(" ", ""))

                    # Write formatted genotypes to file
                    outLine = "{0} {1} {2}\n".format(chrom, ongoingCount, " ".join(cfGenotypes))
                    fileOut.write(outLine)
                    
                    # Iterate position counter
                    ongoingCount += 1
    
    def genome_mode(args, vcf, genotypeDict):
        # Calculate nsites for use in header format
        nsites = sum( seqLen for chrom, seqLen in vcf.lengths.items() )
        
        # Parse genome file
        genomeRecords = SeqIO.parse(open(args.genomeFile, 'r'), "fasta")
        
        # Write output
        with open(args.outputFileName, "w") as fileOut:
            # Write header
            fileOut.write("COUNTSFILE NPOP {0} NSITES {1}\n".format(len(vcf.samples), nsites))
            fileOut.write("CHROM POS {0}\n".format(" ".join(vcf.samples)))
            
            # Format variant genotypes
            for record in genomeRecords:
                ongoingCount = 1
                chrom = record.id
                for nucleotide in record.seq:
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
                        genotype = "2,0,0,0" # impute gaps as 'a' allele
                    else:
                        raise ValueError(f"Unknown nucleotide encountered ('{nucleotide}')")
                    
                    # Produce genotypes for non-variant position
                    if (chrom not in genotypeDict) or (ongoingCount not in genotypeDict[chrom]):
                        cfGenotypes = [genotype for i in range(len(vcf.samples))]
                    # Produce genotypes for variant position
                    else:
                        refAlt, genotypes = genotypeDict[chrom][ongoingCount]
                        cfGenotypes = []
                        for g in genotypes:
                            newGenotype = [0,0,0,0]
                            g = g.replace("|", "/").replace("0", refAlt[0]).replace("1", refAlt[1]) # refAlt is always len 2 since multialleles=False
                            
                            # Count alleles in genotype
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
                                    raise ValueError(f"Unknown allele encountered ('{allele}')")
                            cfGenotypes.append(str(newGenotype).replace("[", "").replace("]", "").replace(" ", ""))
                    
                    # Write formatted genotypes to file
                    outLine = "{0} {1} {2}\n".format(chrom, ongoingCount, " ".join(cfGenotypes))
                    fileOut.write(outLine)
                    
                    # Iterate position counter
                    ongoingCount += 1
    
    # Parse VCF into necessary data structure
    vcf = VCFTopia(args.vcfFile)
    genotypeDict = vcf.as_genotype_dict(multialleles=False, indels=False, mnps=False)
    
    # Split behaviour based on onlySNPs argument
    if args.onlySNPs:
        onlysnps_mode(args, vcf, genotypeDict)
    else:
        genome_mode(args, vcf, genotypeDict)

def vcf_to_geno(args):
    '''
    Handles "reformat geno" mode of variantopia.
    '''
    # Parse VCF into necessary data structure
    vcf = VCFTopia(args.vcfFile)
    
    # Parse VCF while writing output to file
    with gzip.open(args.outputFileName, "wt") as fileOut:
        for variant in vcf:
            ref = variant.REF
            alt = variant.ALT
            
            # Skip variants if applicable
            if (not args.keepMultiallelic) and len(variant.ALT) > 1: # skip multiallelic
                continue
            if (not args.keepIndels) and variant.is_indel: # skip indels
                continue
            
            # Format geno line for this variant
            geno = []
            for genotype in VCFTopia.format_genotype_codes(variant):
                # Edit genotype to have a consistently predictable separator
                "We don't care if the VCF is phased or not for this function"
                if genotype == ".":
                    genotype = "./."
                
                genotype = genotype.replace("/", "|")
                alleles = genotype.split("|")
                
                # Get the genotype character for this sample
                genoNum = sum([ 1 for x in alleles if x == "0" ]) # count reference alleles
                geno.append(str(genoNum))
            
            # Write the geno line to the output file
            fileOut.write("\t".join(geno) + "\n")

def vcf_to_pos(args):
    '''
    Handles "reformat pos" mode of variantopia.
    '''
    # Parse VCF into necessary data structure
    vcf = VCFTopia(args.vcfFile)
    
    # Parse VCF while writing output to file
    with open(args.outputFileName, "w") as fileOut:
        for variant in vcf:
            # Skip variants if applicable
            if (not args.keepMultiallelic) and len(variant.ALT) > 1: # skip multiallelic
                continue
            if (not args.keepIndels) and variant.is_indel: # skip indels
                continue
            
            # Write the pos line to the output file
            fileOut.write(f"{variant.CHROM}\t{variant.POS}\n")

def vcf_to_table(args):
    '''
    Handles "reformat table" mode of variantopia.
    '''
    # Parse VCF into necessary data structure
    vcf = VCFTopia(args.vcfFile)

    # Parse sample order
    if args.sampleOrderFile is None:
        sampleOrder = vcf.samples # use VCF sample order if no file provided
    else:
        sampleOrder = parse_sample_order(args.sampleOrderFile)
    
    # Write output file
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#CHROM\tPOS\tREF\tALT\t" + "\t".join(sampleOrder) + "\n")
        
        # Write each variant
        for variant in vcf:
            alt = ",".join(variant.ALT)
            sampleAlleles = VCFTopia.format_sample_alleles(vcf.samples, variant)
            
            # Reorder sample alleles
            formattedAlleles = [ sampleAlleles[sample] for sample in sampleOrder ]
            
            # Write the variant line to the output file
            fileOut.write(f"{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{alt}\t" + "\t".join(formattedAlleles) + "\n")

def vcf_to_msa(args):
    '''
    Handles "reformat msa" mode of variantopia.
    '''
    # Parse VCF into necessary data structure
    vcf = VCFTopia(args.vcfFile)
    
    # Parse VCF into MSA list
    msaDict = {}
    for variant in vcf:
        refAlt = VCFTopia.format_refAlt(variant)
        alleleLength = max([ len(allele) for allele in refAlt ])
        sampleAlleles = VCFTopia.format_sample_alleles(vcf.samples, variant)
        
        # Skip if all alleles are missing
        if all("." in alleles for alleles in sampleAlleles.values()):
            continue
        
        # Format nucleotide(s) for MSA at this position
        for sample, alleles in sampleAlleles.items():
            msaDict.setdefault(sample, [])
            
            # Recode the genotype
            if "." in alleles:
                msaDict[sample].append("-"*alleleLength*args.ploidy) # missing genotype
            else:
                msaDict[sample].append("".join([ allele.ljust(alleleLength, "-") for allele in alleles.replace("|", "/").split("/") ]))
    
    # Write the MSA to the output file
    with open(args.outputFileName, "w") as fileOut:
        for sampleID, msaGTs in msaDict.items():
            seq = "".join(msaGTs)
            fileOut.write(f">{sampleID}\n{seq}\n")
