import os

from .vcf import VCFTopia

def stats_to_tsv(args):
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
