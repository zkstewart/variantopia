import os

# Reformat mode
def validate_r(args):
    '''
    Validation for arguments common to all "reformat" mode commands.
    '''
    # Validate VCF file
    args.vcfFile = os.path.abspath(args.vcfFile)
    if not os.path.isfile(args.vcfFile):
        raise FileNotFoundError(f"VCF file (-i {args.vcfFile}) does not exist!")

    # Validate output file name
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_r_cf(args):
    '''
    Validation for arguments used in "reformat cf" mode.
    '''
    # Validate genome file (if relevant)
    if not args.onlySNPs:
        if args.genomeFile is None:
            raise ValueError("Genome file (-g) is required to reformat VCF -> CF without --snps option.")
        
        args.genomeFile = os.path.abspath(args.genomeFile)
        if not os.path.isfile(args.genomeFile):
            raise FileNotFoundError(f"Genome file (--genomeFile {args.genomeFile}) does not exist!")
