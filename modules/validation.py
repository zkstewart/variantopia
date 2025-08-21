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
    args.outputFileName = os.path.abspath(args.outputFileName)
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

def validate_r_geno(args):
    '''
    Validation for arguments used in "reformat geno" mode.
    '''
    # Validate uncalled character argument
    if len(args.uncalledCharacter) < 1:
        raise ValueError("--uncalled must be given one or more characters.")
    
    # Validate output file location
    if not args.outputFileName.endswith(".gz"):
        args.outputFileName += ".gz"
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_r_pos(args):
    '''
    Validation for arguments used in "reformat pos" mode.
    '''
    pass # no specific validation needed for pos mode

def validate_r_table(args):
    '''
    Validation for arguments used in "reformat table" mode.
    '''
    if (args.sampleOrderFile != None) and (not os.path.isfile(args.sampleOrderFile)):
        raise FileNotFoundError(f"Unable to locate the sample order file (--sampleOrder {args.sampleOrderFile})")

def validate_r_msa(args):
    '''
    Validation for arguments used in "reformat msa" mode.
    '''
    # Validate numeric arguments
    if args.ploidy < 1:
        raise ValueError("--ploidy must be a positive integer greater than 0.")
