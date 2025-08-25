import os

# Plot mode
def parse_text_list(fileName):
    ids = []
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            ids.append(line.rstrip("\r\n\t "))
    return ids

def validate_m(args):
    '''
    Validation for arguments common to all "msa" mode commands.
    '''
    # Validate MSA file
    msaFiles = []
    for location in args.msaFiles:
        location = os.path.abspath(location)
        
        if os.path.isfile(location):
            msaFiles.append(location)
        elif os.path.isdir(location):
            foundAny = False
            for f in os.listdir(location):
                if any([ f.endswith(suffix) for suffix in args.suffixes ]):
                    msaFiles.append(os.path.join(location, f))
                    foundAny = True
            if not foundAny:
                raise FileNotFoundError(f"No FASTA files found in directory '{location}' ending with a --suffix value")
        # Error out if location does not exist
        else:
            raise FileNotFoundError(f"-i file or directory '{location}' not found!")
    args.msaFiles = msaFiles
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_m_plot(args):
    '''
    Validation for arguments used in "msa plot" mode.
    '''
    ALLOWED_EXTENSIONS = [".png", ".pdf", ".svg"]
    
    # Validate numeric arguments
    if args.width < 1:
        raise ValueError("--width must be a positive integer greater than 0.")
    if args.height < 1:
        raise ValueError("--height must be a positive integer greater than 0.")
    
    # Validate output file name
    if not any(args.outputFileName.endswith(ext) for ext in ALLOWED_EXTENSIONS):
        raise ValueError(f"Output file (-o {args.outputFileName}) must have one of the following extensions: {', '.join(ALLOWED_EXTENSIONS)}")

def validate_p(args):
    '''
    Validation for arguments common to all "plot" mode commands.
    '''
    ALLOWED_EXTENSIONS = [".png", ".pdf", ".svg"]
    
    # Validate VCF file
    args.vcfFile = os.path.abspath(args.vcfFile)
    if not os.path.isfile(args.vcfFile):
        raise FileNotFoundError(f"VCF file (-i {args.vcfFile}) does not exist!")
    
    # Validate optional genome and/or GFF3 files
    if args.feature == "genes" and args.gff3File == None:
        raise ValueError("'-f genes' necessitates that --gff3 be provided.")
    if args.gff3File != None:
        args.gff3File = os.path.abspath(args.gff3File)
        if not os.path.isfile(args.gff3File):
            raise FileNotFoundError(f"GFF3 file (--gff3 {args.gff3File}) does not exist!")
    
    if args.feature == "chromosomes" and args.genomeFile == None:
        raise ValueError("'-f chromosomes' necessitates that --genome be provided.")
    if args.genomeFile != None:
        args.genomeFile = os.path.abspath(args.genomeFile)
        if not os.path.isfile(args.genomeFile):
            raise FileNotFoundError(f"Genome file (--genome {args.genomeFile}) does not exist!")
    
    # Validate --ids values
    args.ids = []
    if args.idsToPlot != None:
        for value in args.idsToPlot:
            if not os.path.isfile(value):
                args.ids.append(value)
            else:
                args.ids.extend(parse_text_list(value))
    if args.ids == []:
        args.ids = None # None will flag that no ID filtering should occur
    
    # Validate numeric arguments
    if args.windowSize < 1:
        raise ValueError("-w must be a positive integer greater than 0.")
    if args.width < 1:
        raise ValueError("--width must be a positive integer greater than 0.")
    if args.height < 1:
        raise ValueError("--height must be a positive integer greater than 0.")
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if not any(args.outputFileName.endswith(ext) for ext in ALLOWED_EXTENSIONS):
        raise ValueError(f"Output file (-o {args.outputFileName}) must have one of the following extensions: {', '.join(ALLOWED_EXTENSIONS)}")
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

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
