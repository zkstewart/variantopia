import os, sys, importlib

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file

def parse_text_list(fileName):
    ids = []
    with read_gz_file(fileName) as fileIn:
        for line in fileIn:
            ids.append(line.rstrip("\r\n\t "))
    return ids

# msa mode
class MissingArgumentError(Exception):
    pass

def validate_m(args):
    '''
    Validation for arguments common to all "msa" mode commands.
    '''
    # Validate MSA file(s)
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

def validate_m_plot(args):
    '''
    Validation for arguments common to all "msa plot" mode commands.
    '''
    # Validate numeric arguments
    if args.width < 1:
        raise ValueError("--width must be a positive integer greater than 0.")
    if args.height < 1:
        raise ValueError("--height must be a positive integer greater than 0.")
    
    # Validate metadataGroups file
    if args.statistic == "uniqueness":
        if args.metadataGroups == None:
            raise ValueError(f"'-s uniqueness' requires that you specify a file to --metadata")
        
        args.metadataGroups = os.path.abspath(args.metadataGroups)
        if not os.path.isfile(args.metadataGroups):
            raise FileNotFoundError(f"Metadata file (--metadata {args.metadataGroups}) does not exist!")

def validate_m_plot_stats(args):
    '''
    Validation for arguments used in "msa plot stats" mode.
    '''
    ALLOWED_EXTENSIONS = [".png", ".pdf", ".svg"]
    
    # Validate output file location and name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")
    
    if not any(args.outputFileName.endswith(ext) for ext in ALLOWED_EXTENSIONS):
        raise ValueError(f"Output file (-o {args.outputFileName}) must have one of the following extensions: {', '.join(ALLOWED_EXTENSIONS)}")

def validate_m_plot_alignment(args):
    '''
    Validation for arguments used in "msa plot alignment" mode.
    '''
    # Validate numeric arguments
    if args.wrapLength < 1:
        raise ValueError("--wrap must be a positive integer greater than 0.")
    
    # Validate domtblout file
    if args.domtbloutFileName != None:
        args.domtbloutFileName = os.path.abspath(args.domtbloutFileName)
        if not os.path.isfile(args.domtbloutFileName):
            raise FileNotFoundError(f"domtblout file (--domtblout {args.domtbloutFileName}) does not exist!")
        
        # Validate annotarium location
        if args.annotariumDir == None:
            raise MissingArgumentError("--annotarium is mandatory since you have specified --domtblout")
        
        args.annotariumDir = os.path.abspath(args.annotariumDir)
        if not os.path.isdir(args.annotariumDir):
            raise FileNotFoundError(f"annotarium location (--annotarium {args.annotariumDir}) does not exist!")
        
        # Validate that necessary modules are discoverable
        try:
            sys.path.append(os.path.dirname(args.annotariumDir))
            from annotarium import Domains, OverlapResolver
        except ModuleNotFoundError:
            raise ModuleNotFoundError(f"Could not import Domains and OverlapResolver from '{args.annotariumDir}'")
    
    # Handle output location
    args.outputDirectory = os.path.abspath(args.outputDirectory)
    if not os.path.exists(args.outputDirectory):
        os.makedirs(args.outputDirectory, exist_ok=True)
        print(f"# Created '{args.outputDirectory}' as part of argument validation")

def validate_m_report(args):
    '''
    Validation for arguments common to all "msa report" mode commands.
    '''
    # Validate output file location and name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_m_report_pv(args):
    '''
    Validation for arguments used in "msa report per_variant" mode.
    '''
    pass # no specific validation needed for 'msa report per_variant' mode

def validate_m_report_ps(args):
    '''
    Validation for arguments used in "msa report per_sequence" mode.
    '''
    pass # no specific validation needed for 'msa report per_sequence' mode

# vcf mode
def validate_v(args):
    '''
    Validation for arguments common to all "vcf" mode commands.
    '''
    # Validate VCF file
    args.vcfFile = os.path.abspath(args.vcfFile)
    if not os.path.isfile(args.vcfFile):
        raise FileNotFoundError(f"VCF file (-i {args.vcfFile}) does not exist!")
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_v_cn(args):
    '''
    Validation for arguments common to all "vcf copynum" mode commands.
    '''
    pass # no specific validation needed for 'vcf copynum' mode

def validate_v_cn_plot(args):
    '''
    Validation for arguments used in "copynum plot" mode.
    '''
    # Validate numeric arguments
    if args.windowSize < 0:
        raise ValueError("--window must be >= 0")
    
    # Validate output file name
    if not args.outputFileName.endswith(".html"):
        raise ValueError(f"Output file (-o {args.outputFileName}) must end in .html for 'copynum plot' mode")

def validate_v_plot(args):
    '''
    Validation for arguments common used in "vcf plot" mode.
    '''
    ALLOWED_EXTENSIONS = [".png", ".pdf", ".svg"]
    
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
    
    # Validate metadataGroups file
    if args.statistic == "uniqueness":
        if args.metadataGroups == None:
            raise ValueError(f"'-s uniqueness' requires that you specify a file to --metadata")
        
        args.metadataGroups = os.path.abspath(args.metadataGroups)
        if not os.path.isfile(args.metadataGroups):
            raise FileNotFoundError(f"Metadata file (--metadata {args.metadataGroups}) does not exist!")
    
    # Validate output file name
    if not any(args.outputFileName.endswith(ext) for ext in ALLOWED_EXTENSIONS):
        raise ValueError(f"Output file (-o {args.outputFileName}) must have one of the following extensions: {', '.join(ALLOWED_EXTENSIONS)}")

def validate_v_stats(args):
    '''
    Validation for arguments used by "vcf stats" mode.
    '''
    pass # no specific validation needed for 'vcf stats' mode

def validate_v_to(args):
    '''
    Validation for arguments common to all "vcf to" mode commands.
    '''
    pass # no specific validation needed for 'vcf stats' mode

def validate_v_to_cf(args):
    '''
    Validation for arguments used in "vcf to cf" mode.
    '''
    # Validate genome file (if relevant)
    if not args.onlySNPs:
        if args.genomeFile is None:
            raise ValueError("Genome file (-g) is required to reformat VCF -> CF without --snps option.")
        
        args.genomeFile = os.path.abspath(args.genomeFile)
        if not os.path.isfile(args.genomeFile):
            raise FileNotFoundError(f"Genome file (--genomeFile {args.genomeFile}) does not exist!")

def validate_v_to_geno(args):
    '''
    Validation for arguments used in "vcf to geno" mode.
    '''
    # Validate uncalled character argument
    if len(args.uncalledCharacter) < 1:
        raise ValueError("--uncalled must be given one or more characters.")
    
    # Validate output file location
    if not args.outputFileName.endswith(".gz"):
        args.outputFileName += ".gz"
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_v_to_pos(args):
    '''
    Validation for arguments used in "vcf to pos" mode.
    '''
    pass # no specific validation needed for 'vcf to pos' mode

def validate_v_to_table(args):
    '''
    Validation for arguments used in "vcf to table" mode.
    '''
    if (args.sampleOrderFile != None) and (not os.path.isfile(args.sampleOrderFile)):
        raise FileNotFoundError(f"Unable to locate the sample order file (--sampleOrder {args.sampleOrderFile})")

def validate_v_to_msa(args):
    '''
    Validation for arguments used in "vcf to msa" mode.
    '''
    # Validate numeric arguments
    if args.ploidy < 1:
        raise ValueError("--ploidy must be a positive integer greater than 0.")
