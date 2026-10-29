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

import os, sys, shutil, importlib
from matplotlib.colors import is_color_like

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, parse_text_list, parse_memory_str
from subprocess_runners import check_r_packages, r_package_version
from saige import SAIGE, SAIGE_SUBMODULE_DIR, parse_pc_set, find_gwas_manifest
from annotarium_importers import import_annotarium_gff3, import_annotarium_domains

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
            raise FileNotFoundError(f"domtblout file (--domtblout {args.domtbloutFileName}) is not a file or does not exist.")
        
        # Validate annotarium location
        if args.annotariumDir == None:
            raise MissingArgumentError("--annotarium is mandatory since you have specified --domtblout")
        
        args.annotariumDir = os.path.abspath(args.annotariumDir)
        if not os.path.isdir(args.annotariumDir):
            raise FileNotFoundError(f"annotarium location (--annotarium {args.annotariumDir}) does not exist!")
        
        # Validate that necessary modules are discoverable
        try:
            Domains, OverlapResolver = import_annotarium_domains(args.annotariumDir)
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

# pangenome mode
def validate_pan(args):
    '''
    Validation for arguments common to all "pangenome" mode commands.
    '''
    # Validate pangenome VCF file
    args.vcfFile = os.path.abspath(args.vcfFile)
    if not os.path.isfile(args.vcfFile):
        raise FileNotFoundError(f"VCF file (-i {args.vcfFile}) does not exist!")
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_pan_paint(args):
    '''
    Validation for arguments used in "pangenome paint" mode.
    '''
    # Validate sample arguments
    if args.parents[0] == args.parents[1]:
        raise ValueError(f"Parents (-p {args.parents}) must be two different samples for 'pangenome paint' mode")
    if not args.reference in args.parents:
        raise ValueError(f"Reference (-r {args.reference}) must be one of the parents (-p {args.parents})")
    
    # Validate numeric arguments
    if args.recombinationRate <= 0:
        raise ValueError("--rate must be greater than 0")
    if not 0 < args.errorRate < 0.5:
        raise ValueError("--error must be greater than 0 and less than 0.5")
    if args.minGQ < 0:
        raise ValueError("--minGQ must be >= 0")
    if not 0 <= args.maxSiteError <= 1:
        raise ValueError("--maxSiteError must be between 0 and 1")
    if args.minBlockSize < 0:
        raise ValueError("--minBlock must be >= 0")
    if args.minMarkers < 0:
        raise ValueError("--minMarkers must be >= 0")
    if not 0 <= args.phaseShared <= 1:
        raise ValueError("--phaseShared must be between 0 and 1")

def validate_pan_plot(args):
    '''
    Validation for arguments common to all "pangenome plot" mode commands.
    '''
    pass # no specific validation needed for 'pangenome progeny' mode

def validate_pan_plot_het(args):
    '''
    Validation for arguments used in "pangenome plot heterozygosity" mode.
    '''
    pass # no specific validation needed for 'pangenome plot heterozygosity' mode

def validate_pan_plot_inherit(args):
    '''
    Validation for arguments used in "pangenome plot inheritance" mode.
    '''
    pass # no specific validation needed for 'pangenome plot inheritance' mode

def validate_pan_plot_paint(args):
    '''
    Validation for arguments used in "pangenome plot paint" mode.
    '''
    # Validate paint file
    args.paintFile = os.path.abspath(args.paintFile)
    if not os.path.isfile(args.paintFile):
        raise FileNotFoundError(f"Paint file (-t {args.paintFile}) does not exist!")
    
    # Validate colours
    for colour in args.colours:
        if not is_color_like(colour):
            raise ValueError(f"--colours value '{colour}' is not a recognised colour")
    
    # Validate numeric arguments
    if args.tickInterval < 1:
        raise ValueError("--ticks must be a positive integer greater than 0")
    if args.figSize <= 0:
        raise ValueError("--size must be greater than 0")
    if args.labelSize is not None and args.labelSize <= 0:
        raise ValueError("--labelSize must be greater than 0")
    
    # Validate output file name
    ALLOWED_EXTENSIONS = [".pdf", ".png", ".svg"]
    if not any(args.outputFileName.endswith(ext) for ext in ALLOWED_EXTENSIONS):
        raise ValueError(f"Output file (-o {args.outputFileName}) must have one of the following extensions: {', '.join(ALLOWED_EXTENSIONS)}")

# plink mode
def validate_plink_runtime(args):
    '''
    Validation for plink2 runtime arguments shared by "plink" and "vcf to plink"
    mode commands (--threads, --mem, --plink2, and --tmpDir if present).
    '''
    # Validate numeric arguments
    if args.threads < 1:
        raise ValueError("--threads must be a positive integer greater than 0.")

    # Validate memory argument
    args.memory = parse_memory_str(args.memory, "Mb") # plink2 accepts memory in megabyte format
    args.memory = args.memory.strip("Mb") # but plink2 wants it just as an integer, not with 'Mb' suffix

    # Validate temporary directory location
    if getattr(args, "tmpDir", None) != None:
        args.tmpDir = os.path.abspath(args.tmpDir)
        if not os.path.isdir(args.tmpDir):
            raise FileNotFoundError(f"Temporary directory location (--tmpDir {args.tmpDir}) does not exist!")

    # Validate plink2 executable
    if args.plink2 == None:
        args.plink2 = shutil.which("plink2")
        if args.plink2 == None:
            raise FileNotFoundError("plink2 is not discoverable in your system PATH and --plink2 was not specified")
    else:
        args.plink2 = os.path.abspath(args.plink2)
        if not os.path.isfile(args.plink2):
            raise FileNotFoundError(f"plink2 executable (--plink2 {args.plink2}) does not exist!")

def validate_rscript(args, packages):
    '''
    Validation for the --rscript argument, also checking that the R packages
    needed by a utilities script are installed.
    '''
    # Validate Rscript executable
    if args.rscript == None:
        args.rscript = shutil.which("Rscript")
        if args.rscript == None:
            raise FileNotFoundError("Rscript is not discoverable in your system PATH and --rscript was not specified")
    else:
        args.rscript = os.path.abspath(args.rscript)
        if not os.path.isfile(args.rscript):
            raise FileNotFoundError(f"Rscript executable (--rscript {args.rscript}) does not exist!")

    # Validate R packages
    missing = check_r_packages(packages, args.rscript)
    if missing != []:
        raise ModuleNotFoundError(f"R package(s) {missing} are not installed for Rscript '{args.rscript}'")

def validate_plink_fileset(prefix, fileFormat, flag):
    '''
    Validation for a PLINK fileset argument, returning its absolute prefix.
    '''
    prefix = os.path.abspath(prefix)
    if fileFormat == "pgen":
        requiredSuffixes = [".pgen", ".psam"]
        if not (os.path.isfile(prefix + ".pvar") or os.path.isfile(prefix + ".pvar.zst")):
            raise FileNotFoundError(f"PLINK fileset ({flag} {prefix}) has no .pvar or .pvar.zst file!")
    elif fileFormat == "bgen":
        requiredSuffixes = [".bgen", ".sample"]
    else:
        requiredSuffixes = [".bed", ".bim", ".fam"]
    for suffix in requiredSuffixes:
        if not os.path.isfile(prefix + suffix):
            raise FileNotFoundError(f"PLINK fileset ({flag} {prefix}) has no {suffix} file!")
    return prefix

def validate_pl(args):
    '''
    Validation for arguments common to all "plink" mode commands.
    '''
    # Validate PLINK fileset
    args.plinkPrefix = validate_plink_fileset(args.plinkPrefix, args.inputFormat, "-i")

    # Validate plink2 runtime arguments
    validate_plink_runtime(args)

def validate_pl_filter(args):
    '''
    Validation for arguments used in "plink filter" mode, which "plink prune"
    mode shares.
    '''
    # Validate optional sample removal file
    if args.removeFile is not None:
        args.removeFile = os.path.abspath(args.removeFile)
        if not os.path.isfile(args.removeFile):
            raise FileNotFoundError(f"Sample removal file (--remove {args.removeFile}) does not exist!")

    # Validate output location
    args.outputPrefix = os.path.abspath(args.outputPrefix)
    if not os.path.isdir(os.path.dirname(args.outputPrefix)):
        raise FileNotFoundError(f"Output directory for -o {args.outputPrefix} does not exist!")

    # Validate numeric arguments
    if not 0 <= args.geno <= 1:
        raise ValueError("--geno must be a number from 0 to 1.")
    if not 0 <= args.maf <= 0.5:
        raise ValueError("--maf must be a number from 0 to 0.5.")

def validate_pl_output_directory(args):
    '''
    Validation for the -o output directory of "plink" mode commands, creating
    it if it does not exist.
    '''
    args.outputDirectory = os.path.abspath(args.outputDirectory)
    if not os.path.exists(args.outputDirectory):
        os.makedirs(args.outputDirectory, exist_ok=True)
        print(f"# Created '{args.outputDirectory}' as part of argument validation")
    elif not os.path.isdir(args.outputDirectory):
        raise FileExistsError(f"Output directory (-o {args.outputDirectory}) exists but is not a directory!")

def validate_pl_king(args):
    '''
    Validation for arguments used in "plink king" mode.
    '''
    # Validate output directory
    validate_pl_output_directory(args)

    # Validate optional variant file
    if args.extractFile is not None:
        args.extractFile = os.path.abspath(args.extractFile)
        if not os.path.isfile(args.extractFile):
            raise FileNotFoundError(f"Variant file (--extract {args.extractFile}) does not exist!")

    # Validate numeric arguments
    if not 0 < args.minKinship <= 0.5:
        raise ValueError("--minKinship must be a number greater than 0 and at most 0.5.")

    # Validate plot title
    if args.title == None:
        args.title = os.path.basename(args.plinkPrefix)

    # Validate Rscript and the packages needed by plot_king.R
    validate_rscript(args, ["optparse", "data.table", "ggplot2"])

def validate_pl_pca(args):
    '''
    Validation for arguments used in "plink pca" mode.
    '''
    # Validate output directory
    validate_pl_output_directory(args)

    # Validate numeric arguments
    if args.npcs < 1:
        raise ValueError("--npcs must be a positive integer greater than 0.")
    if args.plotPCs < 2:
        raise ValueError("--plotPCs must be an integer >= 2.")

    # Validate plot title
    if args.title == None:
        args.title = os.path.basename(args.plinkPrefix)

    # Validate Rscript and the packages needed by plot_pca.R
    validate_rscript(args, ["optparse", "data.table", "ggplot2"])

def validate_pl_prune(args):
    '''
    Validation for arguments used in "plink prune" mode.
    '''
    # Validate sample and variant selection arguments
    validate_pl_filter(args)

    # Validate pruning arguments
    if not 0 < args.r2 <= 1:
        raise ValueError("--r2 must be a number greater than 0 and at most 1.")
    if args.windowKb == 0:
        args.windowKb = None # skip the second pruning stage
    elif args.windowKb < 0:
        raise ValueError("--windowKb must be 0 (skip) or a positive integer.")
    if args.densityCap is not None and args.densityCap <= 0:
        raise ValueError("--densityCap must be a number greater than 0.")

def validate_pl_qc(args):
    '''
    Validation for arguments used in "plink qc" mode.
    '''
    # Validate output directory
    validate_pl_output_directory(args)

# saige mode
def validate_saige_executables(args):
    '''
    Validation for the --saigeDir and --saigeRscript arguments of "saige" mode
    commands, checking that the installed SAIGE package matches the version of
    the SAIGE scripts.
    '''
    # Validate SAIGE directory; default to the submodule, then $SAIGEDIR
    if args.saigeDir == None:
        for candidate in [SAIGE_SUBMODULE_DIR, os.environ.get("SAIGEDIR")]:
            if candidate != None and os.path.isfile(os.path.join(candidate, "extdata", "step1_fitNULLGLMM.R")):
                args.saigeDir = candidate
                break
        else:
            raise FileNotFoundError("SAIGE scripts were not found; initialise variantopia's SAIGE submodule with " +
                                    f"'git -C {os.path.dirname(os.path.dirname(SAIGE_SUBMODULE_DIR))} submodule " +
                                    "update --init external/SAIGE', or specify --saigeDir")
    args.saigeDir = os.path.abspath(args.saigeDir)
    for script in SAIGE.SCRIPTS.values():
        if not os.path.isfile(os.path.join(args.saigeDir, "extdata", script)):
            raise FileNotFoundError(f"SAIGE directory (--saigeDir {args.saigeDir}) has no extdata/{script} file!")

    # Validate SAIGE's Rscript; default to the SAIGE directory's pixi environment, then PATH
    if args.saigeRscript == None:
        pixiRscript = os.path.join(args.saigeDir, ".pixi", "envs", "default", "bin", "Rscript")
        args.saigeRscript = pixiRscript if os.path.isfile(pixiRscript) else shutil.which("Rscript")
        if args.saigeRscript == None:
            raise FileNotFoundError(f"No pixi environment was found in '{args.saigeDir}', Rscript is not " +
                                    "discoverable in your system PATH, and --saigeRscript was not specified")
    else:
        args.saigeRscript = os.path.abspath(args.saigeRscript)
        if not os.path.isfile(args.saigeRscript):
            raise FileNotFoundError(f"Rscript executable (--saigeRscript {args.saigeRscript}) does not exist!")

    # Validate SAIGE R packages and version
    missing = check_r_packages(["SAIGE", "optparse", "data.table"], args.saigeRscript)
    if missing != []:
        raise ModuleNotFoundError(f"R package(s) {missing} are not installed for Rscript '{args.saigeRscript}'")
    with open(os.path.join(args.saigeDir, "DESCRIPTION"), "r") as fileIn:
        scriptVersion = next( line.split(":", 1)[1].strip() for line in fileIn if line.startswith("Version:") )
    installedVersion = r_package_version("SAIGE", args.saigeRscript)
    if installedVersion != scriptVersion:
        raise ValueError(f"SAIGE scripts in '{args.saigeDir}' are version {scriptVersion}, but the SAIGE package " +
                         f"for Rscript '{args.saigeRscript}' is version {installedVersion}")

def validate_sa(args):
    '''
    Validation for arguments common to all "saige" mode commands.
    '''
    # Validate PLINK fileset and plink2 runtime arguments
    validate_pl(args)

    # Validate PCA file
    args.eigenvecFile = os.path.abspath(args.eigenvecFile)
    if not os.path.isfile(args.eigenvecFile):
        raise FileNotFoundError(f"PCA eigenvector file (-e {args.eigenvecFile}) does not exist!")

    # Validate output directory
    validate_pl_output_directory(args)

    # Validate numeric arguments
    if args.jobs < 1:
        raise ValueError("--jobs must be a positive integer greater than 0.")
    if not 0 <= args.relCutoff <= 1:
        raise ValueError("--relCutoff must be a number from 0 to 1.")
    if args.grmMarkers < 1:
        raise ValueError("--grmMarkers must be a positive integer greater than 0.")
    if not 0 <= args.minMAF <= 0.5:
        raise ValueError("--minMAF must be a number from 0 to 0.5.")
    if args.minMAC < 0:
        raise ValueError("--minMAC must be a number >= 0.")
    if args.noInvNormalize and args.traitType != "quantitative":
        raise ValueError("--noInvNormalize only applies to --traitType quantitative.")

    # Validate SAIGE scripts and R
    validate_saige_executables(args)

def validate_sa_calibrate(args):
    '''
    Validation for arguments used in "saige calibrate" mode.
    '''
    # Validate PC sets; e.g., '2' and '1,2' are the same set
    pcSets = [ frozenset(parse_pc_set(pcSet)[0]) for pcSet in args.pcSets ]
    if len(set(pcSets)) != len(pcSets):
        raise ValueError(f"--pcSets values must be unique sets of PCs; got {args.pcSets}")

    # Validate plot title
    if args.title == None:
        args.title = os.path.basename(args.plinkPrefix)

    # Validate Rscript and the packages needed by plot_calibration.R
    validate_rscript(args, ["optparse", "data.table", "ggplot2"])

def validate_sa_gwas(args):
    '''
    Validation for arguments used in "saige gwas" mode.
    '''
    # Validate GRM fileset
    args.grmPrefix = validate_plink_fileset(args.grmPrefix, args.grmFormat, "-g")

    # Validate PC covariates
    if args.calibrationFile != None:
        args.calibrationFile = os.path.abspath(args.calibrationFile)
        if not os.path.isfile(args.calibrationFile):
            raise FileNotFoundError(f"PC calibration file (--calibration {args.calibrationFile}) does not exist!")
    elif args.pcs != None:
        parse_pc_set(args.pcs, asList=True)
    elif args.npcs < 0:
        raise ValueError("--npcs must be an integer >= 0.")

    # Validate covariate variants
    if args.covarSNPs != None:
        args.covarSNPs = [ v.strip() for v in args.covarSNPs.split(",") if v.strip() != "" ]
        if args.covarSNPs == [] or len(set(args.covarSNPs)) != len(args.covarSNPs):
            raise ValueError("--covarSNPs must be a comma-separated list of unique variant IDs.")

def validate_sa_summarise(args):
    '''
    Validation for arguments used in "saige summarise" mode. Options not given
    are resolved later from an earlier summary of the same run, or defaults.
    '''
    # Validate 'saige gwas' output directory
    args.gwasDir = os.path.abspath(args.gwasDir)
    if not os.path.isdir(args.gwasDir):
        raise FileNotFoundError(f"'saige gwas' output directory (-i {args.gwasDir}) does not exist!")
    find_gwas_manifest(args.gwasDir)

    # Validate numeric arguments
    for value, flag, isValid, description in [
        [args.suggestive, "--suggestive", lambda x: 0 < x <= 1, "a number from 0 to 1"],
        [args.clumpP2, "--clumpP2", lambda x: 0 < x <= 1, "a number from 0 to 1"],
        [args.clumpR2, "--clumpR2", lambda x: 0 <= x <= 1, "a number from 0 to 1"],
        [args.annotateP, "--annotateP", lambda x: 0 < x <= 1, "a number from 0 to 1"],
        [args.regionGap, "--regionGap", lambda x: x >= 0, "an integer >= 0"],
        [args.topN, "--topN", lambda x: x >= 1, "a positive integer"],
        [args.clumpKb, "--clumpKb", lambda x: x >= 1, "a positive integer"],
        [args.supportKb, "--supportKb", lambda x: x >= 0, "an integer >= 0"],
        [args.radius, "--radius", lambda x: x >= 0, "an integer >= 0"]]:
        if value != None and not isValid(value):
            raise ValueError(f"{flag} must be {description}.")
    if args.linkR2 != None and args.linkR2 != "auto":
        try:
            args.linkR2 = float(args.linkR2)
        except ValueError:
            raise ValueError("--linkR2 must be 'auto' or a number from 0 to 1.")
        if not 0 <= args.linkR2 <= 1:
            raise ValueError("--linkR2 must be 'auto' or a number from 0 to 1.")
    if args.suggestive != None and args.clumpP2 != None and args.clumpP2 < args.suggestive:
        raise ValueError("--clumpP2 must be >= --suggestive.")

    # Validate annotation files
    for attr, flag in [["gff3File", "--gff3"], ["annotableFile", "--annotable"], ["chromMap", "--chromMap"]]:
        if getattr(args, attr) != None:
            setattr(args, attr, os.path.abspath(getattr(args, attr)))
            if not os.path.isfile(getattr(args, attr)):
                raise FileNotFoundError(f"File ({flag} {getattr(args, attr)}) does not exist!")
    if args.annotariumDir != None:
        args.annotariumDir = os.path.abspath(args.annotariumDir)
        if not os.path.isdir(args.annotariumDir):
            raise FileNotFoundError(f"annotarium location (--annotarium {args.annotariumDir}) does not exist!")

    # Validate plink2 runtime arguments and the Rscript for plot_manhattan.R
    validate_plink_runtime(args)
    validate_rscript(args, ["optparse", "data.table", "ggplot2"])

# vcf mode
def validate_v(args):
    '''
    Validation for arguments common to all "vcf" mode commands.
    '''
    # Validate VCF file
    args.vcfFile = os.path.abspath(args.vcfFile)
    if not os.path.isfile(args.vcfFile):
        raise FileNotFoundError(f"VCF file (-i {args.vcfFile}) does not exist!")
    
    # Validate output file name (modes writing multiple files from a prefix validate their own)
    if hasattr(args, "outputFileName"):
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

def validate_v_filter(args):
    '''
    Validation for arguments used by "vcf filter" mode.
    '''
    # Validate optional GFF3 file
    if args.gff3File != None:
        args.gff3File = os.path.abspath(args.gff3File)
        if not os.path.isfile(args.gff3File):
            raise FileNotFoundError(f"GFF3 file (--gff3 {args.gff3File}) is not a file or does not exist.")
        
        # Validate annotarium location
        if args.annotariumDir == None:
            raise MissingArgumentError("--annotarium is mandatory since you have specified --gff3")
        
        args.annotariumDir = os.path.abspath(args.annotariumDir)
        if not os.path.isdir(args.annotariumDir):
            raise FileNotFoundError(f"annotarium location (--annotarium {args.annotariumDir}) is not a directory or does not exist.")
        
        # Validate that necessary modules are discoverable
        try:
            GFF3Feature, GFF3Tarium = import_annotarium_gff3(args.annotariumDir)
        except ModuleNotFoundError:
            raise ModuleNotFoundError(f"Could not import GFF3Tarium from '{args.annotariumDir}'")
    
    # Validate that a filtering option was provided
    if not any([args.gff3File != None]): # eventually this any() check can be expanded to include other filter options
        raise ValueError("--gff3 must be provided or 'vcf filter' has nothing to do!")

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
        
        # Validate annotarium location
        if args.annotariumDir == None:
            raise MissingArgumentError("--annotarium is mandatory since you have specified --gff3")
        
        args.annotariumDir = os.path.abspath(args.annotariumDir)
        if not os.path.isdir(args.annotariumDir):
            raise FileNotFoundError(f"annotarium location (--annotarium {args.annotariumDir}) is not a directory or does not exist.")
        
        # Validate that necessary modules are discoverable
        try:
            GFF3Feature, GFF3Tarium = import_annotarium_gff3(args.annotariumDir)
        except ModuleNotFoundError:
            raise ModuleNotFoundError(f"Could not import GFF3Tarium from '{args.annotariumDir}'")
    
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

def validate_v_relabel(args):
    '''
    Validation for arguments used by "vcf relabel" mode.
    '''
    # Validate logic of behavioural arguments
    if args.sampleMap == None and args.chromMap == None:
        raise ValueError("Must provide either --sampleMap or --chromMap for this mode to have something to do!")
    
    # Validate ID mapping file(s)
    if args.sampleMap != None and ( not os.path.isfile(args.sampleMap) ):
        raise FileNotFoundError(f"Samples ID mapping file (--sampleMap {args.sampleMap}) does not exist!")
    if args.chromMap != None and ( not os.path.isfile(args.chromMap) ):
        raise FileNotFoundError(f"Chromosomes ID mapping file (--chromMap {args.chromMap}) does not exist!")

def validate_v_stats(args):
    '''
    Validation for arguments used by "vcf stats" mode.
    '''
    pass # no specific validation needed for 'vcf stats' mode

def validate_v_panel(args):
    '''
    Validation for arguments used by "vcf panel" mode.
    '''
    # Validate FASTA file
    args.fastaFile = os.path.abspath(args.fastaFile)
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError(f"Genome FASTA file (-f {args.fastaFile}) does not exist!")
    
    # Validate numeric arguments
    if args.numVariants < 1:
        raise ValueError("-n must be a positive integer greater than 0.")
    if args.threads < 1:
        raise ValueError("--threads must be a positive integer greater than 0.")

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

def validate_v_to_plink(args):
    '''
    Validation for arguments used in "vcf to plink" mode.
    '''
    # Validate fam file
    args.famFile = os.path.abspath(args.famFile)
    if not os.path.isfile(args.famFile):
        raise FileNotFoundError(f"Fam file (-f {args.famFile}) does not exist!")

    # Validate optional mapping files
    if args.chromMap is not None:
        args.chromMap = os.path.abspath(args.chromMap)
        if not os.path.isfile(args.chromMap):
            raise FileNotFoundError(f"Chromosome map file (--chromMap {args.chromMap}) does not exist!")
    if args.sampleMap is not None:
        args.sampleMap = os.path.abspath(args.sampleMap)
        if not os.path.isfile(args.sampleMap):
            raise FileNotFoundError(f"Sample map file (--sampleMap {args.sampleMap}) does not exist!")

    # Validate output location
    args.outputPrefix = os.path.abspath(args.outputPrefix)
    if not os.path.isdir(os.path.dirname(args.outputPrefix)):
        raise FileNotFoundError(f"Output directory for -o {args.outputPrefix} does not exist!")

    # Validate numeric arguments
    if args.maxAlleles == 0:
        args.maxAlleles = None # no limit
    elif args.maxAlleles < 2:
        raise ValueError("--maxAlleles must be 0 (no limit) or an integer >= 2.")
    if args.outputFormat == "bed" and args.maxAlleles != 2:
        raise ValueError("-of bed can only hold biallelic variants, so requires --maxAlleles 2.")

    # Validate plink2 runtime arguments
    validate_plink_runtime(args)

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
