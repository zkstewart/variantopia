#! python3
# variantopia.py
# Front-end interface for variant analysis tools
# that have been developed over several years in the
# Various_scripts Z.K.S repository, but which need
# to be consolidated and made more user-friendly.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import validate_r, validate_r_cf
from modules.reformat import vcf_to_cf
from _version import __version__

def main():
    usage = """%(prog)s encapsulates a variety of variant analysis tools for
    predicting or working with variants in VCF format.
    """
    # Establish main parser
    p = argparse.ArgumentParser()
    
    # Set arguments shared by subparsers
    p.add_argument("-v", "--version",
                   action="version",
                   version="variantopia.py {version}".format(version=__version__))
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser(description=usage)
    subParentParser.add_argument("-v", "--version",
                                 action="version",
                                 version="variantopia.py {version}".format(version=__version__))
    
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    # BayeScan mode
    bparser = subparsers.add_parser("bayescan",
                                    parents=[p],
                                    add_help=False,
                                    help="Handle BayeScan data")
    bparser.set_defaults(func=bmain)
    
    # Filter mode
    fparser = subparsers.add_parser("filter",
                                    parents=[p],
                                    add_help=False,
                                    help="Filter VCF data")
    fparser.set_defaults(func=fmain)
    
    # Haplotype mode
    hparser = subparsers.add_parser("haplotype",
                                    parents=[p],
                                    add_help=False,
                                    help="Haplotype analysis")
    hparser.set_defaults(func=hmain)
    
    # Plot mode
    pparser = subparsers.add_parser("plot",
                                    parents=[p],
                                    add_help=False,
                                    help="Plot VCF data")
    pparser.set_defaults(func=pmain)
    
    # Reformat mode
    rparser = subparsers.add_parser("reformat",
                                    parents=[p],
                                    add_help=False,
                                    help="Reformat VCF data")
    rparser.set_defaults(func=rmain)
    
    subReformatParsers = rparser.add_subparsers(dest="reformatMode",
                                                required=True)
    
    cfparser = subReformatParsers.add_parser("cf",
                                             parents=[p],
                                             add_help=False,
                                             help="Convert to CF format e.g., for IQ-TREE2")
    cfparser.add_argument("-i", dest="vcfFile",
                          required=True,
                          help="Location of VCF file to reformat")
    cfparser.add_argument("-o", dest="outputFileName",
                          required=True,
                          help="Location to write reformatted file")
    cfparser.add_argument("--genomeFile", dest="genomeFile",
                          required=False,
                          help="Specify genome file unless --snps is used",
                          default=None)
    cfparser.add_argument("--snps", dest="onlySNPs",
                          required=False,
                          action="store_true",
                          help="Optionally, produce a .cf with only SNPs present",
                          default=False)
    
    genoparser = subReformatParsers.add_parser("geno",
                                               parents=[p],
                                               add_help=False,
                                               help="Convert to geno format e.g., for ngsLD")
    genoparser.add_argument("-i", dest="vcfFile",
                            required=True,
                            help="Location of VCF file to reformat")
    genoparser.add_argument("-o", dest="outputFileName",
                            required=True,
                            help="Location to write reformatted file")
    
    posparser = subReformatParsers.add_parser("pos",
                                              parents=[p],
                                              add_help=False,
                                              help="Convert to pos format e.g., for tree ngsLD")
    posparser.add_argument("-i", dest="vcfFile",
                           required=True,
                           help="Location of VCF file to reformat")
    posparser.add_argument("-o", dest="outputFileName",
                           required=True,
                           help="Location to write reformatted file")
    
    tableparser = subReformatParsers.add_parser("table",
                                                parents=[p],
                                                add_help=False,
                                                help="Convert to table format e.g., for manual inspection")
    tableparser.add_argument("-i", dest="vcfFile",
                             required=True,
                             help="Location of VCF file to reformat")
    tableparser.add_argument("-o", dest="outputFileName",
                              required=True,
                              help="Location to write reformatted file")
    
    vcfmsaparser = subReformatParsers.add_parser("msa",
                                                 parents=[p],
                                                 add_help=False,
                                                 help="Convert to MSA format of just variants e.g., for tree building")
    vcfmsaparser.add_argument("-i", dest="vcfFile",
                              required=True,
                              help="Location of VCF file to reformat")
    vcfmsaparser.add_argument("-o", dest="outputFileName",
                              required=True,
                              help="Location to write reformatted file")
    
    # Stats mode
    sparser = subparsers.add_parser("stats",
                                    parents=[p],
                                    add_help=False,
                                    help="Generate statistics for a VCF file")
    sparser.set_defaults(func=rmain)
    
    # Filter-subparser arguments
    ## Required arguments
    fparser.add_argument("-m", dest="measurementType",
                         required=True,
                         choices=["ed-call", "ed-depth", "splsda"],
                         help="""Specify whether you are analysing 'ed-call' (Euclidean distance
                         of 'call' variants), 'ed-depth' (Euclidean distance of 'depth'
                         CNVs), or 'splsda' (Sparse Partial Least Squares Discriminant Analysis)
                         measurements""")
    
    # Plot-subparser arguments
    ## Required arguments
    pparser.add_argument("-i", dest="inputType",
                         required=True,
                         nargs="+",
                         choices=["call", "depth"],
                         help="""Specify one or both of 'call' and 'depth' to indicate which
                         types of results to process.""")
    
    args = subParentParser.parse_args()
    
    # Split into mode-specific functions
    if args.mode == "bayescan":
        print("## variantopia.py - bayescan ##")
        validate_b(args)
    elif args.mode == "filter":
        print("## variantopia.py - filter ##")
        validate_f(args)
    elif args.mode == "haplotype":
        print("## variantopia.py - haplotype ##")
        validate_h(args)
    elif args.mode == "plot":
        print("## variantopia.py - plot ##")
        validate_p(args)
    elif args.mode == "reformat":
        print("## variantopia.py - reformat ##")
        validate_r(args)
        rmain(args)
    elif args.mode == "stats":
        print("## variantopia.py - stats ##")
        validate_s(args)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def bmain(args):
    print("BayeScan handling complete!")

def fmain(args):
    print("Filtering complete!")

def hmain(args):
    print("Haplotype analysis complete!")

def pmain(args):
    print("Plotting complete!")

def rmain(args):
    # Split into sub-mode-specific functions
    if args.reformatMode == "cf":
        print("## VCF -> CF ##")
        validate_r_cf(args)
        vcf_to_cf(args)
    elif args.reformatMode == "geno":
        print("## VCF -> Geno ##")
        validate_r_geno(args)
    elif args.reformatMode == "pos":
        print("## VCF -> Pos ##")
        validate_r_pos(args)
    elif args.reformatMode == "table":
        print("## VCF -> Table ##")
        validate_r_table(args)
    elif args.reformatMode == "msa":
        print("## VCF -> MSA ##")
        validate_r_msa(args)
    
    print("Reformatting complete!")

if __name__ == "__main__":
    main()
