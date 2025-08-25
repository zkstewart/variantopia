#! python3
# variantopia.py
# Front-end interface for variant analysis tools
# that have been developed over several years in the
# Various_scripts Z.K.S repository, but which need
# to be consolidated and made more user-friendly.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import validate_p, \
    validate_m, validate_m_plot, validate_m_report, \
    validate_r, validate_r_cf, validate_r_geno, validate_r_pos, validate_r_table, validate_r_msa
from modules.msa import msa_to_plot, msa_to_variant_report, msa_to_sequence_report
from modules.vcf import VCFTopia
from modules.gff3 import GFF3Topia
from modules.plot import GenesPlot, ChromosomesPlot
from modules.reformat import vcf_to_cf, vcf_to_geno, vcf_to_pos, vcf_to_table, vcf_to_msa
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
    
    # MSA mode
    mparser = subparsers.add_parser("msa",
                                    parents=[p],
                                    add_help=False,
                                    help="MSA variant analysis")
    mparser.set_defaults(func=mmain)
    
    subMSAParsers = mparser.add_subparsers(dest="msaMode",
                                           required=True)
    
    msaplotparser = subMSAParsers.add_parser("plot",
                                             parents=[p],
                                             add_help=False,
                                             help="Plot MSA statistics")
    msaplotparser.add_argument("-i", dest="msaFiles",
                               required=True,
                               nargs="+",
                               help="""Specify the location(s) of MSA FASTA file(s) to plot
                               and/or directories containing MSA files which end with
                               any of the indicated --suffix values""")
    msaplotparser.add_argument("-o", dest="outputFileName",
                               required=True,
                               help="Location to write output file")
    msaplotparser.add_argument("-s", dest="statistic",
                               required=True,
                               choices=["gc", "mac", "maf", "gaprate"],
                               help="Specify which statistic to plot")
    msaplotparser.add_argument("--suffix", dest="suffixes",
                               required=False,
                               nargs="+",
                               help="""Optionally, specify one or more suffixes
                               to identify MSA files in any directories
                               specified with -i; default is '.fa .fasta .fna .faa .fas'""",
                               default=[".fa", ".fasta", ".fna", ".faa", ".fas"])
    msaplotparser.add_argument("--colour", dest="colourMap",
                               required=False,
                               choices=["viridis", "Greys", "GnBu", "RdBu"],
                               help="""Optionally, specify the colour scheme to use for the plot;
                               default is 'viridis'; refer to
                               https://matplotlib.org/stable/users/explain/colors/colormaps.html
                               for examples""",
                               default="viridis")
    msaplotparser.add_argument("--width", dest="width",
                               type=int,
                               required=False,
                               help="""Optionally, specify the output plot width (default=10)""",
                               default=10)
    msaplotparser.add_argument("--height", dest="height",
                               type=int,
                               required=False,
                               help="""Optionally, specify the output plot height (default=6)""",
                               default=6)
    
    msareportparser = subMSAParsers.add_parser("report",
                                             parents=[p],
                                             add_help=False,
                                             help="Report MSA variants")
    msareportparser.add_argument("-i", dest="msaFiles",
                                 required=True,
                                 nargs="+",
                                 help="""Specify the location(s) of MSA FASTA file(s) to plot
                                 and/or directories containing MSA files which end with
                                 any of the indicated --suffix values""")
    msareportparser.add_argument("-f", dest="reportFormat",
                                 required=True,
                                 choices=["per_variant", "per_sequence"],
                                 help="Specify the format of the output report file(s)")
    msareportparser.add_argument("-o", dest="outputFileName",
                                 required=True,
                                 help="Location to write output file")
    msareportparser.add_argument("--suffix", dest="suffixes",
                                 required=False,
                                 nargs="+",
                                 help="""Optionally, specify one or more suffixes
                                 to identify MSA files in any directories
                                 specified with -i; default is '.fa .fasta .fna .faa .fas'""",
                                 default=[".fa", ".fasta", ".fna", ".faa", ".fas"])
    msareportparser.add_argument("--reportUntilStop", dest="reportUntilStop",
                                 required=False,
                                 action="store_true",
                                 help="""Optionally provide this argument if you do not want to see variants
                                 reported for a sequence after a stop codon is encountered in that sequence
                                 (variants will be reported in other sequences if they don't encounter that
                                 stop codon)""",
                                 default=False)
    
    # Plot mode
    pparser = subparsers.add_parser("plot",
                                    parents=[p],
                                    add_help=False,
                                    help="Plot VCF data")
    pparser.set_defaults(func=pmain)
    
    pparser.add_argument("-i", dest="vcfFile",
                         required=True,
                         help="Location of VCF file to plot")
    pparser.add_argument("-o", dest="outputFileName",
                         required=True,
                         help="""Location to write plot output; file extension
                         must be one of .png, .pdf, or .svg""")
    pparser.add_argument("-s", dest="statistic",
                         required=True,
                         choices=["snpnumber", "mac", "maf", "callrate", "het"],
                         help="Specify which statistic to plot")
    pparser.add_argument("-f", dest="feature",
                         required=True,
                         choices=["genes", "chromosomes"],
                         help="""Specify which type of feature to plot""")
    pparser.add_argument("-w", dest="windowSize",
                         required=True,
                         type=int,
                         help="""Specify the window size for statistics summarisation;
                         a size of 1 will plot the statistic for each nucleotide position
                         individually, while larger sizes will summarise the statistic
                         over that many basepairs in length.""")
    pparser.add_argument("--ids", dest="idsToPlot",
                         required=False,
                         nargs="+",
                         help="""Optionally, specify the IDs of genes or chromosomes
                         to limit plotting to; alternatively, you may specify file
                         name(s) which list these IDs, one per line.""",
                         default=None)
    pparser.add_argument("--genome", dest="genomeFile",
                         required=False,
                         help="Specify genome file if '-f chromosomes' is used",
                         default=None)
    pparser.add_argument("--gff3", dest="gff3File",
                         required=False,
                         help="""Specify GFF3 file if '-f genes' is used, or if
                         you want gene locations to be annotated on a
                         '-f chromosomes' plot""",
                         default=None)
    pparser.add_argument("--colour", dest="colourMap",
                         required=False,
                         choices=["viridis", "Greys", "GnBu", "RdBu"],
                         help="""Optionally, specify the colour scheme to use for the plot;
                         default is 'viridis'; refer to
                         https://matplotlib.org/stable/users/explain/colors/colormaps.html
                         for examples""",
                         default="viridis")
    pparser.add_argument("--width", dest="width",
                         type=int,
                         required=False,
                         help="""Optionally, specify the output plot width (default=10)""",
                         default=10)
    pparser.add_argument("--height", dest="height",
                         type=int,
                         required=False,
                         help="""Optionally, specify the output plot height (default=6)""",
                         default=6)
    
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
    genoparser.add_argument("--uncalled", dest="uncalledCharacter",
                            required=False,
                            help="""Optionally, set the character to use for uncalled
                            genotypes. Default is -1.""",
                            default="-1")
    genoparser.add_argument("--keepM", dest="keepMultiallelic",
                            required=False,
                            action="store_true",
                            help="Optionally, keep multiallelic sites in the output file",
                            default=False)
    genoparser.add_argument("--keepI", dest="keepIndels",
                            required=False,
                            action="store_true",
                            help="Optionally, keep indels in the output file",
                            default=False)
    
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
    posparser.add_argument("--header", dest="header",
                           required=False,
                           action="store_true",
                           help="Optionally, add a header to the output file",
                           default=False)
    posparser.add_argument("--keepM", dest="keepMultiallelic",
                           required=False,
                           action="store_true",
                           help="Optionally, keep multiallelic sites in the output file",
                           default=False)
    posparser.add_argument("--keepI", dest="keepIndels",
                           required=False,
                           action="store_true",
                           help="Optionally, keep indels in the output file",
                           default=False)
    
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
    tableparser.add_argument("--sampleOrder", dest="sampleOrderFile",
                             required=False,
                             help="""Optionally, specify a text file listing the order of samples
                             (as columns) in the output table; default is to use the VCF header order""",
                             default=None)
    
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
    vcfmsaparser.add_argument("--ploidy", dest="ploidy",
                              required=False,
                              type=int,
                              help="Optionally, indicate the ploidy of the samples (default: 2)",
                              default=2)
    
    # Stats mode
    sparser = subparsers.add_parser("stats",
                                    parents=[p],
                                    add_help=False,
                                    help="Generate statistics for a VCF file")
    sparser.set_defaults(func=smain)
    
    args = subParentParser.parse_args()
    
    # Split into mode-specific functions
    if args.mode == "bayescan":
        print("## variantopia.py - bayescan ##")
        raise NotImplementedError("Bayescan mode is not yet implemented")
        validate_b(args)
    elif args.mode == "filter":
        print("## variantopia.py - filter ##")
        raise NotImplementedError("Filter mode is not yet implemented")
        validate_f(args)
    elif args.mode == "haplotype":
        print("## variantopia.py - haplotype ##")
        raise NotImplementedError("Haplotype mode is not yet implemented")
        validate_h(args)
    elif args.mode == "msa":
        print("## variantopia.py - msa ##")
        validate_m(args)
        mmain(args)
    elif args.mode == "plot":
        print("## variantopia.py - plot ##")
        validate_p(args) # sets args.ids
        pmain(args)
    elif args.mode == "reformat":
        print("## variantopia.py - reformat ##")
        validate_r(args)
        rmain(args)
    elif args.mode == "stats":
        print("## variantopia.py - stats ##")
        raise NotImplementedError("Stats mode is not yet implemented")
        validate_s(args)
        smain(args)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def bmain(args):
    print("BayeScan handling complete!")

def fmain(args):
    print("Filtering complete!")

def hmain(args):
    print("Haplotype analysis complete!")

def mmain(args):
    # Split into sub-mode-specific functions
    if args.msaMode == "plot":
        print("## MSA plotting ##")
        validate_m_plot(args)
        msa_to_plot(args)
    elif args.msaMode == "report":
        validate_m_report(args)
        if args.reportFormat == "per_variant":
            print("## per-variant report ##")
            msa_to_variant_report(args)
        elif args.reportFormat == "per_sequence":
            print("## per-sequence report ##")
            msa_to_sequence_report(args)
    
    print("MSA analysis complete!")

def pmain(args):
    # Load VCF and GFF3 files
    vcf = VCFTopia(args.vcfFile)
    if args.gff3File != None:
        gff3 = GFF3Topia(args.gff3File)
        gff3.create_ncls_index(["gene"])
    else:
        gff3 = None
    
    # Initialise plot object
    if args.feature == "genes":
        print("## gene statistic genegrams ##")
        plot = GenesPlot(args.statistic, args.feature, args.windowSize,
                         vcf, gff3, args.genomeFile,
                         args.width, args.height
        )
    elif args.feature == "chromosomes":
        print("## chromosome statistic ideograms ##")
        plot = ChromosomesPlot(args.statistic, args.feature, args.windowSize,
                               vcf, gff3, args.genomeFile,
                               args.width, args.height
        )
    
    # Generate plot
    plot.colourMap = args.colourMap
    plot.plot(args.outputFileName, idsToPlot=args.ids)
    
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
        vcf_to_geno(args)
    elif args.reformatMode == "pos":
        print("## VCF -> Pos ##")
        validate_r_pos(args)
        vcf_to_pos(args)
    elif args.reformatMode == "table":
        print("## VCF -> Table ##")
        validate_r_table(args)
        vcf_to_table(args)
    elif args.reformatMode == "msa":
        print("## VCF -> MSA ##")
        validate_r_msa(args)
        vcf_to_msa(args)
    
    print("Reformatting complete!")

def smain(args):
    print("Statistics analysis complete!")

if __name__ == "__main__":
    main()
