#! python3
# variantopia.py
# Front-end interface for variant analysis tools
# that have been developed over several years in the
# Various_scripts Z.K.S repository, but which need
# to be consolidated and made more user-friendly.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import validate_m, validate_m_plot, validate_m_report, \
    validate_v, validate_v_plot, validate_v_stats, validate_v_to, \
    validate_v_to_cf, validate_v_to_geno, validate_v_to_pos, validate_v_to_table, validate_v_to_msa, \
    validate_v_cn, validate_v_cn_plot
from modules.copynum import copynum_plot
from modules.msa import msa_plot_stats, msa_to_variant_report, msa_to_sequence_report
from modules.plot import vcf_plot
from modules.reformat import vcf_to_cf, vcf_to_geno, vcf_to_pos, vcf_to_table, vcf_to_msa
from modules.stats import stats_to_tsv
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
    
    # MSA subparser
    mparser = subparsers.add_parser("msa",
                                    parents=[p],
                                    add_help=False,
                                    help="MSA variant analysis")
    mparser.set_defaults(func=mmain)
    
    subMSAParsers = mparser.add_subparsers(dest="msaMode",
                                           required=True)
    
    # MSA > plot subparser
    mplotparser = subMSAParsers.add_parser("plot",
                                           parents=[p],
                                           add_help=False,
                                           help="Plot MSA file data")
    mplotparser.set_defaults(func=mmain)
    
    subMsaPlotParsers = mplotparser.add_subparsers(dest="msaPlotMode",
                                                   required=True)
    
    # MSA > plot > stats mode
    msaplotstatsparser = subMsaPlotParsers.add_parser("stats",
                                                      parents=[p],
                                                      add_help=False,
                                                      help="Plot MSA statistics")
    msaplotstatsparser.add_argument("-i", dest="msaFiles",
                                    required=True,
                                    nargs="+",
                                    help="""Specify the location(s) of MSA FASTA file(s) to plot
                                    and/or directories containing MSA files which end with
                                    any of the indicated --suffix values""")
    msaplotstatsparser.add_argument("-o", dest="outputFileName",
                                    required=True,
                                    help="Location to write output file")
    msaplotstatsparser.add_argument("-s", dest="statistic",
                                    required=True,
                                    choices=["gc", "mac", "maf", "gaprate", "uniqueness"],
                                    help="Specify which statistic to plot")
    msaplotstatsparser.add_argument("--metadata", dest="metadataGroups",
                                    required=False,
                                    help="""If you are choosing to produce '-s uniqueness'
                                    statistics, you must provide a headerless two-column tab-separated
                                    file with a unique prefix (left) followed by its group number
                                    (1 or 2; or 0 to exclude sample from calculation)""",
                                    default=None)
    msaplotstatsparser.add_argument("--suffix", dest="suffixes",
                                    required=False,
                                    nargs="+",
                                    help="""Optionally, specify one or more suffixes
                                    to identify MSA files in any directories
                                    specified with -i; default is '.fa .fasta .fna .faa .fas'""",
                                    default=[".fa", ".fasta", ".fna", ".faa", ".fas"])
    msaplotstatsparser.add_argument("--colour", dest="colourMap",
                                    required=False,
                                    choices=["viridis", "Greys", "GnBu", "RdBu"],
                                    help="""Optionally, specify the colour scheme to use for the plot;
                                    default is 'viridis'; refer to
                                    https://matplotlib.org/stable/users/explain/colors/colormaps.html
                                    for examples""",
                                    default="viridis")
    msaplotstatsparser.add_argument("--width", dest="width",
                                    type=int,
                                    required=False,
                                    help="""Optionally, specify the output plot width (default=10)""",
                                    default=10)
    msaplotstatsparser.add_argument("--height", dest="height",
                                    type=int,
                                    required=False,
                                    help="""Optionally, specify the output plot height (default=6)""",
                                    default=6)
    
    # MSA > report subparser
    msareportparser = subMSAParsers.add_parser("report",
                                             parents=[p],
                                             add_help=False,
                                             help="Report MSA variants")
    
    msareportparser.set_defaults(func=mmain)
    
    subMsaReportParsers = msareportparser.add_subparsers(dest="msaReportMode",
                                                         required=True)
    
    # MSA > report > per_variant mode
    msareppv = subMsaReportParsers.add_parser("per_variant",
                                              parents=[p],
                                              add_help=False,
                                              help="Per-variant MSA report")
    msareppv.add_argument("-i", dest="msaFiles",
                          required=True,
                          nargs="+",
                          help="""Specify the location(s) of MSA FASTA file(s) to plot
                          and/or directories containing MSA files which end with
                          any of the indicated --suffix values""")
    msareppv.add_argument("-o", dest="outputFileName",
                          required=True,
                          help="Location to write output file")
    msareppv.add_argument("--suffix", dest="suffixes",
                          required=False,
                          nargs="+",
                          help="""Optionally, specify one or more suffixes
                          to identify MSA files in any directories
                          specified with -i; default is '.fa .fasta .fna .faa .fas'""",
                          default=[".fa", ".fasta", ".fna", ".faa", ".fas"])
    msareppv.add_argument("--reportUntilStop", dest="reportUntilStop",
                          required=False,
                          action="store_true",
                          help="""Optionally provide this argument if you do not want to see variants
                          reported for a sequence after a stop codon is encountered in that sequence
                          (variants will be reported in other sequences if they don't encounter that
                          stop codon)""",
                          default=False)
    
    # MSA > report > per_sequence mode
    msarepps = subMsaReportParsers.add_parser("per_sequence",
                                              parents=[p],
                                              add_help=False,
                                              help="Per-sequence MSA report")
    msarepps.add_argument("-i", dest="msaFiles",
                          required=True,
                          nargs="+",
                          help="""Specify the location(s) of MSA FASTA file(s) to plot
                          and/or directories containing MSA files which end with
                          any of the indicated --suffix values""")
    msarepps.add_argument("-o", dest="outputFileName",
                          required=True,
                          help="Location to write output file")
    msarepps.add_argument("--suffix", dest="suffixes",
                          required=False,
                          nargs="+",
                          help="""Optionally, specify one or more suffixes
                          to identify MSA files in any directories
                          specified with -i; default is '.fa .fasta .fna .faa .fas'""",
                          default=[".fa", ".fasta", ".fna", ".faa", ".fas"])
    msarepps.add_argument("--reportUntilStop", dest="reportUntilStop",
                          required=False,
                          action="store_true",
                          help="""Optionally provide this argument if you do not want to see variants
                          reported for a sequence after a stop codon is encountered in that sequence
                          (variants will be reported in other sequences if they don't encounter that
                          stop codon)""",
                          default=False)
    
    # VCF subparser
    vparser = subparsers.add_parser("vcf",
                                    parents=[p],
                                    add_help=False,
                                    help="VCF handling")
    vparser.set_defaults(func=vmain)
    
    subVcfParsers = vparser.add_subparsers(dest="vcfMode",
                                           required=True)
    
    # VCF > copynum subparser
    vcnparser = subVcfParsers.add_parser("copynum",
                                         parents=[p],
                                         add_help=False,
                                         help="Copy number analysis")
    vcnparser.set_defaults(func=vmain)
    
    subVcfCopynumParsers = vcnparser.add_subparsers(dest="vcfCnMode",
                                                    required=True)
    
    # VCF > copynum > plot mode
    vcnplotparser = subVcfCopynumParsers.add_parser("plot",
                                                   parents=[p],
                                                   add_help=False,
                                                   help="Plot copy number inferences")
    vcnplotparser.add_argument("-i", dest="vcfFile",
                               required=True,
                               help="Location of VCF file to plot copy numbers for")
    vcnplotparser.add_argument("-o", dest="outputFileName",
                               required=True,
                               help="Location to write output file; must end in .html")
    vcnplotparser.add_argument("--window", dest="windowSize",
                               required=False,
                               type=int,
                               help="""Optionally specify a window size within which to obtain
                               a median copy number value to smooth the data; default==10000
                               but set this to 0 or 1 to turn off window smoothing""",
                               default=10000)
    
    # VCF > plot mode
    vplotparser = subVcfParsers.add_parser("plot",
                                           parents=[p],
                                           add_help=False,
                                           help="Plot VCF data")
    vplotparser.set_defaults(func=vmain)
    
    vplotparser.add_argument("-i", dest="vcfFile",
                             required=True,
                             help="Location of VCF file to plot")
    vplotparser.add_argument("-o", dest="outputFileName",
                             required=True,
                             help="""Location to write plot output; file extension
                             must be one of .png, .pdf, or .svg""")
    vplotparser.add_argument("-s", dest="statistic",
                             required=True,
                             choices=["snpnumber", "mac", "maf", "callrate", "het"],
                             help="Specify which statistic to plot")
    vplotparser.add_argument("-f", dest="feature",
                             required=True,
                             choices=["genes", "chromosomes"],
                             help="""Specify which type of feature to plot""")
    vplotparser.add_argument("-w", dest="windowSize",
                             required=True,
                             type=int,
                             help="""Specify the window size for statistics summarisation;
                             a size of 1 will plot the statistic for each nucleotide position
                             individually, while larger sizes will summarise the statistic
                             over that many basepairs in length.""")
    vplotparser.add_argument("--ids", dest="idsToPlot",
                             required=False,
                             nargs="+",
                             help="""Optionally, specify the IDs of genes or chromosomes
                             to limit plotting to; alternatively, you may specify file
                             name(s) which list these IDs, one per line.""",
                             default=None)
    vplotparser.add_argument("--genome", dest="genomeFile",
                             required=False,
                             help="Specify genome file if '-f chromosomes' is used",
                             default=None)
    vplotparser.add_argument("--gff3", dest="gff3File",
                             required=False,
                             help="""Specify GFF3 file if '-f genes' is used, or if
                             you want gene locations to be annotated on a
                             '-f chromosomes' plot""",
                             default=None)
    vplotparser.add_argument("--colour", dest="colourMap",
                             required=False,
                             choices=["viridis", "Greys", "GnBu", "RdBu"],
                             help="""Optionally, specify the colour scheme to use for the plot;
                             default is 'viridis'; refer to
                             https://matplotlib.org/stable/users/explain/colors/colormaps.html
                             for examples""",
                             default="viridis")
    vplotparser.add_argument("--width", dest="width",
                             type=int,
                             required=False,
                             help="""Optionally, specify the output plot width (default=10)""",
                             default=10)
    vplotparser.add_argument("--height", dest="height",
                             type=int,
                             required=False,
                             help="""Optionally, specify the output plot height (default=6)""",
                             default=6)
    
    # VCF > stats mode
    vstatsparser = subVcfParsers.add_parser("stats",
                                            parents=[p],
                                            add_help=False,
                                            help="Generate statistics for a VCF file")
    vstatsparser.set_defaults(func=vmain)
    vstatsparser.add_argument("-i", dest="vcfFile",
                              required=True,
                              help="Location of VCF file")
    vstatsparser.add_argument("-o", dest="outputFileName",
                              required=True,
                              help="Location to write statistics output")
    
    # VCF > filter mode
    vfilterparser = subVcfParsers.add_parser("filter",
                                             parents=[p],
                                             add_help=False,
                                             help="Filter VCF data")
    vfilterparser.set_defaults(func=vmain)
    
    # VCF > haplotype mode
    vhapparser = subVcfParsers.add_parser("haplotype",
                                          parents=[p],
                                          add_help=False,
                                          help="Haplotype analysis")
    vhapparser.set_defaults(func=vmain)
    
    # VCF > to subparser
    vtoparser = subVcfParsers.add_parser("to",
                                         parents=[p],
                                         add_help=False,
                                         help="Reformat VCF data")
    vtoparser.set_defaults(func=vmain)
    
    subVcfToParsers = vtoparser.add_subparsers(dest="vcfToMode",
                                               required=True)
    
    # VCF > to > cf mode
    vtocparser = subVcfToParsers.add_parser("cf",
                                            parents=[p],
                                            add_help=False,
                                            help="Convert to CF format e.g., for IQ-TREE2")
    vtocparser.add_argument("-i", dest="vcfFile",
                            required=True,
                            help="Location of VCF file to reformat")
    vtocparser.add_argument("-o", dest="outputFileName",
                            required=True,
                            help="Location to write reformatted file")
    vtocparser.add_argument("--genomeFile", dest="genomeFile",
                            required=False,
                            help="Specify genome file unless --snps is used",
                            default=None)
    vtocparser.add_argument("--snps", dest="onlySNPs",
                            required=False,
                            action="store_true",
                            help="Optionally, produce a .cf with only SNPs present",
                            default=False)
    
    # VCF > to > geno mode
    vtogparser = subVcfToParsers.add_parser("geno",
                                            parents=[p],
                                            add_help=False,
                                            help="Convert to geno format e.g., for ngsLD")
    vtogparser.add_argument("-i", dest="vcfFile",
                            required=True,
                            help="Location of VCF file to reformat")
    vtogparser.add_argument("-o", dest="outputFileName",
                            required=True,
                            help="Location to write reformatted file")
    vtogparser.add_argument("--uncalled", dest="uncalledCharacter",
                            required=False,
                            help="""Optionally, set the character to use for uncalled
                            genotypes. Default is -1.""",
                            default="-1")
    vtogparser.add_argument("--keepM", dest="keepMultiallelic",
                            required=False,
                            action="store_true",
                            help="Optionally, keep multiallelic sites in the output file",
                            default=False)
    vtogparser.add_argument("--keepI", dest="keepIndels",
                            required=False,
                            action="store_true",
                            help="Optionally, keep indels in the output file",
                            default=False)
    
    # VCF > to > pos mode
    vtopparser = subVcfToParsers.add_parser("pos",
                                              parents=[p],
                                              add_help=False,
                                              help="Convert to pos format e.g., for tree ngsLD")
    vtopparser.add_argument("-i", dest="vcfFile",
                           required=True,
                           help="Location of VCF file to reformat")
    vtopparser.add_argument("-o", dest="outputFileName",
                           required=True,
                           help="Location to write reformatted file")
    vtopparser.add_argument("--header", dest="header",
                           required=False,
                           action="store_true",
                           help="Optionally, add a header to the output file",
                           default=False)
    vtopparser.add_argument("--keepM", dest="keepMultiallelic",
                           required=False,
                           action="store_true",
                           help="Optionally, keep multiallelic sites in the output file",
                           default=False)
    vtopparser.add_argument("--keepI", dest="keepIndels",
                           required=False,
                           action="store_true",
                           help="Optionally, keep indels in the output file",
                           default=False)
    
    # VCF > to > table mode
    vtotparser = subVcfToParsers.add_parser("table",
                                            parents=[p],
                                            add_help=False,
                                            help="Convert to table format e.g., for manual inspection")
    vtotparser.add_argument("-i", dest="vcfFile",
                            required=True,
                            help="Location of VCF file to reformat")
    vtotparser.add_argument("-o", dest="outputFileName",
                            required=True,
                            help="Location to write reformatted file")
    vtotparser.add_argument("--sampleOrder", dest="sampleOrderFile",
                            required=False,
                            help="""Optionally, specify a text file listing the order of samples
                            (as columns) in the output table; default is to use the VCF header order""",
                            default=None)
    
    # VCF > to > msa mode
    vtomparser = subVcfToParsers.add_parser("msa",
                                            parents=[p],
                                            add_help=False,
                                            help="Convert to MSA format of just variants e.g., for tree building")
    vtomparser.add_argument("-i", dest="vcfFile",
                            required=True,
                            help="Location of VCF file to reformat")
    vtomparser.add_argument("-o", dest="outputFileName",
                            required=True,
                            help="Location to write reformatted file")
    vtomparser.add_argument("--ploidy", dest="ploidy",
                            required=False,
                            type=int,
                            help="Optionally, indicate the ploidy of the samples (default: 2)",
                            default=2)
    
    args = subParentParser.parse_args()
    
    # Split into mode-specific functions
    if args.mode == "bayescan":
        print("## variantopia.py - BayeScan results handling ##")
        raise NotImplementedError("Bayescan mode is not yet implemented")
        validate_b(args)
    elif args.mode == "msa":
        print("## variantopia.py - MSA file handling ##")
        validate_m(args)
        mmain(args)
    elif args.mode == "vcf":
        print("## variantopia.py - VCF file handling ##")
        validate_v(args)
        rmain(args)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def bmain(args):
    print("BayeScan handling complete!")

def mmain(args):
    # Split into sub-mode-specific functions
    if args.msaMode == "plot":
        validate_m_plot(args)
        if args.msaPlotMode == "stats":
            print("## Plot MSA statistics ##")
            msa_plot_stats(args)
        if args.msaPlotMode == "alignment":
            print("## Plot MSA alignment ##")
            msa_plot_alignment(args)
    
    elif args.msaMode == "report":
        validate_m_report(args)
        if args.msaReportMode == "per_variant":
            print("## Per-variant MSA variant report ##")
            msa_to_variant_report(args)
        elif args.msaReportMode == "per_sequence":
            print("## Per-sequence MSA variant report ##")
            msa_to_sequence_report(args)
    
    print("MSA analysis complete!")

def vmain(args):
    # Split into sub-mode-specific functions
    if args.vcfMode == "copynum":
        validate_v_cn(args)
        if args.vcfCnMode == "plot":
            print("## VCF copy number inference plots ##")
            validate_v_cn_plot(args)
            copynum_plot(args)
    if args.vcfMode == "plot":
        print("## Plot VCF variant information ##")
        validate_v_plot(args)
        vcf_plot(args)
    if args.vcfMode == "stats":
        print("## Generate VCF file statistics ##")
        validate_v_stats(args)
        stats_to_tsv(args)
    if args.vcfMode == "filter":
        print("## Filter VCF file ##")
        raise NotImplementedError("'vcf filter' mode not yet implemented")
    if args.vcfMode == "haplotype":
        print("## VCF haplotype analysis ##")
        raise NotImplementedError("'vcf haplotype' mode not yet implemented")
    if args.vcfMode == "to":
        validate_v_to(args)
        if args.vcfToMode == "cf":
            print("## VCF to CF conversion ##")
            validate_v_to_cf(args)
            vcf_to_cf(args)
        elif args.vcfToMode == "geno":
            print("## VCF to Geno conversion ##")
            validate_v_to_geno(args)
            vcf_to_geno(args)
        elif args.vcfToMode == "pos":
            print("## VCF to Pos conversion ##")
            validate_v_to_pos(args)
            vcf_to_pos(args)
        elif args.vcfToMode == "table":
            print("## VCF to Table conversion ##")
            validate_v_to_table(args)
            vcf_to_table(args)
        elif args.vcfToMode == "msa":
            print("## VCF to MSA conversion ##")
            validate_v_to_msa(args)
            vcf_to_msa(args)
    
    print("VCF handling complete!")

if __name__ == "__main__":
    main()
