#!/usr/bin/env Rscript
# plot_manhattan.R
# Produce Manhattan and QQ plots from GWAS summary statistics (SAIGE by default,
# but column names are configurable so PLINK2 .glm output also works). Run by
# variantopia 'saige summarise', or standalone.
#
# Example:
#   Rscript plot_manhattan.R -i reticulata.saige.txt -o reticulata --title "C. reticulata"
#
# Requires: R packages optparse, data.table, ggplot2 (and png, for --format pdf)
#
# Outputs:
#   <prefix>.manhattan.png, <prefix>.qq.png (.pdf with --format pdf), <prefix>.top_hits.tsv,
#   <prefix>.chrom_lambda.tsv, and <prefix>.regions.tsv if any variant has
#   P < --region_p

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
})

###########################################################################
##                                                                       ##
##                         ARGUMENT PARSING                              ##
##                                                                       ##
###########################################################################

optionList <- list(
    make_option(c("-i", "--input"), type="character",
                help="GWAS results file (tab/space-delimited with header)"),
    make_option(c("-o", "--output"), type="character",
                help="Output file prefix"),
    make_option("--chr_col", type="character", default="CHR",
                help="Chromosome column [default %default]"),
    make_option("--pos_col", type="character", default="POS",
                help="Position column [default %default]"),
    make_option("--p_col", type="character", default="p.value",
                help="P-value column [default %default]"),
    make_option("--id_col", type="character", default="MarkerID",
                help="Variant ID column [default %default]"),
    make_option("--sig", type="character", default="bonferroni",
                help=paste("Genome-wide significance threshold; a P-value, 'bonferroni'",
                           "(0.05 / number of tests), or 'none' [default %default]")),
    make_option("--suggestive", type="character", default="none",
                help="Suggestive threshold; a P-value or 'none' [default %default]"),
    make_option("--meff", type="integer", default=NA,
                help=paste("If given, draw an additional line at 0.05 / meff, i.e., a Bonferroni",
                           "correction for an effective number of independent tests")),
    make_option("--chroms", type="character", default=NA,
                help="Comma-separated chromosomes to plot, in order [default all, natural sort]"),
    make_option("--title", type="character", default="",
                help="Plot title"),
    make_option("--annotate", type="double", default=NA,
                help="Label the top variant per chromosome with P below this value"),
    make_option("--top_n", type="integer", default=50,
                help="Number of top variants to write to <prefix>.top_hits.tsv [default %default]"),
    make_option("--region_p", type="double", default=1e-5,
                help=paste("Variants with P below this are merged into associated regions",
                           "written to <prefix>.regions.tsv [default %default]")),
    make_option("--region_gap", type="double", default=2e6,
                help="Maximum gap (bp) between variants merged into one region [default %default]"),
    make_option("--format", type="character", default="png",
                help="Plot file format; 'png' or 'pdf' (vector graphics, e.g., for editing) [default %default]"),
    make_option("--raster_dpi", type="integer", default=300,
                help=paste("With --format pdf, the resolution at which dense point layers are embedded",
                           "as images, keeping the PDF small enough to edit; 0 keeps them as vector",
                           "points [default %default]")),
    make_option("--width", type="double", default=12,
                help="Manhattan plot width in inches [default %default]"),
    make_option("--height", type="double", default=5,
                help="Manhattan plot height in inches [default %default]")
)
opt <- parse_args(OptionParser(option_list=optionList))
if (is.null(opt$input) || is.null(opt$output)) {
    stop("-i/--input and -o/--output are required")
}
if (!opt$format %in% c("png", "pdf")) stop("--format must be 'png' or 'pdf'")
rasterise <- opt$format == "pdf" && opt$raster_dpi > 0
if (rasterise && !requireNamespace("png", quietly=TRUE)) {
    stop("--format pdf needs the R package 'png' to embed dense points as images; install it, or use --raster_dpi 0")
}

parse_threshold <- function(value, nTests) {
    if (is.na(value) || value == "none") return(NA_real_)
    if (value == "bonferroni") return(0.05 / nTests)
    as.numeric(value)
}

###########################################################################
##                                                                       ##
##                          RASTERISED POINTS                            ##
##                                                                       ##
###########################################################################

# Dense point layers make very large PDFs (millions of vector points) that are
# slow to open and edit. rasterised_points() is a drop-in for geom_point() that
# draws its points as one embedded image at --raster_dpi, sized to its panel
# when the plot is drawn, while axes, labels, lines and text remain vectors.
# Requires the 'png' R package.
rasterised_points <- function(..., dpi=300) {
    layer <- geom_point(...)
    parent <- layer$geom
    layer$geom <- ggproto(NULL, parent,
        draw_panel=function(self, ...) {
            grob <- ggproto_parent(parent, self)$draw_panel(...)
            grid::gTree(children=grid::gList(grob), dpi=dpi, cl="rasterised_points")
        })
    layer
}
# grid draws a gTree via its makeContent() method, found here by S3 dispatch
makeContent.rasterised_points <- function(x) {
    width <- grid::convertWidth(grid::unit(1, "npc"), "in", valueOnly=TRUE)
    height <- grid::convertHeight(grid::unit(1, "npc"), "in", valueOnly=TRUE)
    imageFile <- tempfile(fileext=".png")
    currentDevice <- grDevices::dev.cur()
    grDevices::png(imageFile, width=width, height=height, units="in", res=x$dpi, bg="transparent", type="cairo")
    # Points are in native units, which span 0-1 in a panel as in a new viewport
    grid::grid.newpage()
    grid::pushViewport(grid::viewport())
    grid::grid.draw(x$children)
    grDevices::dev.off()
    grDevices::dev.set(currentDevice)
    image <- png::readPNG(imageFile)
    unlink(imageFile)
    grid::setChildren(x, grid::gList(grid::rasterGrob(image, width=grid::unit(width, "in"), height=grid::unit(height, "in"))))
}

###########################################################################
##                                                                       ##
##                           DATA LOADING                                ##
##                                                                       ##
###########################################################################

gwas <- fread(opt$input, select=c(opt$chr_col, opt$pos_col, opt$p_col, opt$id_col))
setnames(gwas, c(opt$chr_col, opt$pos_col, opt$p_col, opt$id_col), c("CHR", "POS", "P", "ID"))
gwas[, CHR := as.character(CHR)]
gwas[, P := suppressWarnings(as.numeric(P))]
gwas <- gwas[!is.na(P) & P > 0 & P <= 1]
if (nrow(gwas) == 0) stop("No valid P-values found in input")

# Order chromosomes (natural sort, or user-specified)
if (!is.na(opt$chroms)) {
    chromOrder <- strsplit(opt$chroms, ",")[[1]]
    gwas <- gwas[CHR %in% chromOrder]
} else {
    chromKeys <- unique(gwas$CHR)
    chromNum <- suppressWarnings(as.numeric(gsub("\\D", "", chromKeys)))
    chromOrder <- chromKeys[order(chromNum, chromKeys)]
}
gwas[, CHR := factor(CHR, levels=chromOrder)]
setorder(gwas, CHR, POS)
gwas[, LOGP := -log10(P)]

nTests <- nrow(gwas)
sigLine <- parse_threshold(opt$sig, nTests)
suggestiveLine <- parse_threshold(opt$suggestive, nTests)
meffLine <- if (is.na(opt$meff)) NA_real_ else 0.05 / opt$meff

# Genomic inflation factor from median chi-square (1 df)
lambdaGC <- median(qchisq(gwas$P, df=1, lower.tail=FALSE)) / qchisq(0.5, df=1)
cat(sprintf("Variants plotted: %d\n", nTests))
cat(sprintf("Genomic inflation (lambda GC): %.3f\n", lambdaGC))
cat(sprintf("Min P: %.3g\n", min(gwas$P)))
if (!is.na(sigLine)) cat(sprintf("Significance threshold: %.3g (%d variants pass)\n", sigLine, sum(gwas$P < sigLine)))
if (!is.na(meffLine)) cat(sprintf("Meff threshold: %.3g (%d variants pass)\n", meffLine, sum(gwas$P < meffLine)))
if (!is.na(suggestiveLine)) cat(sprintf("Suggestive threshold: %.3g (%d variants pass)\n", suggestiveLine, sum(gwas$P < suggestiveLine)))

# Per-chromosome inflation: a large-effect locus in a family population drags
# whole chromosome segments into association, raising lambda on its own
# chromosome(s), whereas confounding raises it on all of them
chromLambda <- gwas[, .(N_VARIANTS=.N,
                        LAMBDA_GC=round(median(qchisq(P, df=1, lower.tail=FALSE)) / qchisq(0.5, df=1), 4),
                        MIN_P=signif(min(P), 3),
                        N_SIG=if (is.na(sigLine)) NA_integer_ else sum(P < sigLine),
                        N_SUGGESTIVE=if (is.na(suggestiveLine)) NA_integer_ else sum(P < suggestiveLine)),
                    keyby=CHR]
fwrite(chromLambda, paste0(opt$output, ".chrom_lambda.tsv"), sep="\t")
cat(sprintf("Per-chromosome lambda GC: %.3f (%s) to %.3f (%s)\n",
            min(chromLambda$LAMBDA_GC), chromLambda$CHR[which.min(chromLambda$LAMBDA_GC)],
            max(chromLambda$LAMBDA_GC), chromLambda$CHR[which.max(chromLambda$LAMBDA_GC)]))

###########################################################################
##                                                                       ##
##                         MANHATTAN PLOT                                ##
##                                                                       ##
###########################################################################

# Cumulative positions along the genome
chromLengths <- gwas[, .(LEN=max(POS)), by=CHR]
chromLengths[, OFFSET := shift(cumsum(as.numeric(LEN)), fill=0)]
gwas <- merge(gwas, chromLengths[, .(CHR, OFFSET)], by="CHR", sort=FALSE)
gwas[, BPCUM := POS + OFFSET]
axisDf <- gwas[, .(CENTER=(min(BPCUM) + max(BPCUM)) / 2), by=CHR]

# Thin the uninformative bulk of points (-log10 P < 2) for speed and file size;
# points are binned on a fine grid, so the visual result is unchanged
genomeSpan <- max(gwas$BPCUM)
plotDf <- rbind(
    gwas[LOGP >= 2],
    unique(gwas[LOGP < 2][, `:=`(XBIN=round(BPCUM / genomeSpan * 4000),
                                 YBIN=round(LOGP, 2))], by=c("XBIN", "YBIN"))[, !c("XBIN", "YBIN")]
)

# For a PDF, points below the lowest threshold line are embedded as an image,
# while those above it remain vectors so that hits can be edited
if (rasterise) {
    vectorLimit <- -log10(max(c(sigLine, meffLine, suggestiveLine, 1e-5), na.rm=TRUE))
    pointLayers <- list(rasterised_points(data=plotDf[LOGP < vectorLimit], size=0.6, alpha=0.8, dpi=opt$raster_dpi),
                        geom_point(data=plotDf[LOGP >= vectorLimit], size=0.6, alpha=0.8))
} else {
    pointLayers <- geom_point(size=0.6, alpha=0.8)
}
p <- ggplot(plotDf, aes(x=BPCUM, y=LOGP, colour=CHR)) +
    pointLayers +
    scale_colour_manual(values=rep(c("#1F4E79", "#7FA7D1"), length.out=length(chromOrder)), guide="none") +
    scale_x_continuous(breaks=axisDf$CENTER, labels=axisDf$CHR, expand=c(0.01, 0)) +
    scale_y_continuous(expand=expansion(mult=c(0, 0.08))) +
    labs(x="Chromosome", y=expression(-log[10](italic(P))),
         title=opt$title,
         subtitle=sprintf("n variants = %s; lambda GC = %.3f", format(nTests, big.mark=","), lambdaGC)) +
    coord_cartesian(clip="off") +
    theme_bw() +
    theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          plot.margin=margin(5.5, 40, 5.5, 5.5))

if (!is.na(sigLine)) p <- p + geom_hline(yintercept=-log10(sigLine), colour="#C0392B", linetype="solid")
if (!is.na(meffLine)) p <- p + geom_hline(yintercept=-log10(meffLine), colour="#C0392B", linetype="dashed")
if (!is.na(suggestiveLine)) p <- p + geom_hline(yintercept=-log10(suggestiveLine), colour="grey40", linetype="dotted")

if (!is.na(opt$annotate)) {
    labelDf <- gwas[P < opt$annotate][order(P)][, .SD[1], by=CHR]
    if (nrow(labelDf) > 0) {
        p <- p + geom_text(data=labelDf, aes(label=ID), colour="black", size=2.5, vjust=-0.6)
    }
}

ggsave(paste0(opt$output, ".manhattan.", opt$format), p, width=opt$width, height=opt$height, dpi=300)

###########################################################################
##                                                                       ##
##                             QQ PLOT                                   ##
##                                                                       ##
###########################################################################

qqDf <- data.table(OBS=sort(gwas$LOGP, decreasing=TRUE),
                   EXP=-log10(ppoints(nTests)))
qqDf[, CI_UPPER := -log10(qbeta(0.025, seq_len(nTests), nTests - seq_len(nTests) + 1))]
qqDf[, CI_LOWER := -log10(qbeta(0.975, seq_len(nTests), nTests - seq_len(nTests) + 1))]
qqPlot <- rbind(qqDf[EXP >= 2], unique(qqDf[EXP < 2][, EBIN := round(EXP, 3)], by="EBIN")[, !"EBIN"])

q <- ggplot(qqPlot, aes(x=EXP, y=OBS)) +
    geom_ribbon(aes(ymin=CI_LOWER, ymax=CI_UPPER), fill="grey85") +
    geom_abline(intercept=0, slope=1, colour="#C0392B") +
    geom_point(size=0.6, colour="#1F4E79") +
    labs(x=expression(Expected~-log[10](italic(P))), y=expression(Observed~-log[10](italic(P))),
         title=opt$title, subtitle=sprintf("lambda GC = %.3f", lambdaGC)) +
    theme_bw()
ggsave(paste0(opt$output, ".qq.", opt$format), q, width=5, height=5, dpi=300)

###########################################################################
##                                                                       ##
##                            TOP HITS                                   ##
##                                                                       ##
###########################################################################

topHits <- gwas[order(P)][seq_len(min(opt$top_n, nrow(gwas))), .(ID, CHR, POS, P)]
fwrite(topHits, paste0(opt$output, ".top_hits.tsv"), sep="\t")

# Associated regions: in family-based populations, LD from few meioses spreads
# association over Mb-scale blocks, so hits are summarised as regions in which
# consecutive variants with P < region_p lie within region_gap bp of each other
hits <- gwas[P < opt$region_p][order(CHR, POS)]
if (nrow(hits) > 0) {
    hits[, REGION := cumsum(c(TRUE, diff(POS) > opt$region_gap | CHR[-1] != CHR[-.N]))]
    regions <- hits[, .(START=min(POS), END=max(POS), N_VARIANTS=.N,
                        LEAD_ID=ID[which.min(P)], LEAD_POS=POS[which.min(P)], LEAD_P=min(P),
                        N_SIG=if (is.na(sigLine)) NA_integer_ else sum(P < sigLine)),
                    by=.(CHR, REGION)][, REGION := NULL]
    regions[, SPAN_MB := round((END - START) / 1e6, 2)]
    fwrite(regions[order(LEAD_P)], paste0(opt$output, ".regions.tsv"), sep="\t")
    cat(sprintf("Associated regions (P < %g, gap <= %g bp): %d\n", opt$region_p, opt$region_gap, nrow(regions)))
}

cat("Program completed successfully!\n")
