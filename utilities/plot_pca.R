#!/usr/bin/env Rscript
# plot_pca.R
# Diagnostics to help choose the number of PCs to use as GWAS covariates, from
# PLINK2 --pca output. Run by variantopia 'plink pca', or standalone.
#
# Example:
#   Rscript plot_pca.R -e x.pca.eigenvec -v x.pca.eigenval -f x.psam \
#       -l x.pca.eigenvec.var -o x.pca
#
# Requires: R packages optparse, data.table, ggplot2 (and png, for --format pdf)
#
# Outputs:
#   <prefix>.scree.png      eigenvalues / approx. % variance per PC
#   <prefix>.pcs.png        PC pairs (1v2, 3v4, ...) coloured by family, shaped by phenotype
#   <prefix>.loadings.png   |SNP loading| along the genome per PC (if -l given)
#   (plots are .pdf rather than .png with --format pdf)
#   <prefix>.summary.tsv    per PC: eigenvalue, % variance, R^2 explained by family,
#                           phenotype association P, and loading concentration

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
})

optionList <- list(
    make_option(c("-e", "--eigenvec"), type="character", help="PLINK2 .eigenvec file"),
    make_option(c("-v", "--eigenval"), type="character", help="PLINK2 .eigenval file"),
    make_option(c("-f", "--fam"), type="character",
                help=".fam or .psam file for the PCA samples (family ID and phenotype)"),
    make_option(c("-l", "--loadings"), type="character", default=NA,
                help="Optional PLINK2 .eigenvec.var file (from --pca biallelic-var-wts)"),
    make_option(c("-o", "--output"), type="character", help="Output file prefix"),
    make_option("--plot_pcs", type="integer", default=8,
                help="Number of PCs to show in the pairs and loadings plots [default %default]"),
    make_option("--format", type="character", default="png",
                help="Plot file format; 'png' or 'pdf' (vector graphics, e.g., for editing) [default %default]"),
    make_option("--raster_dpi", type="integer", default=300,
                help=paste("With --format pdf, the resolution at which dense point layers are embedded",
                           "as images, keeping the PDF small enough to edit; 0 keeps them as vector",
                           "points [default %default]")),
    make_option("--title", type="character", default="", help="Plot title")
)
opt <- parse_args(OptionParser(option_list=optionList))
for (req in c("eigenvec", "eigenval", "fam", "output")) {
    if (is.null(opt[[req]])) stop(sprintf("--%s is required", req))
}
if (!opt$format %in% c("png", "pdf")) stop("--format must be 'png' or 'pdf'")
rasterise <- opt$format == "pdf" && opt$raster_dpi > 0
if (rasterise && !requireNamespace("png", quietly=TRUE)) {
    stop("--format pdf needs the R package 'png' to embed dense points as images; install it, or use --raster_dpi 0")
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

eigenvec <- fread(opt$eigenvec)
setnames(eigenvec, 1, sub("^#", "", names(eigenvec)[1]))
pcCols <- grep("^PC[0-9]+$", names(eigenvec), value=TRUE)
numPCs <- length(pcCols)
eigenval <- fread(opt$eigenval, header=FALSE)$V1

# A .psam has a '#'-prefixed header naming its columns, of which FID and PHENO1
# are optional; a .fam has six fixed columns and no header
if (startsWith(readLines(opt$fam, n=1), "#")) {
    fam <- fread(opt$fam, colClasses="character")
    setnames(fam, 1, sub("^#", "", names(fam)[1]))
    fam[, FAMID := if ("FID" %in% names(fam)) FID else "none"]
    fam[, PHENO := if ("PHENO1" %in% names(fam)) suppressWarnings(as.integer(PHENO1)) else NA_integer_]
} else {
    fam <- fread(opt$fam, header=FALSE, col.names=c("FAMID", "IID", "PAT", "MAT", "SEX", "PHENO"),
                 colClasses=list(character=1:4))
}
pcaDf <- merge(eigenvec[, c("IID", pcCols), with=FALSE], fam[, .(IID, FAMID, PHENO)], by="IID")
if (nrow(pcaDf) != nrow(eigenvec)) warning("Some PCA samples were not found in the .fam file")
pcaDf[, PHENOTYPE := factor(fcase(PHENO == 2, "case", PHENO == 1, "control", default="missing"),
                            levels=c("case", "control", "missing"))]

# Approximate proportion of total genetic variance: the trace of the
# standardised GRM is ~N, so eigenvalue / N
numSamples <- nrow(pcaDf)
summaryDf <- data.table(PC=seq_len(numPCs), EIGENVAL=eigenval[seq_len(numPCs)])
summaryDf[, PCT_VAR := 100 * EIGENVAL / numSamples]

###########################################################################
##                                                                       ##
##                    PER-PC STRUCTURE DIAGNOSTICS                       ##
##                                                                       ##
###########################################################################

# R^2 of each PC explained by family membership, and the association of each
# PC with phenotype (logistic regression Wald P) among phenotyped samples
hasFamilies <- length(unique(pcaDf$FAMID)) > 1
phenoDf <- pcaDf[PHENOTYPE != "missing"]
phenoDf[, Y := as.integer(PHENOTYPE == "case")]
summaryDf[, FAMILY_R2 := sapply(pcCols, function(pc) {
    if (!hasFamilies) return(NA_real_)
    summary(lm(pcaDf[[pc]] ~ factor(pcaDf$FAMID)))$r.squared
})]
summaryDf[, PHENO_P := sapply(pcCols, function(pc) {
    if (length(unique(phenoDf$Y)) < 2) return(NA_real_)
    coef(summary(glm(phenoDf$Y ~ phenoDf[[pc]], family=binomial)))[2, 4]
})]

###########################################################################
##                                                                       ##
##                         LOADING DIAGNOSTICS                           ##
##                                                                       ##
###########################################################################

# A PC driven by genome-wide structure has loadings spread across chromosomes;
# one driven by a single LD block / inversion concentrates them in a region.
# MAX_CHR_SHARE: largest share of squared loadings on one chromosome.
# TOP1PCT_SHARE: share of squared loadings carried by the top 1% of SNPs.
if (!is.na(opt$loadings)) {
    loadings <- fread(opt$loadings)
    setnames(loadings, 1, sub("^#", "", names(loadings)[1]))
    loadPCs <- intersect(pcCols, names(loadings))
    # biallelic-var-wts gives one row per variant; keep the first if duplicated
    loadings <- unique(loadings, by="ID")
    loadings[, POS := as.integer(sub("^[^:]+:([0-9]+):.*$", "\\1", ID))]

    concentration <- rbindlist(lapply(seq_along(loadPCs), function(i) {
        sq <- loadings[[loadPCs[i]]]^2
        chrShare <- tapply(sq, loadings$CHROM, sum) / sum(sq)
        topN <- max(1, ceiling(length(sq) * 0.01))
        data.table(PC=i,
                   MAX_CHR=names(chrShare)[which.max(chrShare)],
                   MAX_CHR_SHARE=max(chrShare),
                   TOP1PCT_SHARE=sum(sort(sq, decreasing=TRUE)[seq_len(topN)]) / sum(sq))
    }))
    summaryDf <- merge(summaryDf, concentration, by="PC", all.x=TRUE)

    # Plot |loading| along the genome for the first PCs
    plotPCs <- loadPCs[seq_len(min(opt$plot_pcs, length(loadPCs)))]
    chromKeys <- unique(loadings$CHROM)
    chromOrder <- chromKeys[order(suppressWarnings(as.numeric(gsub("\\D", "", chromKeys))), chromKeys)]
    loadings[, CHROM := factor(CHROM, levels=chromOrder)]
    setorder(loadings, CHROM, POS)
    chromLengths <- loadings[, .(LEN=max(POS)), by=CHROM]
    chromLengths[, OFFSET := shift(cumsum(as.numeric(LEN)), fill=0)]
    loadings <- merge(loadings, chromLengths[, .(CHROM, OFFSET)], by="CHROM", sort=FALSE)
    loadings[, BPCUM := POS + OFFSET]
    axisDf <- loadings[, .(CENTER=(min(BPCUM) + max(BPCUM)) / 2), by=CHROM]

    longLoad <- melt(loadings[, c("CHROM", "BPCUM", plotPCs), with=FALSE],
                     id.vars=c("CHROM", "BPCUM"), variable.name="PC", value.name="LOADING")
    longLoad[, ABSLOAD := abs(LOADING)]
    longLoad[, PC := factor(PC, levels=plotPCs)]
    lp <- ggplot(longLoad, aes(x=BPCUM, y=ABSLOAD, colour=CHROM)) +
        (if (rasterise) rasterised_points(size=0.2, alpha=0.5, dpi=opt$raster_dpi) else geom_point(size=0.2, alpha=0.5)) +
        facet_wrap(~PC, ncol=1, scales="free_y", strip.position="right") +
        scale_colour_manual(values=rep(c("#1F4E79", "#7FA7D1"), length.out=length(chromOrder)), guide="none") +
        scale_x_continuous(breaks=axisDf$CENTER, labels=axisDf$CHROM, expand=c(0.01, 0)) +
        labs(x="Chromosome", y="|SNP loading|", title=opt$title) +
        theme_bw() +
        theme(panel.grid=element_blank())
    ggsave(paste0(opt$output, ".loadings.", opt$format), lp, width=12, height=1.4 * length(plotPCs) + 1, dpi=200)
}

fwrite(summaryDf, paste0(opt$output, ".summary.tsv"), sep="\t")
print(summaryDf, digits=3)

###########################################################################
##                                                                       ##
##                               PLOTS                                   ##
##                                                                       ##
###########################################################################

# Scree plot
sp <- ggplot(summaryDf, aes(x=PC, y=PCT_VAR)) +
    geom_line(colour="grey50") +
    geom_point(colour="#1F4E79", size=2) +
    scale_x_continuous(breaks=seq_len(numPCs)) +
    labs(x="Principal component", y="Approx. % of genetic variance (eigenvalue / N)", title=opt$title) +
    theme_bw() +
    theme(panel.grid.minor=element_blank())
ggsave(paste0(opt$output, ".scree.", opt$format), sp, width=7, height=4, dpi=200)

# PC pairs coloured by family and shaped by phenotype
numPairs <- min(opt$plot_pcs, numPCs) %/% 2
pairsDf <- rbindlist(lapply(seq_len(numPairs), function(i) {
    x <- 2 * i - 1; y <- 2 * i
    data.table(FACET=sprintf("PC%d (%.1f%%) vs PC%d (%.1f%%)", x, summaryDf$PCT_VAR[x], y, summaryDf$PCT_VAR[y]),
               X=pcaDf[[pcCols[x]]], Y=pcaDf[[pcCols[y]]],
               FAMID=pcaDf$FAMID, PHENOTYPE=pcaDf$PHENOTYPE)
}))
pairsDf[, FACET := factor(FACET, levels=unique(FACET))]
pp <- ggplot(pairsDf, aes(x=X, y=Y, colour=FAMID, shape=PHENOTYPE)) +
    geom_point(size=1.8, alpha=0.8) +
    facet_wrap(~FACET, scales="free", ncol=2) +
    scale_shape_manual(values=c(case=17, control=16, missing=4)) + # legend lists only phenotypes present
    labs(x=NULL, y=NULL, colour="Family ID", shape="Phenotype", title=opt$title) +
    theme_bw() +
    theme(panel.grid.minor=element_blank())
ggsave(paste0(opt$output, ".pcs.", opt$format), pp, width=10, height=4 * ceiling(numPairs / 2) + 0.5, dpi=200)

cat("Program completed successfully!\n")
