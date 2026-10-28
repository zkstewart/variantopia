#!/usr/bin/env Rscript
# plot_king.R
# Heatmap of PLINK2 KING-robust kinship between all samples, hierarchically
# clustered, with each sample's recorded parents alongside. If the recorded
# pedigree is correct, offspring should cluster with their siblings, and the
# parent symbols should form consistent blocks alongside those clusters. Run by
# variantopia 'plink king', or standalone.
#
# Example:
#   Rscript plot_king.R -k x.king.kin0 -f x.psam -c x.king.pedigree_check.tsv -o x.king
#
# Outputs:
#   <prefix>.heatmap.png    clustered kinship heatmap with recorded parents
#
# Requires: R packages optparse, data.table, ggplot2

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
})

optionList <- list(
    make_option(c("-k", "--kin0"), type="character", help="PLINK2 .kin0 file (--make-king-table)"),
    make_option(c("-f", "--fam"), type="character",
                help=".fam or .psam file recording each sample's parents (PAT/MAT)"),
    make_option(c("-c", "--check"), type="character", default=NA,
                help="Optional pedigree check TSV (from variantopia 'plink king'); failed pairs are marked"),
    make_option(c("-o", "--output"), type="character", help="Output file prefix"),
    make_option("--linkage", type="character", default="average",
                help="hclust linkage method [default %default]"),
    make_option("--title", type="character", default="", help="Plot title")
)
opt <- parse_args(OptionParser(option_list=optionList))
for (req in c("kin0", "fam", "output")) {
    if (is.null(opt[[req]])) stop(sprintf("--%s is required", req))
}

###########################################################################
##                                                                       ##
##                           DATA LOADING                                ##
##                                                                       ##
###########################################################################

# Kinship matrix; self-kinship is 0.5 by definition
kin <- fread(opt$kin0, colClasses=list(character=c("IID1", "IID2")))
setnames(kin, 1, sub("^#", "", names(kin)[1]))
samples <- sort(unique(c(kin$IID1, kin$IID2)))
numSamples <- length(samples)
kinMatrix <- matrix(NA_real_, numSamples, numSamples, dimnames=list(samples, samples))
kinMatrix[cbind(kin$IID1, kin$IID2)] <- kin$KINSHIP
kinMatrix[cbind(kin$IID2, kin$IID1)] <- kin$KINSHIP
diag(kinMatrix) <- 0.5
if (anyNA(kinMatrix)) stop("The .kin0 file does not contain every pair of samples")

# Recorded parents; a .psam has a '#'-prefixed header naming its columns, and
# a .fam has six fixed columns and no header
if (startsWith(readLines(opt$fam, n=1), "#")) {
    fam <- fread(opt$fam, colClasses="character")
    setnames(fam, 1, sub("^#", "", names(fam)[1]))
    for (col in c("PAT", "MAT")) if (!col %in% names(fam)) fam[, (col) := "0"]
} else {
    fam <- fread(opt$fam, header=FALSE, col.names=c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"),
                 colClasses="character")
}
fam <- fam[IID %in% samples]

# Failed parent-offspring pairs, if a pedigree check was given
failed <- data.table(IID=character(0), ROLE=character(0))
if (!is.na(opt$check)) {
    failed <- fread(opt$check, colClasses="character")[STATUS == "FAIL", .(IID, ROLE)]
}

###########################################################################
##                                                                       ##
##                            CLUSTERING                                 ##
##                                                                       ##
###########################################################################

# Distance is the shortfall from self-kinship, so closer relatives are nearer
hc <- hclust(as.dist(0.5 - kinMatrix), method=opt$linkage)
leafOrder <- samples[hc$order]
rowY <- setNames(rev(seq_len(numSamples)), leafOrder) # first leaf at the top

# Dendrogram segments: each merge joins two children at the merge height
leafY <- rowY[samples]
nodeY <- numeric(nrow(hc$merge))
childY <- function(k) if (k < 0) leafY[[-k]] else nodeY[k]
childH <- function(k) if (k < 0) 0 else hc$height[k]
segments <- rbindlist(lapply(seq_len(nrow(hc$merge)), function(i) {
    a <- hc$merge[i, 1]; b <- hc$merge[i, 2]
    ya <- childY(a); yb <- childY(b); h <- hc$height[i]
    nodeY[i] <<- (ya + yb) / 2
    data.table(H1=c(childH(a), childH(b), h), H2=c(h, h, h), Y1=c(ya, yb, ya), Y2=c(ya, yb, yb))
}))

###########################################################################
##                                                                       ##
##                              LAYOUT                                   ##
##                                                                       ##
###########################################################################

# x positions: dendrogram (negative), heatmap (1..n), parent annotation columns
# (each recorded parent column is followed by its pedigree check column, if a
# check was given), then sample labels
dendWidth <- 0.15 * numSamples
annotStep <- 0.05 * numSamples
if (!is.na(opt$check)) {
    annotOffsets <- c(PAT=1, PAT_FAIL=1.5, MAT=2.2, MAT_FAIL=2.7, SELF=3.4)
    annotLabels <- c(PAT="Father", PAT_FAIL="Failed", MAT="Mother", MAT_FAIL="Failed", SELF="As parent")
} else {
    annotOffsets <- c(PAT=1, MAT=1.7, SELF=2.4)
    annotLabels <- c(PAT="Father", MAT="Mother", SELF="As parent")
}
annotX <- numSamples + annotStep * annotOffsets
labelX <- max(annotX) + annotStep * 0.6
segments[, `:=`(X1=-H1 / max(hc$height) * dendWidth, X2=-H2 / max(hc$height) * dendWidth)]

heatDf <- data.table(IID1=rep(leafOrder, times=numSamples), IID2=rep(leafOrder, each=numSamples))
heatDf[, `:=`(X=match(IID2, leafOrder), Y=rowY[IID1], KINSHIP=kinMatrix[cbind(IID1, IID2)])]

# Recorded parents as annotation points; samples which are themselves recorded
# parents are marked with their own parent symbol
annotDf <- rbind(
    fam[PAT != "0", .(IID, ROLE="PAT", PARENT=PAT)],
    fam[MAT != "0", .(IID, ROLE="MAT", PARENT=MAT)]
)
parents <- sort(unique(annotDf$PARENT))
annotDf <- rbind(annotDf, data.table(IID=intersect(parents, samples), ROLE="SELF",
                                     PARENT=intersect(parents, samples)))
annotDf[, `:=`(X=annotX[ROLE], Y=rowY[IID])]
failedDf <- failed[, .(X=annotX[paste0(ROLE, "_FAIL")], Y=rowY[IID])]

# Each parent gets a unique colour and shape combination, so many parents can
# be told apart
palette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")
shapes <- c(16, 17, 15, 18, 8, 3, 1, 2)
parentColours <- setNames(palette[(seq_along(parents) - 1) %% length(palette) + 1], parents)
parentShapes <- setNames(shapes[((seq_along(parents) - 1) %/% length(palette)) %% length(shapes) + 1], parents)

# Sizes scale with the number of samples so rows stay legible
plotHeight <- max(6, 0.06 * numSamples + 2)
rowPt <- (plotHeight - 2) * 72 / numSamples
textSize <- min(3, 0.8 * rowPt / .pt)
pointSize <- min(3, 0.9 * rowPt / .pt)

# Room above the heatmap for the (vertical) column labels: ~0.075 inches per
# character at size 3, converted to rows
labelRoom <- (max(nchar(annotLabels)) * 0.075 + 0.1) * numSamples / (plotHeight - 1.5)

###########################################################################
##                                                                       ##
##                               PLOT                                    ##
##                                                                       ##
###########################################################################

# Kinship colours are pinned to the KING relationship degree thresholds across
# the full range of kinship (0 to 0.5 for self), interpolating between them, so
# colours mean the same relatedness in every plot; negative kinship (unrelated)
# is shown as 0. Pinned hues are from the colour-blind-friendly Okabe-Ito palette,
# chosen so neighbouring degrees differ in hue. Lightness mostly decreases with
# kinship (L* 100, 89, 58, 61, 46, 0); 2nd degree green is marginally darker than
# 1st degree pink, accepted for its stronger contrast with the other hues
kingBreaks <- c(0, 0.0442, 0.0884, 0.177, 0.354, 0.5)
kingColours <- c("#FFFFFF", "#F0E442", "#009E73", "#CC79A7", "#0072B2", "#000000")
kingLabels <- c("0 (unrelated)", "0.044 (3rd degree)", "0.088 (2nd degree)", "0.177 (1st degree)",
                "0.354 (duplicate)", "0.5 (self)")
p <- ggplot() +
    geom_raster(data=heatDf, aes(x=X, y=Y, fill=KINSHIP)) + # unlike geom_tile, no seams between cells
    geom_segment(data=segments, aes(x=X1, xend=X2, y=Y1, yend=Y2), linewidth=0.2, colour="grey30") +
    geom_point(data=annotDf, aes(x=X, y=Y, colour=PARENT, shape=PARENT), size=pointSize) +
    geom_point(data=failedDf, aes(x=X, y=Y), shape=4, colour="red", size=pointSize, stroke=0.8) +
    annotate("text", x=annotX, y=numSamples + 1, label=annotLabels[names(annotX)], angle=90, hjust=0, size=3) +
    annotate("text", x=labelX, y=rowY, label=names(rowY), hjust=0, size=textSize) +
    scale_fill_gradientn(colours=kingColours, values=scales::rescale(kingBreaks), limits=range(kingBreaks),
                         breaks=kingBreaks, labels=kingLabels, oob=scales::squish, name="Kinship",
                         guide=guide_colourbar(barheight=unit(16, "lines"))) +
    scale_colour_manual(values=parentColours, name="Recorded parent") +
    scale_shape_manual(values=parentShapes, name="Recorded parent") +
    guides(colour=guide_legend(ncol=ceiling(length(parents) / 30)),
           shape=guide_legend(ncol=ceiling(length(parents) / 30))) +
    scale_x_continuous(expand=expansion(add=c(1, annotStep * 2.5))) +
    scale_y_continuous(expand=expansion(add=c(1, labelRoom))) +
    labs(title=opt$title,
         subtitle=if (nrow(failed) > 0) "Red \u00d7 (Failed): recorded parent below the pedigree check's kinship threshold" else NULL) +
    theme_void() +
    theme(plot.margin=margin(10, 10, 10, 10))

ggsave(paste0(opt$output, ".heatmap.png"), p, width=plotHeight * 1.35 + 2, height=plotHeight, dpi=200,
       limitsize=FALSE, bg="white")
cat(sprintf("%d samples; %d recorded parents; %d failed parent-offspring pairs marked\n",
            numSamples, length(parents), nrow(failedDf)))
cat("Program completed successfully!\n")
