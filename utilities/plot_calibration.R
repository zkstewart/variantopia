#!/usr/bin/env Rscript
# plot_calibration.R
# Genomic inflation (lambda GC) and the GRM variance component (tau) of SAIGE
# null models fitted with each candidate set of PC covariates, from the table
# written by variantopia 'saige calibrate'. Lambda should approach 1 once
# confounding structure is adequately modelled; the suggested set (lambda
# closest to 1) is highlighted. Run by variantopia 'saige calibrate', or
# standalone.
#
# Example:
#   Rscript plot_calibration.R -i x.pc_calibration.tsv -o x.pc_calibration
#
# Outputs:
#   <prefix>.png    lambda GC and tau for each PC set (.pdf with --format pdf)
#
# Requires: R packages optparse, data.table, ggplot2

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(ggplot2)
})

optionList <- list(
    make_option(c("-i", "--input"), type="character", help="PC calibration TSV (from variantopia 'saige calibrate')"),
    make_option(c("-o", "--output"), type="character", help="Output file prefix"),
    make_option("--format", type="character", default="png",
                help="Plot file format; 'png' or 'pdf' (vector graphics, e.g., for editing) [default %default]"),
    make_option("--title", type="character", default="", help="Plot title")
)
opt <- parse_args(OptionParser(option_list=optionList))
for (req in c("input", "output")) {
    if (is.null(opt[[req]])) stop(sprintf("--%s is required", req))
}
if (!opt$format %in% c("png", "pdf")) stop("--format must be 'png' or 'pdf'")

###########################################################################
##                                                                       ##
##                                 PLOT                                  ##
##                                                                       ##
###########################################################################

# One panel per statistic, with PC sets in table order
calib <- fread(opt$input, colClasses=list(character=c("TAG", "COVARIATES", "SUGGESTED")))
calib[, TAG := factor(TAG, levels=TAG)]
long <- rbind(
    calib[, .(TAG, SUGGESTED, STAT="Lambda GC (pruned SNPs)", VALUE=LAMBDA_GC)],
    calib[, .(TAG, SUGGESTED, STAT="Tau (GRM variance)", VALUE=TAU)]
)
long[, STAT := factor(STAT, levels=unique(STAT))]
reference <- data.table(STAT=factor("Lambda GC (pruned SNPs)", levels=levels(long$STAT)), Y=1)

p <- ggplot(long, aes(x=TAG, y=VALUE, group=1)) +
    geom_hline(data=reference, aes(yintercept=Y), linetype="dashed", colour="grey50") +
    geom_line(colour="grey40") +
    geom_point(aes(colour=SUGGESTED == "yes"), size=2.5) +
    scale_colour_manual(values=c("FALSE"="grey40", "TRUE"="#D55E00"), labels=c("FALSE"="Other", "TRUE"="Suggested"),
                        name=NULL) +
    facet_wrap(~STAT, ncol=1, scales="free_y") +
    labs(title=opt$title, subtitle="Suggested set: lambda GC closest to 1", x="PC covariate set", y=NULL) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top")

ggsave(paste0(opt$output, ".", opt$format), p, width=max(5, 0.45 * nrow(calib) + 2), height=6, dpi=200, bg="white")
cat(sprintf("%d PC sets plotted; suggested set is %s\n", nrow(calib), as.character(calib[SUGGESTED == "yes", TAG])))
cat("Program completed successfully!\n")
