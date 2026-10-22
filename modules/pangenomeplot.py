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

import os, sys
from matplotlib.patches import Patch
from pycirclize import Circos

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from vcftopia import VCFTopia

# Colours ordered as parent 1 haplotypes 1 and 2, then parent 2 haplotypes 1 and 2
PAINT_COLOURS = ["#1f5fa8", "#8ec1ee", "#c75a00", "#fdbf80"]
UNPAINTED_COLOUR = "#e6e6e6"

def parse_paint_file(paintFile):
    '''
    Parses the TSV file produced by "pangenome paint" mode.
    
    Parameters:
        paintFile -- a string indicating the location of the paint TSV file
    Returns:
        paintDict -- a dictionary with structure like:
                     {
                         sampleID: {
                             chrom: {
                                 parentID: [[start, end, hap], ...],
                                 ...
                             },
                             ...
                         },
                         ...
                     }
                     Where hap is the integer (1 or 2) haplotype of that parent
        parents -- a list of the two parent IDs in their order of appearance
        chroms -- a set of the chromosome IDs that have blocks
    '''
    HEADER = ["sample", "chrom", "start", "end", "parent", "haplotype", "supporting", "conflicting"]
    
    paintDict = {}
    parents = []
    chroms = set()
    with open(paintFile, "r") as fileIn:
        try:
            header = fileIn.readline().rstrip("\r\n").split("\t")
        except UnicodeDecodeError:
            raise ValueError(f"Paint file '{paintFile}' is not a plain text file; is it the output of 'pangenome paint' mode?")
        if header != HEADER:
            raise ValueError(f"Paint file '{paintFile}' does not have the expected header '{HEADER}'; " +
                             "is it the output of 'pangenome paint' mode?")
        
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if sl == [""]:
                continue
            sampleID, chrom, start, end, parentID, hap = sl[:6]
            
            if not parentID in parents:
                parents.append(parentID)
            chroms.add(chrom)
            paintDict.setdefault(sampleID, {})
            paintDict[sampleID].setdefault(chrom, {})
            paintDict[sampleID][chrom].setdefault(parentID, [])
            paintDict[sampleID][chrom][parentID].append([int(start), int(end), int(hap)])
    
    if len(parents) != 2:
        raise ValueError(f"Paint file '{paintFile}' should contain exactly two parents, but found '{parents}'")
    return paintDict, parents, chroms

def get_chromosome_lengths(vcfFile, chroms, paintDict):
    '''
    Obtains the lengths of painted chromosomes from the VCF header, ordered as they occur
    in the header. If the header lacks lengths, the end of the last painted block is used.
    
    Parameters:
        vcfFile -- a string indicating the location of the VCF file
        chroms -- a set of chromosome IDs to obtain lengths for
        paintDict -- a dictionary as returned by parse_paint_file()
    Returns:
        chromLengths -- a dictionary with structure like:
                        {
                            chrom: length,
                            ...
                        }
    '''
    vcf = VCFTopia(vcfFile)
    lengths = vcf.lengths
    
    missingChroms = chroms.difference(vcf.contigs)
    if len(missingChroms) > 0:
        raise ValueError(f"Paint file chromosomes '{missingChroms}' do not occur within the VCF file contigs")
    
    chromLengths = {}
    for chrom in vcf.contigs:
        if not chrom in chroms:
            continue
        if lengths is not None and lengths.get(chrom, 0) > 0:
            chromLengths[chrom] = lengths[chrom]
        else:
            print(f"WARNING: VCF header lacks a length for '{chrom}'; the end of its last painted block will be used")
            chromLengths[chrom] = max([ end for chromDict in paintDict.values()
                                        for blocks in chromDict.get(chrom, {}).values()
                                        for _, end, _ in blocks ])
    return chromLengths

def pan_plot_paint(args):
    paintDict, parents, chroms = parse_paint_file(args.paintFile)
    chromLengths = get_chromosome_lengths(args.vcfFile, chroms, paintDict)
    
    # Determine which samples to plot and in what order
    if args.samples == []:
        samples = list(paintDict.keys())
    else:
        missingSamples = [ sampleID for sampleID in args.samples if not sampleID in paintDict ]
        if len(missingSamples) > 0:
            raise ValueError(f"--samples '{missingSamples}' do not occur within the paint file samples '{list(paintDict.keys())}'")
        samples = args.samples
    
    # Derive the radial layout of tracks; each sample has an outer and inner track for its
    # haplotypes inherited from the first and second parent, respectively
    OUTER_RADIUS = 92
    INNER_RADIUS = 38
    PAIR_GAP_RATIO = 0.08 # gap between a sample's two tracks relative to track height
    SAMPLE_GAP_RATIO = 1.0 # gap between samples relative to track height
    trackHeight = (OUTER_RADIUS - INNER_RADIUS) / ((len(samples) * (2 + PAIR_GAP_RATIO)) + ((len(samples) - 1) * SAMPLE_GAP_RATIO))
    pairHeight = trackHeight * (2 + PAIR_GAP_RATIO)
    
    trackRadii = [] # list of [(p1Lower, p1Upper), (p2Lower, p2Upper)] for each sample
    currentPosition = OUTER_RADIUS
    for _ in samples:
        p1Radii = (currentPosition - trackHeight, currentPosition)
        p2Radii = (currentPosition - pairHeight, currentPosition - pairHeight + trackHeight)
        trackRadii.append([p1Radii, p2Radii])
        currentPosition -= pairHeight + (trackHeight * SAMPLE_GAP_RATIO)
    
    # Initialise circos instance, leaving a gap at the top for sample labels
    LABEL_GAP = 24
    circos = Circos(chromLengths, space=0 if len(chromLengths) == 1 else 2,
                    start=LABEL_GAP / 2, end=360 - (LABEL_GAP / 2), endspace=False)
    
    for sector in circos.sectors:
        # Plot outer chromosome track
        sector.text(sector.name.rsplit("#", 1)[-1], size=10) # drop any PanSN prefix
        outerTrack = sector.add_track((94.7, 95))
        outerTrack.axis(fc="black")
        minorInterval = args.tickInterval / 10
        if sector.size > minorInterval:
            # Omit labels near the sector end to avoid overlapping the next sector's first label
            outerTrack.xticks_by_interval(args.tickInterval,
                                          label_formatter=lambda v, size=sector.size: "" if size - v < args.tickInterval * 0.3 else f"{v / 1000000:.0f} Mb")
            outerTrack.xticks_by_interval(minorInterval, tick_length=1, show_label=False)
        
        # Plot each sample's inherited haplotype blocks
        for sampleID, radii in zip(samples, trackRadii):
            for parentNum, (parentID, rLim) in enumerate(zip(parents, radii)):
                track = sector.add_track(rLim)
                track.axis(fc=UNPAINTED_COLOUR, ec="none")
                for start, end, hap in paintDict[sampleID].get(sector.name, {}).get(parentID, []):
                    track.rect(start - 1, end, fc=args.colours[(parentNum * 2) + (hap - 1)], ec="none", lw=0) # 1-based to 0-based
    
    # Label samples within the top gap
    labelSize = min(8, max(3, pairHeight * 1.6)) if args.labelSize is None else args.labelSize
    for sampleID, (p1Radii, p2Radii) in zip(samples, trackRadii):
        circos.text(sampleID, r=(p1Radii[1] + p2Radii[0]) / 2, deg=0, size=labelSize, va="center")
    
    # Generate the plot with a legend at its centre
    fig = circos.plotfig(figsize=(args.figSize, args.figSize))
    handles = [ Patch(color=args.colours[(parentNum * 2) + hapIndex], label=f"{parentID} hap{hapIndex + 1}")
               for parentNum, parentID in enumerate(parents) for hapIndex in range(2) ]
    circos.ax.legend(handles=handles,
                     bbox_to_anchor=(0.5, 0.5),
                     loc="center",
                     ncols=2,
                     frameon=True, # pymsaviz disables legend frames globally on import
                     title=f"Outer bar: from {parents[0]}\nInner bar: from {parents[1]}")
    
    fig.savefig(args.outputFileName, dpi=300)
