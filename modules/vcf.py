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

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from vcftopia import VCFTopia
from parsing import read_gz_file, BgzCapableWriter, parse_2col_tsv_as_dict
from annotarium_importers import import_annotarium_gff3

def vcf_relabel(args):
    '''
    Handles "vcf relabel" mode of variantopia.
    '''
    if args.samplesMetadataTsv != None:
        sampleDict = parse_2col_tsv_as_dict(args.samplesMetadataTsv)
    if args.chromosomesMetadataTsv != None:
        chromosomesDict = parse_2col_tsv_as_dict(args.chromosomesMetadataTsv)
    
    warnedChromosomes = set()
    with read_gz_file(args.vcfFile) as fileIn, BgzCapableWriter(args.outputFileName) as fileOut:
        for line in fileIn:
            # Header handling contingent on --samples
            if line.startswith("#CHROM"):
                sl = line.rstrip().split("\t")
                samples = sl[9:]
                
                # Run conversion of oldid to newid based on metadata file
                newSamples = []
                for oldid in samples:
                    if args.samplesMetadataTsv != None:
                        if oldid in sampleDict: # replacement is specified
                            newid = sampleDict[oldid]
                        elif not args.allowNoSampleMatch: # replacement NOT specified, and not allowed
                            raise ValueError(f"VCF sample '{oldid}' has no match in your metadata file '{args.samplesMetadataTsv}'")
                        else: # replacement NOT specified, but is allowed
                            print(f"# --relaxSamples is allowing '{oldid}' to be unchanged in the resulting VCF header")
                            newid = oldid
                        
                        newSamples.append(newid)
                    else: # --samples unspecified
                        newSamples.append(oldid)
                
                # Format and write updated line
                newsl = sl[0:9] + newSamples
                fileOut.write("\t".join(newsl) + ("\r\n" if line.endswith("\r\n") else "\n"))
            # Header handling contingent on --chromosomes
            elif line.startswith("##contig=<ID="):
                commentDict = VCFTopia.comment_as_dict(line)
                oldid = commentDict["ID"]
                
                if args.chromosomesMetadataTsv != None:
                    if oldid in chromosomesDict: # replacement is specified
                        newid = chromosomesDict[oldid]
                    elif not args.allowNoChromosomeMatch: # replacement NOT specified, and not allowed
                        raise ValueError(f"VCF chromosome '{oldid}' has no match in your metadata file '{args.chromosomesMetadataTsv}'")
                    else: # replacement NOT specified, but is allowed
                        if not oldid in warnedChromosomes:
                            print(f"# --relaxChromosomes is allowing '{oldid}' to be unchanged in the resulting VCF body")
                            warnedChromosomes.add(oldid)
                        newid = oldid
                    
                    fileOut.write(line.replace(f"ID={oldid}", f"ID={newid}"))
                else: # --chromosomes unspecified
                    fileOut.write(line)
            # Header lines with no specific handling
            elif line.startswith("#"):
                fileOut.write(line)
            # Body handling contigent on --chromosomes
            else:
                sl = line.rstrip().split("\t")
                oldid = sl[0]
                
                if args.chromosomesMetadataTsv != None:
                    if oldid in chromosomesDict: # replacement is specified
                        sl[0] = chromosomesDict[oldid]
                    elif not args.allowNoChromosomeMatch: # replacement NOT specified, and not allowed
                        raise ValueError(f"VCF chromosome '{oldid}' has no match in your metadata file '{args.chromosomesMetadataTsv}'")
                    else: # replacement NOT specified, but is allowed
                        if not oldid in warnedChromosomes:
                            print(f"# --relaxChromosomes is allowing '{oldid}' to be unchanged in the resulting VCF body")
                            warnedChromosomes.add(oldid)
                        pass
                    
                    line = "\t".join(sl) + ("\r\n" if line.endswith("\r\n") else "\n")
                else: # --chromosomes unspecified
                    pass
                
                fileOut.write(line)

def vcf_stats(args):
    '''
    Handles "vcf stats" mode of variantopia.
    '''
    vcf = VCFTopia(args.vcfFile)
    genomeStatsDF, sampleStatDF = vcf.comprehensive_statistics()
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("# genome-level statistics\n")
        genomeStatsDF.to_csv(fileOut, sep="\t", na_rep=".")
        fileOut.write("##\n# sample-level statistics\n")
        sampleStatDF.to_csv(fileOut, sep="\t", na_rep=".")

def vcf_filter(args):
    '''
    Handles "vcf filter" mode of variantopia.
    '''
    vcf = VCFTopia(args.vcfFile)
    if args.gff3File != None:
        GFF3Feature, GFF3Tarium = import_annotarium_gff3(args.annotariumDir)
        gff3 = GFF3Tarium(args.gff3File)
        gff3.create_ncls_index(typeToIndex=args.featureType)
    else:
        gff3 = None
    
    with BgzCapableWriter(args.outputFileName) as fileOut:
        # Write header
        fileOut.write(vcf.vcf.raw_header)
        
        # Iterate through and filter variants
        for variant in vcf:
            # Filter based on GFF3
            if gff3 != None:
                matches = gff3.ncls_finder(variant.POS, variant.POS, "contig", variant.CHROM)
                if len(matches) == 0: # no match, no variant output
                    continue
            
            # Write to file if we passed all previous filters
            fileOut.write(str(variant))
