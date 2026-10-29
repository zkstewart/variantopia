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
from parsing import parse_2col_tsv_as_dict

def parse_lead_variants(regionsFile=None, clumpsFile=None):
    '''
    Collects the lead variants of associated regions (a .regions.tsv from
    plot_manhattan.R) and the index variants of LD clumps (a plink2 .clumps file).

    Parameters:
        regionsFile -- (OPTIONAL) a string pointing to a .regions.tsv file
        clumpsFile -- (OPTIONAL) a string pointing to a .clumps file
    Returns:
        leads -- a list of dictionaries sorted by P, each structured like:
                 {
                     "id": variantID, # string
                     "chrom": chrom, # string
                     "pos": pos, # int
                     "p": p, # string, as written in the source file
                     "sources": ["region", "clump"] # either or both
                 }
    '''
    leadDict = {}
    for sourceFile, source, columns in [[regionsFile, "region", ["LEAD_ID", "CHR", "LEAD_POS", "LEAD_P"]],
                                        [clumpsFile, "clump", ["ID", "CHROM", "POS", "P"]]]:
        if sourceFile is None:
            continue
        with open(sourceFile, "r") as fileIn:
            header = fileIn.readline().rstrip("\r\n").lstrip("#").split("\t")
            indices = [ header.index(c) for c in columns ]
            for line in fileIn:
                sl = line.rstrip("\r\n").split("\t")
                variantID, chrom, pos, p = [ sl[i] for i in indices ]
                lead = leadDict.setdefault(variantID, { "id": variantID, "chrom": chrom, "pos": int(pos),
                                                        "p": p, "sources": [] })
                lead["sources"].append(source)
    return sorted(leadDict.values(), key=lambda lead: (float(lead["p"]), lead["id"]))

class LeadGeneReporter:
    '''
    Relates GWAS lead variants to the genes within a radius of them, using a GFF3
    of gene models and, optionally, an annotable (an annotarium .annotation.tsv
    of best BLAST hits per transcript) to name the genes. Each gene is
    represented by its longest isoform.

    The GFF3 and annotable are given as annotarium GFF3Tarium and Annotable
    objects (see annotarium_importers.py).

    Example:
        reporter = LeadGeneReporter(GFF3Tarium("genes.gff3"), Annotable("genes.annotation.tsv"))
        reporter.chromMap = "chrom_map.tsv"
        reporter.radius = 50000
        reporter.write(parse_lead_variants("x.regions.tsv", "x.clumps"), "x.lead_genes.tsv")
    '''
    HEADER = ["LEAD_ID", "CHR", "POS", "P", "SOURCE", "GENE_ID", "TRANSCRIPT_ID", "SEQID", "STRAND",
              "GENE_START", "GENE_END", "DISTANCE", "LOCATION", "NEAREST",
              "BEST_HIT", "GENE_NAME", "TAXON", "PCT_ID", "EVALUE", "GO"]
    ANNOTABLE_HIT_COLUMNS = ["Target_accession", "Gene_names", "Taxa_names", "Percentage_identity", "Expect_value"]
    ANNOTABLE_GO_COLUMN = "Best_mapped_GOs"

    def __init__(self, gff3Obj, annotable=None):
        self.gff3Obj = gff3Obj
        self.annotable = annotable
        self.chromMap = None
        self.radius = 50000
        self.gff3Obj.create_ncls_index(typeToIndex=["gene"])

    @property
    def chromMap(self):
        return self._chromMap

    @chromMap.setter
    def chromMap(self, value):
        '''
        None (GWAS chromosomes are GFF3 seqids), or a two-column TSV file in the
        format of 'vcf to plink --chromMap', relating GFF3 seqids (left) to GWAS
        chromosomes (right).
        '''
        if value is None:
            self._chromMap = None
        elif isinstance(value, str):
            if not os.path.isfile(value):
                raise FileNotFoundError(f"Chromosome map file '{value}' is not a file")
            self._chromMap = parse_2col_tsv_as_dict(value, keyIsLeft=False)
        else:
            raise TypeError("chromMap must be None or a file path string")

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        '''
        Distance (bp) either side of a lead variant within which genes are reported.
        '''
        if not isinstance(value, int) or value < 0:
            raise ValueError("radius must be an integer >= 0")
        self._radius = value

    def seqid(self, chrom):
        '''
        Returns the GFF3 seqid of a GWAS chromosome.
        '''
        if self.chromMap is None:
            return chrom
        if chrom not in self.chromMap:
            raise ValueError(f"GWAS chromosome '{chrom}' is not in the chromosome map")
        return self.chromMap[chrom]

    @staticmethod
    def position_in_gene(pos, gene, transcript):
        '''
        Locates a variant relative to a gene.

        Parameters:
            pos -- an integer of the variant position
            gene -- a GFF3Feature of the gene
            transcript -- a GFF3Feature of its longest isoform, or None
        Returns:
            distance -- an integer of the distance (bp) to the gene, being 0 within
                        the gene, negative upstream of it and positive downstream
                        (relative to the gene's strand)
            location -- a string of 'exon', 'intron', 'genic' (within a gene with no
                        exon models), 'upstream' or 'downstream'
        '''
        if gene.start <= pos <= gene.end:
            exons = getattr(transcript, "exon", getattr(transcript, "CDS", [])) if transcript is not None else []
            if exons == []:
                return 0, "genic"
            return 0, "exon" if any(e.start <= pos <= e.end for e in exons) else "intron"
        offset = pos - gene.start if pos < gene.start else pos - gene.end
        distance = -offset if gene.strand == "-" else offset
        return distance, "upstream" if distance < 0 else "downstream"

    def annotation(self, gene, transcript):
        '''
        Returns the annotable's first BLAST hit values and GO terms for a gene,
        looked up by its longest isoform's ID (then its own ID), or '.' values.
        '''
        empty = ["."] * (len(LeadGeneReporter.ANNOTABLE_HIT_COLUMNS) + 1)
        if self.annotable is None:
            return empty
        for queryID in ([transcript.ID] if transcript is not None else []) + [gene.ID]:
            if queryID in self.annotable:
                row = self.annotable[queryID]
                break
        else:
            return empty

        # Hits are formatted like 'hit1 [hit2][hit3]'; accessions contain no brackets, so give the hit count
        numHits = row["Target_accession"].count("[") + 1 if row["Target_accession"] not in {".", ""} else 0
        values = []
        for column in LeadGeneReporter.ANNOTABLE_HIT_COLUMNS:
            hits = self.annotable.split_hits(row[column], numHits) if column in row.index else None
            values.append(hits[0] if hits else ".")
        go = row[LeadGeneReporter.ANNOTABLE_GO_COLUMN] if LeadGeneReporter.ANNOTABLE_GO_COLUMN in row.index else "."
        return values + [go if go != "" else "."]

    def gene_rows(self, lead):
        '''
        Returns the output rows of the genes within radius of a lead variant,
        nearest first, or a single row of '.' gene values if there are none.
        '''
        seqid = self.seqid(lead["chrom"])
        leadValues = [lead["id"], lead["chrom"], str(lead["pos"]), lead["p"], "+".join(lead["sources"])]
        genes = self.gff3Obj.ncls_finder(max(1, lead["pos"] - self.radius), lead["pos"] + self.radius, "contig", seqid)
        if genes == []:
            return [ leadValues + [".", ".", seqid] + ["."] * (len(LeadGeneReporter.HEADER) - 8) ]

        located = []
        for gene in genes:
            transcript = self.gff3Obj.longest_feature(gene)
            distance, location = LeadGeneReporter.position_in_gene(lead["pos"], gene, transcript)
            located.append([gene, transcript, distance, location])
        located.sort(key=lambda x: (abs(x[2]), x[0].start, x[0].ID))
        nearest = abs(located[0][2])

        rows = []
        for gene, transcript, distance, location in located:
            rows.append(leadValues + [gene.ID, transcript.ID if transcript is not None else ".", seqid,
                                      gene.strand, str(gene.start), str(gene.end), str(distance), location,
                                      "yes" if abs(distance) == nearest else "no"] +
                        self.annotation(gene, transcript))
        return rows

    def write(self, leads, outputFileName):
        '''
        Writes the genes within radius of each lead variant.

        Parameters:
            leads -- a list of lead variant dictionaries (see parse_lead_variants())
            outputFileName -- a string of the location to write the TSV file
        Returns:
            counts -- a dictionary with structure like:
                      {
                          "leads": int, # lead variants written
                          "withGenes": int, # lead variants with a gene within radius
                          "genes": int, # lead variant and gene rows written
                          "missingSeqids": [seqid, ...] # seqids absent from the GFF3
                      }
        '''
        counts = { "leads": 0, "withGenes": 0, "genes": 0, "missingSeqids": [] }
        with open(outputFileName, "w") as fileOut:
            fileOut.write("\t".join(LeadGeneReporter.HEADER) + "\n")
            for lead in leads:
                rows = self.gene_rows(lead)
                seqid = rows[0][7]
                if seqid not in self.gff3Obj.contigs and seqid not in counts["missingSeqids"]:
                    counts["missingSeqids"].append(seqid)
                counts["leads"] += 1
                if rows[0][5] != ".":
                    counts["withGenes"] += 1
                    counts["genes"] += len(rows)
                for row in rows:
                    fileOut.write("\t".join(row) + "\n")
        return counts
