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

import os, sys, bisect
from collections import Counter
import numpy as np

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

def parse_saige_leads(resultsFile, leads, windowKb, suggestive, traitColumns=[]):
    '''
    Reads the SAIGE results of lead variants, and counts the tested and
    suggestive variants near each.

    Parameters:
        resultsFile -- a string pointing to a SAIGE single-variant results file
        leads -- a list of lead variant dictionaries (see parse_lead_variants())
        windowKb -- a number of the distance (kb) either side of a lead variant
                    within which variants are counted
        suggestive -- a float of the P threshold of suggestive variants
        traitColumns -- (OPTIONAL) a list of trait-specific result columns to read
                        (e.g., BinaryTrait.SAIGE_RESULT_COLUMNS)
    Returns:
        leadStats -- a dictionary with structure like:
                     {
                         variantID: {
                             "Allele1": str, "Allele2": str, "AF_Allele2": str,
                             "BETA": str, "SE": str, # and traitColumns
                             "nWindow": int, # tested variants within windowKb
                             "nWindowSuggestive": int # those with P < suggestive
                         },
                         ...
                     }
    '''
    STAT_COLUMNS = ["Allele1", "Allele2", "AF_Allele2", "BETA", "SE"] + traitColumns
    window = int(windowKb * 1000)
    leadPositions = {} # chrom: [sorted positions, matching lead IDs]
    for lead in sorted(leads, key=lambda lead: lead["pos"]):
        positions, ids = leadPositions.setdefault(lead["chrom"], [[], []])
        positions.append(lead["pos"])
        ids.append(lead["id"])
    leadStats = { lead["id"]: { "nWindow": 0, "nWindowSuggestive": 0 } for lead in leads }

    with open(resultsFile, "r") as fileIn:
        header = fileIn.readline().rstrip("\r\n").split("\t")
        chromIndex, posIndex, idIndex, pIndex = [ header.index(c) for c in ["CHR", "POS", "MarkerID", "p.value"] ]
        statIndices = [ header.index(c) for c in STAT_COLUMNS ]
        for line in fileIn:
            sl = line.split("\t", max(chromIndex, posIndex) + 1)
            if sl[chromIndex] not in leadPositions:
                continue
            positions, ids = leadPositions[sl[chromIndex]]
            pos = int(sl[posIndex])
            start, end = bisect.bisect_left(positions, pos - window), bisect.bisect_right(positions, pos + window)
            if start == end:
                continue
            sl = line.rstrip("\r\n").split("\t")
            isSuggestive = sl[pIndex] != "NA" and float(sl[pIndex]) < suggestive
            for variantID in ids[start:end]:
                leadStats[variantID]["nWindow"] += 1
                leadStats[variantID]["nWindowSuggestive"] += int(isSuggestive)
            if sl[idIndex] in leadStats:
                leadStats[sl[idIndex]].update({ c: sl[i] for c, i in zip(STAT_COLUMNS, statIndices) })

    missing = [ variantID for variantID, stats in leadStats.items() if "Allele2" not in stats ]
    if missing != []:
        raise ValueError(f"{len(missing)} lead variant(s) are not in '{resultsFile}' (e.g., '{missing[0]}')")
    return leadStats

def mean_impute(dosages):
    '''
    Returns a copy of a dosage matrix (samples x variants) with missing values
    replaced by each variant's mean dosage.
    '''
    means = np.nanmean(dosages, axis=0)
    return np.where(np.isnan(dosages), means, dosages)

def inter_chromosomal_r2(dosages, chroms):
    '''
    Computes the r^2 between all pairs of variants on different chromosomes, as
    a null distribution of LD between unlinked loci.

    Parameters:
        dosages -- a numpy array (samples x variants) of dosages, NaN if missing
        chroms -- a list of each variant's chromosome
    Returns:
        nullR2 -- a sorted numpy array of the r^2 values
    '''
    keep = np.nanstd(dosages, axis=0) > 0
    chroms = np.array(chroms)[keep]
    r = np.corrcoef(mean_impute(dosages[:, keep]), rowvar=False)
    upper = np.triu_indices(len(chroms), 1)
    inter = chroms[upper[0]] != chroms[upper[1]]
    return np.sort(r[upper][inter] ** 2)

class LeadInspector:
    '''
    Inspects GWAS lead variants for the direction of their effect, the families
    their effect allele segregates in, its parental origin (where parents are
    genotyped), Mendelian errors, and LD with stronger lead variants on other
    chromosomes. The latter flags candidate 'ghost' associations: variants that
    may only associate through LD with a stronger locus, such as reads from it
    mis-mapped to paralogous or rearranged sequence. LD is judged against a null
    distribution of r^2 between unlinked variants on different chromosomes: a
    lead is linked if the empirical P of its strongest r^2, Bonferroni-adjusted
    for the stronger leads it was compared with, is below alpha (or, if linkR2
    is set, if its r^2 is at least that). Links are followed to their root (the
    first lead in the chain that is not itself linked), and a lead is only
    flagged as linked if its root is on another chromosome; one whose chain
    leads back to its own chromosome is flagged 'own-chromosome', as long-range
    LD within chromosomes is expected in families.

    Effect allele dosages are given for all samples of the fileset (so that
    parents not in the analysis can be used), while statistics are computed from
    the tested samples only. Phenotype summaries (e.g., case rates or mean
    phenotypes of carriers) are provided by a trait handler (see traits.py).

    Example:
        inspector = LeadInspector(leads, leadStats, iids, dosages, sampleDict, phenos, BinaryTrait())
        inspector.set_null(nullR2, alpha=0.05)
        inspector.write("x.leads.tsv", "x.lead_families.tsv")
    '''
    LEAD_COLUMNS = ["LEAD_ID", "CHR", "POS", "P", "SOURCE", "ALLELE1", "ALLELE2", "AF_ALLELE2", "BETA", "SE"]
    LEAD_GENOTYPE_COLUMNS = ["N_WINDOW", "N_WINDOW_SUGGESTIVE", "CARRIER_FAMILIES", "ORIGIN", "MENDEL_ERRORS",
                             "MENDEL_TESTED", "LINK_STATUS", "LINKED_TO", "LINK_ROOT", "LINK_R", "LINK_R2",
                             "LINK_EMPIRICAL_P", "LINK_P_ADJUSTED", "CARRIER_OVERLAP", "EXPECTED_OVERLAP"]
    FAMILY_PARENT_COLUMNS = ["PAT", "MAT", "PAT_GENOTYPE", "MAT_GENOTYPE", "ORIGIN", "LINEAGE",
                             "MENDEL_ERRORS", "MENDEL_TESTED"]
    MISSING_PARENT = { "0", "", "NA", "." }
    MENDEL_RATE_FLAG = 0.05 # Mendelian error rate (errors / offspring tested) at which leads are counted

    def __init__(self, leads, leadStats, iids, dosages, sampleDict, phenos, trait):
        '''
        Parameters:
            leads -- a list of lead variant dictionaries (see parse_lead_variants())
            leadStats -- a dictionary of SAIGE results (see parse_saige_leads())
            iids -- a list of the IIDs of all fileset samples, in dosages row order
            dosages -- a numpy array (samples x leads, in leads order) of SAIGE
                       Allele2 (effect allele) dosages, NaN if missing
            sampleDict -- a dictionary of .psam rows (see plink.parse_sample_info())
            phenos -- a dictionary pairing tested IIDs to their phenotype as a
                      number (for a binary trait, 1 = case and 0 = control)
            trait -- a trait handler (see traits.py) of the phenotype
        '''
        self.leads = leads
        self.leadStats = leadStats
        self.sampleDict = sampleDict
        self.trait = trait
        self.alpha = None
        self.linkR2 = None
        self.nullR2 = None

        # Hard calls of all samples, and of the tested samples as a matrix
        self.genotypes = { iid: row for iid, row in zip(iids, np.round(dosages)) }
        self.tested = [ iid for iid in iids if iid in phenos ]
        if len(self.tested) != len(phenos):
            raise ValueError(f"{len(phenos) - len(self.tested)} tested sample(s) have no genotypes in the fileset")
        self.testedDosages = np.array([ self.genotypes[iid] for iid in self.tested ])
        self.testedImputed = mean_impute(self.testedDosages)
        self.phenotypes = np.array([ phenos[iid] for iid in self.tested ], dtype=float)
        self.families = self._group_families()

    def set_null(self, nullR2, alpha=0.05, linkR2=None):
        '''
        Sets the null distribution of r^2 between unlinked variants, and how a
        lead is called linked to a stronger lead: by its adjusted empirical P
        being below alpha, or (if linkR2 is given) by its r^2 being at least that.
        '''
        if len(nullR2) == 0:
            raise ValueError("The null r^2 distribution is empty")
        self.nullR2 = np.sort(nullR2)
        self.alpha = alpha
        self.linkR2 = linkR2

    def empirical_p(self, r2):
        '''
        Returns the empirical P of an r^2 under the null distribution.
        '''
        numAtLeast = len(self.nullR2) - np.searchsorted(self.nullR2, r2, side="left")
        return (numAtLeast + 1) / (len(self.nullR2) + 1)

    def _group_families(self):
        '''
        Groups tested samples by FID (or, without FIDs, by their parents), noting
        each family's parents as the most common parent pair of its members.

        Returns:
            families -- a list of dictionaries sorted by FID, each like:
                        {
                            "fid": str, "members": [testedIndex, ...],
                            "pat": str or None, "mat": str or None, # parent IIDs
                            "offspring": [testedIndex, ...] # members with that pair
                        }
        '''
        groups = {}
        for index, iid in enumerate(self.tested):
            row = self.sampleDict[iid]
            pair = (row.get("PAT", "0"), row.get("MAT", "0"))
            fid = row["FID"] if "FID" in row else "x".join(pair)
            groups.setdefault(fid, []).append([index, pair])

        families = []
        for fid, members in groups.items():
            pairs = Counter( pair for _, pair in members
                             if not all(p in LeadInspector.MISSING_PARENT for p in pair) )
            pat, mat = pairs.most_common(1)[0][0] if pairs else (None, None)
            families.append({ "fid": fid, "members": [ index for index, _ in members ],
                              "pat": None if pat in LeadInspector.MISSING_PARENT else pat,
                              "mat": None if mat in LeadInspector.MISSING_PARENT else mat,
                              "offspring": [ index for index, pair in members if pairs and pair == (pat, mat) ] })
        return sorted(families, key=lambda f: (not f["fid"].isdigit(), int(f["fid"]) if f["fid"].isdigit() else 0,
                                               f["fid"]))

    @property
    def hasParents(self):
        '''
        Whether any family has a parent with genotypes.
        '''
        return any( parent in self.genotypes for f in self.families for parent in [f["pat"], f["mat"]] )

    def genotype(self, iid, leadIndex):
        '''
        Returns a sample's effect allele hard call (0/1/2) at a lead, or None if
        the sample has no genotypes or is missing a call.
        '''
        if iid is None or iid not in self.genotypes or np.isnan(self.genotypes[iid][leadIndex]):
            return None
        return int(self.genotypes[iid][leadIndex])

    @property
    def leadHeader(self):
        return LeadInspector.LEAD_COLUMNS + self.trait.LEAD_DIRECTION_COLUMNS + ["N_CARRIERS"] + \
               self.trait.LEAD_CARRIER_COLUMNS + LeadInspector.LEAD_GENOTYPE_COLUMNS

    @property
    def familyHeader(self):
        return ["LEAD_ID", "FID", "N"] + self.trait.FAMILY_COLUMNS + ["N_CARRIERS"] + \
               self.trait.FAMILY_CARRIER_COLUMNS + LeadInspector.FAMILY_PARENT_COLUMNS

    def lineage(self, iid, leadIndex):
        '''
        Traces an effect allele carrier's allele back through its genotyped
        parents that also carry it, like 'child<parent<grandparent'; where both
        parents carry it, the chain ends with 'pat|mat'.
        '''
        chain, seen = [iid], {iid}
        while iid in self.sampleDict:
            row = self.sampleDict[iid]
            carriers = [ p for p in [row.get("PAT"), row.get("MAT")]
                         if p not in seen and (self.genotype(p, leadIndex) or 0) > 0 ]
            if len(carriers) != 1:
                if len(carriers) == 2:
                    chain.append("|".join(carriers))
                break
            iid = carriers[0]
            chain.append(iid)
            seen.add(iid)
        return "<".join(chain)

    def family_origin(self, family, leadIndex, offspringCarriers):
        '''
        Determines which parent(s) of a family carry the effect allele.

        Returns:
            origin -- a list of the carrying parents' IIDs (suffixed with
                      '(inferred)' when deduced from offspring carriers and the
                      other parent), ['none'] if neither carries, or [] if it
                      cannot be assessed
        '''
        parents = [ p for p in [family["pat"], family["mat"]] if p is not None ]
        known = { p: self.genotype(p, leadIndex) for p in parents if self.genotype(p, leadIndex) is not None }
        if known == {}:
            return []
        carriers = [ p for p, g in known.items() if g > 0 ]
        if carriers != []:
            return carriers
        unknown = [ p for p in parents if p not in known ]
        if unknown != [] and offspringCarriers > 0:
            return [ f"{p}(inferred)" for p in unknown ]
        return ["none"]

    def mendel(self, family, leadIndex):
        '''
        Counts the family's offspring whose genotype is impossible given their
        parents' genotypes; a parent without genotypes may transmit either allele.

        Returns:
            errors, tested -- integers of offspring with errors, and offspring tested
        '''
        gametes = { 0: {0}, 1: {0, 1}, 2: {1}, None: {0, 1} }
        patG, matG = self.genotype(family["pat"], leadIndex), self.genotype(family["mat"], leadIndex)
        if patG is None and matG is None:
            return 0, 0
        allowed = { a + b for a in gametes[patG] for b in gametes[matG] }
        calls = [ int(c) for c in self.testedDosages[family["offspring"], leadIndex] if not np.isnan(c) ]
        return sum( c not in allowed for c in calls ), len(calls)

    def family_rows(self, leadIndex):
        '''
        Returns the per-family output rows of a lead, and a summary across families.
        '''
        lead = self.leads[leadIndex]
        calls = self.testedDosages[:, leadIndex]
        rows, carrierFamilies, origins, mendelErrors, mendelTested = [], [], [], 0, 0
        for family in self.families:
            memberCalls, memberPhenotypes = calls[family["members"]], self.phenotypes[family["members"]]
            isCarrier, hasCall = memberCalls > 0, ~np.isnan(memberCalls) # NaN compares False
            numCarriers = int(isCarrier.sum())
            if numCarriers > 0:
                carrierFamilies.append([numCarriers, family["fid"], len(family["members"])])

            parentValues, mendelValues = ["."] * 6, [".", "."]
            if self.hasParents:
                origin = self.family_origin(family, leadIndex, int((calls[family["offspring"]] > 0).sum()))
                lineages = [ self.lineage(p, leadIndex) for p in origin if p in self.genotypes ]
                errors, tested = self.mendel(family, leadIndex)
                mendelErrors, mendelTested = mendelErrors + errors, mendelTested + tested
                if numCarriers > 0:
                    origins += [ p for p in origin if p != "none" ]
                parentGenotypes = [ self.genotype(p, leadIndex) for p in [family["pat"], family["mat"]] ]
                parentValues = [family["pat"] or ".", family["mat"] or "."] + \
                               [ str(g) if g is not None else "." for g in parentGenotypes ] + \
                               [",".join(origin) or ".", ";".join(lineages) or "."]
                mendelValues = [str(errors), str(tested)]
            rows.append([lead["id"], family["fid"], str(len(family["members"]))] +
                        self.trait.family_values(memberPhenotypes) + [str(numCarriers)] +
                        self.trait.family_carrier_values(memberPhenotypes, isCarrier, hasCall) +
                        parentValues + mendelValues)

        carrierFamilies.sort(key=lambda x: -x[0])
        summary = { "carrierFamilies": ",".join( f"{fid}:{n}/{size}" for n, fid, size in carrierFamilies ) or ".",
                    "origin": ";".join(sorted(set(origins))) or ".",
                    "mendel": [str(mendelErrors), str(mendelTested)] if self.hasParents else [".", "."] }
        return rows, summary

    def link(self, leadIndex):
        '''
        Finds the stronger lead variant on another chromosome in the strongest LD
        with a lead.

        Returns:
            values -- a list of the LINK_STATUS, LINKED_TO and LINK_R to
                      EXPECTED_OVERLAP output values, with LINK_STATUS 'linked'
                      before links are followed to their root
            partner -- an integer index of the lead in the strongest LD, or None
        '''
        lead = self.leads[leadIndex]
        candidates = [ j for j in range(leadIndex) if self.leads[j]["chrom"] != lead["chrom"] ]
        x = self.testedImputed[:, leadIndex]
        if leadIndex == 0:
            return ["top"] + ["."] * 7, None
        if np.std(x) == 0:
            return ["monomorphic"] + ["."] * 7, None
        best, numCompared = None, 0
        for j in candidates:
            y = self.testedImputed[:, j]
            if np.std(y) == 0:
                continue
            numCompared += 1
            r = float(np.corrcoef(x, y)[0, 1])
            if best is None or r * r > best[1] ** 2:
                best = [j, r]
        if best is None:
            return ["independent"] + ["."] * 7, None
        j, r = best
        other = self.leads[j]

        # Overlap of this lead's carriers with carriers of the other lead's allele in coupling with it
        xCalls, yCalls = self.testedDosages[:, leadIndex], self.testedDosages[:, j]
        hasCalls = ~np.isnan(xCalls) & ~np.isnan(yCalls)
        xCarrier = xCalls[hasCalls] > 0
        yCoupled = yCalls[hasCalls] > 0 if r >= 0 else yCalls[hasCalls] < 2
        overlap = f"{(xCarrier & yCoupled).sum() / xCarrier.sum():.3f}" if xCarrier.sum() > 0 else "."

        empiricalP = self.empirical_p(r * r)
        adjustedP = min(1.0, empiricalP * numCompared)
        isLinked = r * r >= self.linkR2 if self.linkR2 is not None else adjustedP < self.alpha
        return ["linked" if isLinked else "independent", other["id"], f"{r:.3f}", f"{r * r:.3f}",
                f"{empiricalP:.3g}", f"{adjustedP:.3g}", overlap, f"{yCoupled.mean():.3f}"], j

    def resolve_links(self):
        '''
        Finds each lead's link, following links to their root.

        Returns:
            linkValues -- a list (per lead) of the LINK_STATUS to EXPECTED_OVERLAP
                          output values
        '''
        roots, linkValues = [], []
        for leadIndex, lead in enumerate(self.leads):
            values, partner = self.link(leadIndex)
            isLinked = values[0] == "linked"
            roots.append(roots[partner] if isLinked else leadIndex) # partners precede leads
            root = self.leads[roots[-1]]
            if isLinked and root["chrom"] == lead["chrom"]:
                values[0] = "own-chromosome"
            linkValues.append(values[:2] + [root["id"] if isLinked else "."] + values[2:])
        return linkValues

    def write(self, leadsFileName, familiesFileName):
        '''
        Writes the per-lead and per-lead-and-family inspections.

        Returns:
            counts -- a dictionary with structure like:
                      {
                          "leads": int, # lead variants written
                          "linked": int, # those linked to a stronger lead elsewhere
                          "withMendelErrors": int # those with Mendelian errors in at least
                                                  # MENDEL_RATE_FLAG of offspring tested
                      }
        '''
        if self.nullR2 is None:
            raise ValueError("set_null() must be called before write()")
        counts = { "leads": 0, "linked": 0, "withMendelErrors": 0 }
        with open(leadsFileName, "w") as leadsOut, open(familiesFileName, "w") as familiesOut:
            leadsOut.write("\t".join(self.leadHeader) + "\n")
            familiesOut.write("\t".join(self.familyHeader) + "\n")
            allLinkValues = self.resolve_links()
            for leadIndex, lead in enumerate(self.leads):
                stats = self.leadStats[lead["id"]]
                calls = self.testedDosages[:, leadIndex]
                isCarrier, hasCall = calls > 0, ~np.isnan(calls)
                familyRows, summary = self.family_rows(leadIndex)
                linkValues = allLinkValues[leadIndex]
                leadsOut.write("\t".join(
                    [lead["id"], lead["chrom"], str(lead["pos"]), lead["p"], "+".join(lead["sources"])] +
                    [ stats[c] for c in ["Allele1", "Allele2", "AF_Allele2", "BETA", "SE"] ] +
                    self.trait.lead_direction(stats, self.phenotypes, calls) + [str(int(isCarrier.sum()))] +
                    self.trait.carrier_values(self.phenotypes, isCarrier, hasCall) +
                    [str(stats["nWindow"]), str(stats["nWindowSuggestive"]),
                     summary["carrierFamilies"], summary["origin"]] + summary["mendel"] + linkValues) + "\n")
                for row in familyRows:
                    familiesOut.write("\t".join(row) + "\n")
                counts["leads"] += 1
                counts["linked"] += int(linkValues[0] == "linked")
                if summary["mendel"][1] not in {".", "0"}:
                    errors, tested = int(summary["mendel"][0]), int(summary["mendel"][1])
                    counts["withMendelErrors"] += int(errors / tested >= LeadInspector.MENDEL_RATE_FLAG)
        return counts

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
