import os, sys
import pandas as pd
import numpy as np
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, WriteGzFile

class MSATopia:
    def __init__(self, msaFile):
        self.msaFile = msaFile
        self.df = None
        self.load_msa(self.msaFile)
        self.isMSATopia = True # flag for easier type checking
    
    @staticmethod
    def count_residues(column):
        '''
        Takes a pandas Series column and returns a sorted list of tuples
        containing each residue (as a string) and its count.
        
        Parameters:
            column -- a pd.Series object representing a column of an MSA DataFrame.
        Returns:
            residueCounts -- a dictionary where each residue (in uppercase) is a key and each value is
                             its count. Dictionary is sorted from most to least common residue.
        '''
        unique, counts = np.unique(column.values, return_counts=True)
        residueCounts = dict(sorted(zip([ x.upper() for x in unique ], counts), key=lambda x: x[1], reverse=True)) # x == (residue, count)
        
        return residueCounts
    
    @staticmethod
    def consensus_and_alt_residues(column, asSet=False):
        '''
        Takes a pandas Series column and returns a string indicating the consensus residue
        and a list of strings indicating the alternative residue(s).
        
        Parameters:
            column -- a pd.Series object representing a column of an MSA DataFrame.
            asSet -- a boolean indicating whether to return the alternatives as a set (True)
                     or list (False); default is False
        Returns:
            consensus -- a string indicating the consensus residue (in uppercase)
            alternatives -- a list of strings indicating the alternative residues (in uppercase)
        '''
        residueCounts = MSATopia.count_residues(column)
        
        consensus = None
        alternatives = []
        for residue, count in residueCounts.items():
            if consensus is None:
                consensus = residue
            else:
                alternatives.append(residue)
        
        if asSet:
            alternatives = set(alternatives)
        return consensus, alternatives
    
    @staticmethod
    def calculate_gc(column):
        '''
        Calculate the GC content of a pandas Series column.
        
        Parameters:
            column -- a pd.Series object representing a column of an MSA DataFrame.
        Returns:
            gc -- the GC content as a float.
        '''
        residueCounts = MSATopia.count_residues(column)
        numGC = residueCounts.get("G", 0) + residueCounts.get("C", 0)
        return numGC / len(column)
    
    @staticmethod
    def calculate_mac(column):
        '''
        Calculate the minor allele count (MAC; the count of the second most common allele/residue)
        of a pandas Series column.
        
        Parameters:
            column -- a pd.Series object representing a column of an MSA DataFrame.
        Returns:
            mac -- the minor allele count (MAC) for the column as a float.
        '''
        residueCounts = MSATopia.count_residues(column)
        if len(residueCounts) < 2:
            return 0.0
        
        isSecond = False
        for residue, count in residueCounts.items():
            if isSecond:
                return count
            isSecond = True
    
    @staticmethod
    def calculate_maf(column):
        '''
        Calculate the minor allele frequency (MAF; the frequency of the second most common allele/residue)
        of a pandas Series column.
        
        Parameters:
            column -- a pd.Series object representing a column of an MSA DataFrame.
        Returns:
            maf -- the minor allele frequency (MAF) for the column as a float.
        '''
        residueCounts = MSATopia.count_residues(column)
        if len(residueCounts) < 2:
            return 0.0
        
        isSecond = False
        for residue, count in residueCounts.items():
            if isSecond:
                return count / len(column)
            isSecond = True
    
    @staticmethod
    def calculate_gaprate(column):
        '''
        Calculate the gaprate (i.e., the proportion that is '-')
        of a pandas Series column.
        
        Parameters:
            column -- a pd.Series object representing a column of an MSA DataFrame
        Returns:
            gaprate -- the GC content as a float
        '''
        residueCounts = MSATopia.count_residues(column)
        numGaps = residueCounts.get("-", 0)
        return numGaps / len(column)
    
    @staticmethod
    def calculate_uniqueness(column, groups):
        '''
        Calculate the uniqueness (i.e., the difference in how many times a variant
        occurs in two different groups) of a pandas Series column.
        
        Parameters:
            column -- a pd.Series object representing a column of an MSA DataFrame
            groups -- an iterable with equal length and ordering to the column containing
                      values of 1 or 2, indicating whether the sequence column belongs to
                      group 1 or 2. Use 0 or None to ignore a sample from the calculation.
        Returns:
            uniqueness -- the absolute difference in variant occurrence as a ratio from
                          0 (no difference in variant occurrence between groups) to 1
                          (variant exclusively occurs in one group and not the other)
        '''
        if groups == None:
            raise ValueError("calculate_uniqueness requires the 'groups' metadata variable")
        
        consensus, alts = MSATopia.consensus_and_alt_residues(column)
        
        group1Variants, group2Variants = 0, 0
        group1Samples, group2Samples = 0, 0
        for seqID, residue in column.to_dict().items():
            # Find which group this sample belongs to
            prefixMatch = [ (prefix, group) for prefix, group in groups.items() if seqID.startswith(prefix) ]
            if len(prefixMatch) == 0:
                raise ValueError(f"'{seqID}' does not start with a prefix in the groups metadata dictionary")
            if len(prefixMatch) > 1:
                foundPrefixes = [ x for x,y in prefixMatch ]
                raise ValueError(f"'{seqID}' matches multiple group metadata values ({foundPrefixes})")
            
            prefix, group = prefixMatch[0]
            
            # Tally group and variant occurrence
            if str(group) == "1":
                group1Samples += 1
                if residue != consensus:
                    group1Variants += 1
            elif str(group) == "2":
                group2Samples += 1
                if residue != consensus:
                    group2Variants += 1
        
        # Convert into "uniqueness percent"
        try:
            group1Ratio = group1Variants / group1Samples
        except ZeroDivisionError:
            group1Ratio = 0.0
        try:
            group2Ratio = group2Variants / group2Samples
        except ZeroDivisionError:
            group2Ratio = 0.0
        
        return abs(group1Ratio - group2Ratio)
    
    @property
    def msaFile(self):
        return self._msaFile
    
    @msaFile.setter
    def msaFile(self, value):
        if not isinstance(value, str):
            raise ValueError("MSA file must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"MSA file '{value}' is not a file")
        
        self._msaFile = value
    
    @property
    def nrow(self):
        return len(self.df)
    
    @property
    def ncol(self):
        return len(self.df.columns)
    
    @property
    def rowNames(self):
        return self.df.index.tolist()
    
    @property
    def rows(self):
        for index, rowSeries in self.df.iterrows():
            yield rowSeries
    
    @property
    def colNames(self):
        return self.df.columns.tolist()
    
    @property
    def columns(self):
        for col in self.df.columns:
            yield self.df[col]
    
    def load_msa(self, msaFile):
        msaDict = {}
        with read_gz_file(msaFile) as fileIn:
            records = SeqIO.parse(fileIn, "fasta")
            for record in records:
                msaDict[record.id] = list(record.seq)
        self.df = pd.DataFrame.from_dict(msaDict, orient="index")
    
    def __getitem__(self, key):
        return self.df.loc[key]
    
    def __len__(self):
        return len(self.df)
    
    def __iter__(self):
        yield from self.df.iterrows()
    
    def __contains__(self, value):
        return value in self.df.index
    
    def __repr__(self):
        return "<MSATopia object;file='{0}';num_genes={1};num_residues={2}>".format(
            self.msaFile,
            self.nrow,
            self.ncol
        )

def msa_to_variant_report(args):
    with WriteGzFile(args.outputFileName) as fileOut:
        fileOut.write("#gene\tposition_number\tconsensus_residue\tvariant_residue\tseqs_with_variant\n")
        
        for msaFile in args.msaFiles:
            msa = MSATopia(msaFile)
            msaID = os.path.basename(msaFile).rsplit(".", maxsplit=1)[0]
            
            foundStop = set()
            for position, column in enumerate(msa.columns):
                # Get the consensus and alt residue(s)
                consensus, alt = MSATopia.consensus_and_alt_residues(column, asSet=True)
                if len(alt) == 0:
                    continue # no variants at this position
                
                # Write out variant information
                sampleResidueDict = column.to_dict()
                for altResidue in alt:
                    hasVariant = sorted([
                        sampleID
                        for sampleID, residue in sampleResidueDict.items()
                        if (not sampleID in foundStop) and (residue.upper() == altResidue)
                    ])
                    fileOut.write(f"{msaID}\t{position+1}\t{consensus}\t" +
                                  f"{altResidue}\t{', '.join(hasVariant)}\n")
                    
                    # Conditionally handle reportUntilStop
                    if args.reportUntilStop and altResidue == "*":
                        foundStop.update(hasVariant)

def msa_to_sequence_report(args):
    with WriteGzFile(args.outputFileName) as fileOut:
        fileOut.write("#gene\tsequence_id\tvariants\n")
        
        for msaFile in args.msaFiles:
            msa = MSATopia(msaFile)
            msaID = os.path.basename(msaFile).rsplit(".", maxsplit=1)[0]
            orderedSampleIDs = sorted(msa.rowNames)
            
            # Find all variants in the MSA
            sampleVariantsDict = { sampleID : [] for sampleID in msa.rowNames }
            foundStop = set()
            for position, column in enumerate(msa.columns):
                # Get the consensus and alt residue(s)
                consensus, alt = MSATopia.consensus_and_alt_residues(column, asSet=True)
                if len(alt) == 0:
                    continue # no variants at this position
                
                # Write out variant information
                for sampleID, residue in column.to_dict().items():
                    if (not sampleID in foundStop) and (residue.upper() in alt):
                        sampleVariantsDict[sampleID].append(f"{position+1}:{consensus}:{residue}")
                    if args.reportUntilStop and residue == "*":
                        foundStop.add(sampleID)
            
            # Write out sequence information
            for sampleID in orderedSampleIDs:
                variantList = sampleVariantsDict[sampleID]
                if len(variantList) == 0:
                    continue # omit sequences with no variants
                fileOut.write(f"{msaID}\t{sampleID}\t{', '.join(variantList)}\n")
