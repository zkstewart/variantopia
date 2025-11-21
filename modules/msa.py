import os, sys
import pandas as pd
import numpy as np
import scipy.spatial.distance as ssd
from Bio import SeqIO
from scipy.cluster import hierarchy

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter

class DirectoryNotFoundError(Exception):
    pass

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
    
    def calculate_identity(self, seqid1, seqid2):
        '''
        Calculates the identity between a pair of sequences contained within this MSA.
        
        Parameters:
            seqid1 / seqid2 -- strings that identify sequences found within self.rowNames
                               and obtained by self[seqid1] and self[seqid2]
        '''
        if not seqid1 in self:
            raise KeyError(f"'{seqid1}' not found within this MSA")
        if not seqid2 in self:
            raise KeyError(f"'{seqid2}' not found within this MSA")
        
        if seqid1 == seqid2:
            return 1.0
        
        pwDf = self[[seqid1, seqid2]].loc[:, (self[[seqid1, seqid2]] != "-").any(axis=0)]
        if pwDf.shape[1] == 0:
            return 0.0
        else:
            return (pwDf.loc[seqid1] == pwDf.loc[seqid2]).sum() / pwDf.shape[1]
    
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
        for colIndex in self.df.columns:
            yield self.df[colIndex]
    
    def load_msa(self, msaFile):
        msaDict = {}
        with read_gz_file(msaFile) as fileIn:
            records = SeqIO.parse(fileIn, "fasta")
            for record in records:
                msaDict[record.id] = list(record.seq)
        self.df = pd.DataFrame.from_dict(msaDict, orient="index")
    
    def clean_gaps(self):
        '''
        Cleans up this MSA to remove any columns that are solely gaps. This might occur if
        you've been deleting sequences or otherwise modifying things.
        '''
        # Find and drop columns
        toDrop = []
        for colIndex in self.df.columns:
            column = self.df[colIndex]
            if (column == "-").all():
                toDrop.append(colIndex)
        self.df.drop(columns=toDrop, inplace=True)
    
    def keep(self, idsList):
        '''
        Keeps only IDs listed in idsList for this MSA object.
        
        Parameters:
            idsList -- a list containing strings which must occur within this MSA object.
        Modifies:
            self.df --  the underlying MSA representation is changed to retain only the sequences
                        identified by idsList. Runs .clean_gaps() automatically.
        '''
        for seqID in idsList:
            if not seqID in self:
                raise KeyError(f"'{seqID}' cannot be kept as it is not contained within this MSA")
        
        toKeep = set(idsList) # remove potential duplicates
        toKeep = [ x for x in self.rowNames if x in toKeep ] # maintain order of original MSA
        self.df = self.df.loc[toKeep]
        self.clean_gaps()
    
    def remove(self, idsList):
        '''
        Deletes several IDs as listed in idsList all at once from this MSA object. 
        
        Parameters:
            idsList -- a list containing strings which must occur within this MSA object.
        Modifies:
            self.df --  the underlying MSA representation is changed akin to performing
                        'del self[seqid]' for each value in idsList. Runs .clean_gaps()
                        automatically.
        '''
        for seqID in idsList:
            if not seqID in self:
                raise KeyError(f"'{seqID}' cannot be removed as it is not contained within this MSA")
        
        toRemove = set(idsList)
        if len(toRemove) == self.nrow:
            raise ValueError(".remove cannot be performed as it would remove all members of this MSA")
        
        toKeep = set(self.rowNames).difference(toRemove) # invert deletion list to become toKeep list
        toKeep = [ x for x in self.rowNames if x in toKeep ] # maintain order of original MSA
        self.df = self.df.loc[toKeep]
        self.clean_gaps()
    
    def hierarchical_reorder(self):
        '''
        Applies a simple hierarchical clustering algorithm to a pairwise matrix of sequence
        identities and reorders the MSA inplace corresponding to this.
        '''
        # Format a pairwise matrix of distances
        pwDict = { x: {y: None for y in self.rowNames } for x in self.rowNames }
        for x, yDict in pwDict.items():
            for y in yDict.keys():
                pwDict[x][y] = identity = self.calculate_identity(x, y)
        pwDf = pd.DataFrame.from_dict(pwDict)
        pwDf = pwDf.apply(lambda x: 1-x) # go from identity (/similarity) to distance)
        
        # Run hierarchical clustering of sequence distances
        distArray = ssd.squareform(pwDf)
        hier = hierarchy.linkage(distArray, method="ward")  # you can use other methods
        clusterLabels = hierarchy.fcluster(hier, 0, criterion="distance")
        
        # Get the ordered sequences from clustering labels
        sequenceSorting = sorted(zip(pwDf.columns.tolist(), clusterLabels), key = lambda x: x[1])
        sequencesSorted = [ seqID for seqID, clustNum in sequenceSorting ]
        
        # Reorder the underlying msa DataFrame
        self.df = self.df.loc[sequencesSorted]
    
    def write(self, outputFileName):
        '''
        Writes an output FASTA multiple sequence alignment. File overwriting is not allowed.
        
        Parameters:
            outputFileName -- a string pointing to a file location that does not already exist.
                              If full path is given, the parent folder(s) must already exist.
        '''
        outputFileName = os.path.abspath(outputFileName)
        if os.path.exists(outputFileName):
            raise FileExistsError(f"Cannot .write MSA to '{outputFileName}' as it already exists!")
        parentDir = os.path.dirname(outputFileName)
        if not os.path.isdir(parentDir):
            raise DirectoryNotFoundError(f"Cannot .write MSA '{os.path.basename(outputFileName)}' as " + 
                                         f"its parent directory '{parentDir}' does not exist!")
        
        with GzCapableWriter(outputFileName) as fileOut:
            for seqID, seqRow in self:
                seq = "".join(seqRow)
                fileOut.write(f">{seqID}\n{seq}\n")
    
    def __delitem__(self, key):
        if not key in self:
            raise KeyError(f"'{key}' cannot be deleted as it is not contained within this MSA")
        
        # Drop the sequence row
        self.df.drop(key, inplace=True)
        self.clean_gaps() # just in case this introduced gaps
    
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
    with GzCapableWriter(args.outputFileName) as fileOut:
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
    with GzCapableWriter(args.outputFileName) as fileOut:
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
