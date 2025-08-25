import os
import pandas as pd
import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt
from Bio import SeqIO

from .parsing import read_gz_file

#msaFile = "/mnt/f/plant_group/anuradha/upr_annotation/binge/msas/aligned_aa/cluster_0.aa"

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

class MSAPlot:
    STANDARD_DIMENSION = 5
    
    def __init__(self, msaFiles, statistic,
                 width=None, height=None):
        '''
        Parameters:
            msaFiles -- a list of strings indicating the MSA files to plot
            statistic -- a string indicating the statistic to plot; currently supported
                         values are "gc", "mac", "maf", "gaprate"
            width -- an integer indicating the width of the plot in inches; if None,
                     defaults to MSAPlot.STANDARD_DIMENSION * number of regions
            height -- an integer indicating the height of the plot in inches; if None,
                      defaults to MSAPlot.STANDARD_DIMENSION * number of rows
        '''
        # Behaviour parameters
        self.statistic = statistic
        
        # Data files
        self.msaFiles = msaFiles
        
        # Aesthetic parameters
        self.width = width
        self.height = height
        
        # Default values unset by initialization
        self._colourMap = None
    
    @property
    def msaFiles(self):
        return self._msaFiles
    
    @msaFiles.setter
    def msaFiles(self, value):
        if not isinstance(value, list):
            raise TypeError("msaFiles must be a list of strings")
        for item in value:
            if not isinstance(item, str):
                raise TypeError("msaFiles must be a list of strings")
            if not os.path.isfile(item):
                raise FileNotFoundError(f"msaFile '{item}' is not a file")
        
        self._msaFiles = value
    
    @property
    def statistic(self):
        return self._statistic
    
    @statistic.setter
    def statistic(self, value):
        ACCEPTED_STATISTICS = ["gc", "mac", "maf", "gaprate"]
        if value in ACCEPTED_STATISTICS:
            self._statistic = value
        else:
            raise TypeError(f"statistic must be a string from {ACCEPTED_STATISTICS}, not '{value}'")
    
    @property
    def statFunction(self):
        if self.statistic == "gc":
            return MSATopia.calculate_gc
        elif self.statistic == "mac":
            return MSATopia.calculate_mac
        elif self.statistic == "maf":
            return MSATopia.calculate_maf
        elif self.statistic == "gaprate":
            return MSATopia.calculate_gaprate
        else:
            raise ValueError(f"self.statistic=='{self.statistic}' is not handled by .statFunction")
    
    @property
    def colourMap(self):
        if self._colourMap is None:
            return "viridis"
        return self._colourMap
    
    @colourMap.setter
    def colourMap(self, value):
        ACCEPTED_COLOUR_MAPS = ["viridis", "Greys", "GnBu", "RdBu"]
        if value in ACCEPTED_COLOUR_MAPS:
            self._colourMap = value
        else:
            raise TypeError(f"colourMap must be a string from {ACCEPTED_COLOUR_MAPS}, not '{value}'")
    
    @property
    def cmap(self):
        return getattr(plt.cm, self.colourMap)
    
    @property
    def dtype(self):
        if self.statistic == "gc":
            return np.float64
        elif self.statistic == "mac":
            return np.float64
        elif self.statistic == "maf":
            return np.float64
        elif self.statistic == "gaprate":
            return np.float64
        else:
            raise ValueError(f"self.statistic=='{self.statistic}' is not handled by .dtype")
    
    @property
    def width(self):
        if self._width is None:
            return MSAPlot.STANDARD_DIMENSION
        return self._width
    
    @width.setter
    def width(self, value):
        if value == None:
            self._width = None
            return
        else:
            if not isinstance(value, int):
                raise TypeError("width must be an integer")
            if value < 1:
                raise ValueError(f"width must be >= 1")
            
            self._width = value
    
    @property
    def height(self):
        if self._height is None:
            return MSAPlot.STANDARD_DIMENSION
        return self._height
    
    @height.setter
    def height(self, value):
        if value == None:
            self._height = None
            return
        else:
            if not isinstance(value, int):
                raise TypeError("height must be an integer")
            if value < 1:
                raise ValueError(f"height must be >= 1")

            self._height = value
    
    def get_x_y(self, msaFile):
        '''
        Returns the x list and y array needed for broken_barh plotting of MSA
        statistics. Method of operation is dictated by self.statistic.
        Returned objects have a length equal to the number of residues in the alignment,
        where the y array values indicate the relevant statistic for each position.
        
        Parameters:
            msaFile -- a string indicating the MSA file to plot
        Returns:
            x -- a list containing tuples of (x, 1) where x is the position index
                 and 1 is the length of the position as needed by broken_barh
            y -- a np.array containing the statFunction transformed values for
                 each variant/window
        '''
        # Load the MSA
        msa = MSATopia(msaFile)
        
        # Lay out the y values in a numpy array
        y = []
        for i, column in enumerate(msa.columns):
            y.append(self.statFunction(column))
        y = np.array(y, dtype=self.dtype)
        
        # Create the x array
        x = [ (i, 1) for i in np.arange(0, len(y)) ]
        
        return x, y
    
    def plot(self, outputFileName):
        '''
        Plots the specified .statistic for the MSA.
        
        Parameters:
            outputFileName -- a string indicating the location to write plot output to;
                              assumes file name is pre-validated to ensure no overwriting
                              and to have a valid filename suffix to control the file format
        '''
        SPACING = 0.1
        
        # Obtain data for plotting
        msaLabels = []
        msaArrays = []
        for msaFile in self.msaFiles:
            msaLabels.append(os.path.basename(msaFile).rsplit(".", maxsplit=1)[0])
            msaArrays.append(self.get_x_y(msaFile))
        
        # Sort data from longest to shortest
        msaArrays.sort(key = lambda x: len(x[0])) # [x1,y1] each have same length
        
        # Obtain the minimum and maximum values being plotted
        minValue, maxValue, maxLen = np.inf, -np.inf, -np.inf
        for x, y in msaArrays:
            statMin, statMax = np.min(y), np.max(y)
            if statMin < minValue:
                minValue = statMin
            if statMax > maxValue:
                maxValue = statMax
            if len(y) > maxLen:
                maxLen = len(y)
        
        # Establish colour map
        norm = matplotlib.colors.Normalize(vmin=minValue, vmax=maxValue)
        
        # Configure plot
        fig = plt.figure(figsize=(self.width, self.height), tight_layout=True)
        ax = plt.axes()
        
        # Configure x-axis labels
        ax.set_xlabel(f"Length (bp)", fontweight="bold")
        ax.set_xlim(0, maxLen)
        
        # Plot each MSA
        ongoingCount = 0.5 # this centers the contig label
        for x, y in msaArrays:
            ax.broken_barh(x, (ongoingCount+SPACING, 1-(SPACING*2)), facecolors=self.cmap(norm(y)))
            ongoingCount += 1
        ax.set_ylim(0.5+SPACING, ongoingCount-SPACING)
        
        # Indicate y-axis labels
        ax.set_ylabel("") # no y-axis title label
        ax.set_yticks(range(1, len(msaLabels)+1))
        ax.set_yticklabels(msaLabels)
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical",
                     label="GC" if self.statistic == "gc" else \
                           "Minor Allele Count" if self.statistic == "mac" else \
                           "Minor Allele Frequency" if self.statistic == "maf" else \
                           "Gap Rate" if self.statistic == "gaprate" else \
                           "labeltypenothandled")
        
        # Save output file
        plt.savefig(outputFileName)
        plt.close()

def msa_to_plot(args):
    plot = MSAPlot(args.msaFiles, statistic=args.statistic,
                   width=args.width, height=args.height
    )
    plot.colourMap = args.colourMap
    plot.plot(args.outputFileName)
