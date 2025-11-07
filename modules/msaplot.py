import os, sys
import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt
from pymsaviz import MsaViz

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from msa import MSATopia
from parsing import parse_metadata_groups

class MSAPlot:
    STANDARD_DIMENSION = 5
    
    def __init__(self, msaFiles, width=None, height=None):
        '''
        Parameters:
            msaFiles -- a list of strings indicating the MSA files to plot
            width -- an integer indicating the width of the plot in inches; if None,
                     defaults to MSAPlot.STANDARD_DIMENSION * number of regions
            height -- an integer indicating the height of the plot in inches; if None,
                      defaults to MSAPlot.STANDARD_DIMENSION * number of rows
        '''
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
    
    def plot(self):
        raise NotImplementedError("Inheriting class must override")

class MSAPlotStats(MSAPlot):
    def __init__(self, msaFiles, statistic, width=None, height=None):
        '''
        Parameters:
            msaFiles -- a list of strings indicating the MSA files to plot
            statistic -- a string indicating the statistic to plot; currently supported
                         values are "gc", "mac", "maf", "gaprate", "uniqueness"
            width -- an integer indicating the width of the plot in inches; if None,
                     defaults to MSAPlot.STANDARD_DIMENSION * number of regions
            height -- an integer indicating the height of the plot in inches; if None,
                      defaults to MSAPlot.STANDARD_DIMENSION * number of rows
        '''
        super().__init__(msaFiles, width, height)
        
        # Behaviour parameters
        self.statistic = statistic
        
        # Default values unset by initialization
        self._colourMap = None
    
    @property
    def statistic(self):
        return self._statistic
    
    @statistic.setter
    def statistic(self, value):
        ACCEPTED_STATISTICS = ["gc", "mac", "maf", "gaprate", "uniqueness"]
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
        elif self.statistic == "uniqueness":
            return MSATopia.calculate_uniqueness
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
        elif self.statistic == "uniqueness":
            return np.float64
        else:
            raise ValueError(f"self.statistic=='{self.statistic}' is not handled by .dtype")
    
    def get_x_y(self, msaFile, groupDict=None):
        '''
        Returns the x list and y array needed for broken_barh plotting of MSA
        statistics. Method of operation is dictated by self.statistic.
        Returned objects have a length equal to the number of residues in the alignment,
        where the y array values indicate the relevant statistic for each position.
        
        Parameters:
            msaFile -- a string indicating the MSA file to plot
            groupDict -- (OPTIONAL) if self.statistic == "uniqueness" you must specify
                         a dictionary where sequence prefixes are keys and values are
                         1 (group1), 2 (group2), 0 or None (ignore). All sequences must
                         match a single unique prefix.
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
        for column in msa.columns:
            if self.statistic == "uniqueness":
                y.append(self.statFunction(column, groupDict))
            else:
                y.append(self.statFunction(column))
        y = np.array(y, dtype=self.dtype)
        
        # Create the x array
        x = [ (i, 1) for i in np.arange(0, len(y)) ]
        
        return x, y
    
    def plot(self, outputFileName, groupDict=None):
        '''
        Plots the specified .statistic for the MSA.
        
        Parameters:
            outputFileName -- a string indicating the location to write plot output to;
                              assumes file name is pre-validated to ensure no overwriting
                              and to have a valid filename suffix to control the file format
            groupDict -- (OPTIONAL) if self.statistic == "uniqueness" you must specify
                         a dictionary where sequence prefixes are keys and values are
                         1 (group1), 2 (group2), 0 or None (ignore). All sequences must
                         match a single unique prefix.
        '''
        SPACING = 0.1
        
        # Obtain data for plotting
        msaLabels = []
        msaArrays = []
        for msaFile in self.msaFiles:
            msaLabels.append(os.path.basename(msaFile).rsplit(".", maxsplit=1)[0])
            msaArrays.append(self.get_x_y(msaFile, groupDict))
        
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
        ax.set_xlabel(f"Length (residues)", fontweight="bold")
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
                           "Uniqueness" if self.statistic == "uniqueness" else \
                           "labeltypenothandled")
        
        # Save output file
        plt.savefig(outputFileName)
        plt.close()

class DirectoryNotFoundError(Exception):
    pass
class FileFormatError(Exception):
    pass

class MSAPlotAlignment(MSAPlot):
    STANDARD_WRAP = 60
    COLOURS = [ v for v in matplotlib.colors.TABLEAU_COLORS.values() ]
    
    def __init__(self, msaFiles, statistic, wrapLength=None, width=None, height=None):
        '''
        Parameters:
            msaFiles -- a list of strings indicating the MSA files to plot
            statistic -- a string indicating the statistic to plot; currently supported
                         values are "gc", "mac", "maf", "gaprate", "uniqueness"
            wrapLength -- an integer indicating the number of residues to interleave; if None,
                     defaults to MSAPlotAlignment.STANDARD_WRAP
            width -- an integer indicating the width of the plot in inches; if None,
                     defaults to MSAPlot.STANDARD_DIMENSION * number of regions
            height -- an integer indicating the height of the plot in inches; if None,
                      defaults to MSAPlot.STANDARD_DIMENSION * number of rows
        '''
        super().__init__(msaFiles, width, height)
        
        # Behaviour parameters
        self.statistic = statistic
        self.wrapLength = wrapLength
        
        # Default values unset by initialization
        self.domains = None
    
    @property
    def statistic(self):
        return self._statistic
    
    @statistic.setter
    def statistic(self, value):
        ACCEPTED_STATISTICS = ["maf", "gaprate", "uniqueness"]
        if value in ACCEPTED_STATISTICS:
            self._statistic = value
        else:
            raise TypeError(f"statistic must be a string from {ACCEPTED_STATISTICS}, not '{value}'")
    
    @property
    def statFunction(self):
        if self.statistic == "maf":
            return MSATopia.calculate_maf
        elif self.statistic == "gaprate":
            return MSATopia.calculate_gaprate
        elif self.statistic == "uniqueness":
            return MSATopia.calculate_uniqueness
        else:
            raise ValueError(f"self.statistic=='{self.statistic}' is not handled by .statFunction")
    
    @property
    def dtype(self):
        if self.statistic == "maf":
            return np.float64
        elif self.statistic == "gaprate":
            return np.float64
        elif self.statistic == "uniqueness":
            return np.float64
        else:
            raise ValueError(f"self.statistic=='{self.statistic}' is not handled by .dtype")
    
    @property
    def wrapLength(self):
        if self._wrapLength is None:
            return MSAPlotAlignment.STANDARD_WRAP
        return self._wrapLength
    
    @wrapLength.setter
    def wrapLength(self, value):
        if value == None:
            self._wrapLength = None
            return
        else:
            if not isinstance(value, int):
                raise TypeError("wrapLength must be an integer")
            if value < 1:
                raise ValueError(f"wrapLength must be >= 1")

            self._wrapLength = value
    
    def parse_domtblout(self, fileName, evalue=None):
        '''
        Input a HMMER3 domtblout file to be parsed and have domain predictions annotated
        atop the alignment plot. This function will automatically parse the file
        and store data in self.domains to be used by .plot().
        
        The Domains class used to parse this file needs to be imported from the annotarium
        repository or otherwise made available globally.
        
        Parameters:
            fileName -- a string indicating the HMMER .domtblout file with domain
                        predictions for all files listed under self.msaFiles
            evalue -- (OPTIONAL) a float indicating the E-value cutoff to enforce when parsing
                      domtblout results OR None to perform no filtration
        Sets:
            self.domains -- a dictionary pairing sequences to DomainFeature objects
                            which indicate the location to annotate protein domains.
        '''
        # Load in domain predictions
        if evalue == None:
            self.domains = Domains() # must be available globally,
        else:
            self.domains = Domains(evalue=evalue)
        
        self.domains.parse_domtblout(fileName)
        self.domains.resolve_overlaps() # preliminary overlap resolution per-sequence
    
    def _override_consensus_identity(self, msaFile, groupDict=None):
        '''
        Handles the overriding PyMSAViz's default consensus identity bar plot with an
        alternative statistic (dictated by self.statistic).
        
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
        for column in msa.columns:
            if self.statistic == "uniqueness":
                y.append(self.statFunction(column, groupDict))
            else:
                y.append(self.statFunction(column))
        y = np.array(y, dtype=self.dtype)
        
        # Re-scale y value to be 0->100
        y = y*100 # PyMSAViz does not use a 0->1 ratio
        
        # Override the method
        MsaViz._get_consensus_identity_list = lambda x, start, end: y[start:end]
    
    def plot(self, outputFileLocation, fileFormat, groupDict=None):
        '''
        Plots the specified .statistic for the MSA.
        
        If plotting domain locations (self.domains != None), the OverlapResolver class
        needs to be imported from the annotarium repository or otherwise made available globally.
        
        Parameters:
            outputFileLocation -- a string indicating the location to write plot output(s) to;
                                  file names will be set according to the original MSA file name
                                  sans its file suffix.
            fileFormat -- 
            groupDict -- (OPTIONAL) if self.statistic == "uniqueness" you must specify
                         a dictionary where sequence prefixes are keys and values are
                         1 (group1), 2 (group2), 0 or None (ignore). All sequences must
                         match a single unique prefix.
        '''
        ACCEPTED_OUTPUT_FORMATS = ["png", "pdf", "svg"]
        
        if not os.path.isdir(outputFileLocation):
            raise DirectoryNotFoundError(f"'{outputFileLocation}' is not a directory; cannot write outputs here")
        if not fileFormat in ACCEPTED_OUTPUT_FORMATS:
            raise FileFormatError(f"'{fileFormat}' is not recognised as a valid output format; choose amongst {ACCEPTED_OUTPUT_FORMATS}")
        
        for msaFile in self.msaFiles:
            # Derive the file prefix
            msaPrefix = os.path.basename(msaFile)
            if msaPrefix.endswith(".gz"):
                msaPrefix = msaPrefix.rsplit(".", maxsplit=1)[0]
            msaPrefix = msaPrefix.rsplit(".", maxsplit=1)[0]
            
            # Format the output file name; skip if it exists
            outputFileName = os.path.join(outputFileLocation, f"{msaPrefix}.{self.statistic}.{fileFormat}")
            if os.path.exists(outputFileName):
                print(f"# '{outputFileName}' already exists; skipping...")
                continue
            
            # Generate the MsaViz object
            mv = MsaViz(msaFile, wrap_length=self.wrapLength, show_consensus=True)
            
            # Identify domains if applicable
            msaDomains = []
            if self.domains != None:
                # Find all domains associated with sequences in this MSA
                for seqID, seq in zip(mv.id_list, mv.seq_list):
                    if seqID in self.domains:
                        seqDomains = self.domains[seqID]
                        
                        # Relate gene domain coordinates to the MSA
                        seqIndex = 1 # 1-based coordinate
                        indexMap = {} # key and values are stored as 1-based
                        for position in range(mv.alignment_length):
                            residue = seq[position] # must obtain residue as 0-based
                            
                            if residue != "-":
                                indexMap[seqIndex] = position+1 # then store position as 1-based
                                seqIndex += 1
                        
                        # Update domain predictions with translated MSA positions
                        for feature in seqDomains:
                            feature.start = indexMap[feature.start]
                            feature.end = indexMap[feature.end]
                        
                        msaDomains += seqDomains
                
                # Resolve overlapping predictions
                resolver = OverlapResolver() # must be available globally,
                msaDomains = resolver.resolve(msaDomains)
            
            # Modify the statistic to plot underneath the alignment
            if self.statistic != "maf": # consensus visual shows maf by default
                self._override_consensus_identity(msaFile, groupDict=groupDict)
            
            # Annotate domain locations
            for i, feature in enumerate(msaDomains):
                mv.add_text_annotation((feature.start, feature.end), feature.id,
                                       text_color=MSAPlotAlignment.COLOURS[i],
                                       range_color=MSAPlotAlignment.COLOURS[i])
            
            # Write output
            mv.savefig(outputFileName)

def msa_plot_stats(args):
    plot = MSAPlotStats(args.msaFiles, statistic=args.statistic,
                        width=args.width, height=args.height
    )
    plot.colourMap = args.colourMap
    
    # Parse metadata (if applicable)
    if args.statistic == "uniqueness":
        # Parse the metadata file
        groupDict = parse_metadata_groups(args.metadataGroups)
    else:
        groupDict = None
    
    plot.plot(args.outputFileName, groupDict=groupDict)

def msa_plot_alignment(args):
    plot = MSAPlotAlignment(args.msaFiles, statistic=args.statistic,
                            width=args.width, height=args.height,
                            wrapLength=args.wrapLength
    )
    
    # Parse metadata (if applicable)
    if args.statistic == "uniqueness":
        # Parse the metadata file
        groupDict = parse_metadata_groups(args.metadataGroups)
    else:
        groupDict = None
    
    # Parse the domtblout file (if applicable)
    if args.domtbloutFileName != None:
        try:
            sys.path.append(os.path.dirname(args.annotariumDir))
            global Domains
            global OverlapResolver
            from annotarium import Domains, OverlapResolver
        except ModuleNotFoundError:
            raise ModuleNotFoundError(f"Could not import Domains and OverlapResolver from '{args.annotariumDir}'")
        
        plot.parse_domtblout(args.domtbloutFileName)
    
    plot.plot(args.outputDirectory, args.fileFormat, groupDict=groupDict)
