import os, math
import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt

from Bio import SeqIO

from .vcf import VCFTopia
from .gff3 import GFF3Topia

class Plot:
    STANDARD_DIMENSION = 5
    
    INTERVALS = {
        1: [1, "bp"], # 1 bp
        10: [1, "bp"], # 10 bp
        100: [1, "bp"], # 100 bp
        1000: [1000, "Kb"], # 1 Kb
        10000: [1000, "Kb"], # 10 Kb
        100000: [1000, "Kb"], # 100 Kb
        1000000: [1000000, "Mb"], # 1 Mb
        10000000: [1000000, "Mb"], # 10 Mb
        100000000: [1000000, "Mb"], # 100 Mb
        1000000000: [1000000000, "Gb"], # 1 Gb
        10000000000: [1000000000, "Gb"], # 10 Gb
        100000000000: [1000000000, "Gb"], # 100 Gb
        1000000000000: [1000000000000, "Tb"], # 1 Tb
        10000000000000: [1000000000000, "Tb"], # 10 Tb
        100000000000000: [1000000000000, "Tb"], # 100 Tb
        1000000000000000: [1000000000000000, "Pb"] # that's got to future-proof it for a while
    }
    
    @staticmethod
    def bin_values(values, windowSize, summariseMethod=np.sum):
        '''
        Parameters:
            values -- a list of lists containing three values: [position, contigID, statistical value]
            windowSize -- an integer value indicating the size of the windows
            summariseMethod -- a function to summarise the values within each bin; defaults to np.sum
                               but can be set to other functions like np.mean, np.median, etc.
        Returns:
            windows -- a numpy array with length equal to the number of windows
                       that the values array divides into, with each window's value
                       being the result of applying the summariseMethod to the values
                       within that window.
        '''
        # Bin values into their respective windows
        windows = [
            [] for _ in range(math.ceil(len(values) / windowSize))
        ]
        for i, stat in enumerate(values):
            windowIndex = (i) // windowSize
            windows[windowIndex].append(stat)
        
        # Apply the summariseMethod to each window
        windows = np.array([
            summariseMethod(np.array(window)) if window else 0
            for window in windows
        ])
        return windows
    
    def __init__(self, statistic, feature, windowSize, vcf, gff3=None, genomeFile=None,
                 width=None, height=None):
        '''
        Parameters:
            statistic -- a string indicating the statistic to plot; currently supported
                         values are "snpnumber", "mac", "maf", "callrate", and "het"
            feature -- a string indicating the feature to plot; currently supported
                       values are "genes" and "chromosomes"
            windowSize -- a positive integer indicating the size of the window to use
                          for plotting
            vcf -- a VCFTopia object containing the VCF data to plot
            gff3 -- a GFF3Topia object or None if no gene models are to be plotted
            genomeFile -- a string indicating the path to the genome file in FASTA format
                          or None if no genome information is needed
            width -- an integer indicating the width of the plot in inches; if None,
                     defaults to Plot.STANDARD_DIMENSION * number of regions
            height -- an integer indicating the height of the plot in inches; if None,
                      defaults to Plot.STANDARD_DIMENSION * number of rows
        '''
        # Behaviour parameters
        self.statistic = statistic
        self.feature = feature
        self.windowSize = windowSize
        
        # Data files
        self.vcf = vcf
        self.gff3 = gff3
        self.genomeFile = genomeFile
        
        # Aesthetic parameters
        self.width = width
        self.height = height
        
        # Figure-related parameters (not to be set by user)
        self.fig = None
        self.axs = None
        self.rowNum = None
        
        # Default values unset by initialization
        self._colourMap = None
    
    @property
    def statistic(self):
        return self._statistic
    
    @statistic.setter
    def statistic(self, value):
        ACCEPTED_STATISTICS = ["snpnumber", "mac", "maf", "callrate", "het"]
        if value in ACCEPTED_STATISTICS:
            self._statistic = value
        else:
            raise TypeError(f"statistic must be a string from {ACCEPTED_STATISTICS}, not '{value}'")
    
    @property
    def feature(self):
        return self._feature
    
    @feature.setter
    def feature(self, value):
        ACCEPTED_FEATURES = ["genes", "chromosomes"]
        if value in ACCEPTED_FEATURES:
            self._feature = value
        else:
            raise TypeError(f"feature must be a string from {ACCEPTED_FEATURES}, not '{value}'")
    
    @property
    def windowSize(self):
        return self._windowSize
    
    @windowSize.setter
    def windowSize(self, value):
        if not isinstance(value, int) or value <= 0:
            raise ValueError("windowSize must be a positive integer.")
        self._windowSize = value
    
    @property
    def genomeFile(self):
        return self._genomeFile
    
    @genomeFile.setter
    def genomeFile(self, value):
        if value is None:
            self._genomeFile = None
        elif os.path.isfile(value):
            self._genomeFile = value
            self.genomeLengths = value
        else:
            raise FileNotFoundError(f"Genome file '{value}' does not exist.")
    
    @property
    def genomeLengths(self):
        if self.genomeFile is None:
            raise ValueError("Genome not loaded. Please set genomeFile before accessing .genomeLengths.")
        
        if self._genomeLengths != None:
            return self._genomeLengths
        else:
            genomeRecords = SeqIO.parse(open(self.genomeFile, "r"), "fasta")
            self._genomeLengths = { record.id:len(record) for record in genomeRecords }
            return self._genomeLengths
    
    @property
    def vcf(self):
        return self._vcf
    
    @vcf.setter
    def vcf(self, value):
        if not hasattr(value, "isVCFTopia") or not value.isVCFTopia:
            raise TypeError("vcf must be a VCFTopia object.")
        self._vcf = value
    
    @property
    def gff3(self):
        return self._gff3
    
    @gff3.setter
    def gff3(self, value):
        if not hasattr(value, "isGFF3Topia") or not value.isGFF3Topia:
            raise TypeError("gff3 must be a GFF3Topia object.")
        
        if value.ncls == None:
            value.create_ncls_index(typeToIndex=["gene"])
        self._gff3 = value
    
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
    def variantFunction(self):
        if self.statistic == "snpnumber":
            return lambda x: 1
        elif self.statistic == "mac":
            return VCFTopia.calculate_mac
        elif self.statistic == "maf":
            return VCFTopia.calculate_maf
        elif self.statistic == "callrate":
            return VCFTopia.calculate_callrate
        elif self.statistic == "het":
            return VCFTopia.calculate_heterozygosity
        else:
            raise ValueError(f"self.statistic=='{self.statistic}' is not handled by .variantFunction")
    
    @property
    def windowFunction(self):
        if self.statistic == "snpnumber":
            return np.sum
        elif self.statistic == "mac":
            return np.mean
        elif self.statistic == "maf":
            return np.mean
        elif self.statistic == "callrate":
            return np.mean
        elif self.statistic == "het":
            return np.mean
        else:
            raise ValueError(f"self.statistic=='{self.statistic}' is not handled by .windowFunction")
    
    @property
    def dtype(self):
        if self.statistic == "snpnumber":
            return np.uint32
        elif self.statistic == "mac":
            return np.float64
        elif self.statistic == "maf":
            return np.float64
        elif self.statistic == "callrate":
            return np.float64
        elif self.statistic == "het":
            return np.float64
        else:
            raise ValueError(f"self.statistic=='{self.statistic}' is not handled by .dtype")
    
    @property
    def width(self):
        if self._width is None:
            return Plot.STANDARD_DIMENSION
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
            return Plot.STANDARD_DIMENSION
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
        raise NotImplementedError("plot() must be implemented in subclasses")

class GenesPlot(Plot):
    def __init__(self, statistic, feature, windowSize, vcf, gff3, genomeFile=None,
                 width=None, height=None):
        super().__init__(statistic, feature, windowSize, vcf, gff3, genomeFile, width, height)
    
    def gene(self, geneID, structure=["CDS"]):
        '''
        Obtains the relevant structural coordinates of the specified gene.
        
        Parameters:
            geneID -- a string indicating the gene ID
            structure -- a list of strings indicating the gene structure components
                         to consider; defaults to ["CDS"] which includes only coding sequences
                         but can be set to include other components like "intron", and 
                         "CDS" can be replaced with "exon" to include all exons.
        Returns:
            contig -- a string indicating the contig on which the gene is located
            coordinates -- a list of tuples indicating the start and end positions
                           of the relevant structural components within the gene.
        '''
        ACCEPTED_STRUCTUREs = ["CDS", "intron", "exon"]
        for struct in structure:
            if struct not in ACCEPTED_STRUCTUREs:
                raise ValueError(f"structure must be a list of strings from {ACCEPTED_STRUCTUREs}, not '{structure}'")
        if "CDS" in structure and "exon" in structure:
            raise ValueError("Cannot include both 'CDS' and 'exon' in structure.")
        
        geneFeature = self.gff3.features[geneID]
        mrnaFeature = GFF3Topia.longest_isoform(geneFeature)
        if "CDS" in structure:
            if not hasattr(mrnaFeature, "CDS"):
                raise ValueError(f"Gene '{geneID}' longest isoform '{mrnaFeature.ID}' does not have CDS features.")
        elif "exon" in structure or "intron" in structure: # exons define introns
            if not hasattr(mrnaFeature, "exon"):
                raise ValueError(f"Gene '{geneID}' longest isoform '{mrnaFeature.ID}' does not have exon features.")
        
        # Get the coordinates within which to locate variants
        if len(structure) == 2: # "intron" and ("CDS" or "exon")
            start, end = mrnaFeature.start_and_end("CDS" if "CDS" in structure else "exon")
            coordinates = [(start, end)]
        else:
            if "intron" in structure: # len(structure) == 1 if we entered this else condition
                exonCoordinates = [ (exonFeature.start, exonFeature.end) for exonFeature in mrnaFeature.exon ]
                exonCoordinates.sort(key=lambda x: x[0]) # sort by start position
                
                coordinates = [ (exonCoordinates[i-1][1]+1, exonCoordinates[i][0]-1) for i in range(1, len(exonCoordinates)-1) ]
            else:
                coordinates = [ (subFeature.start, subFeature.end) for subFeature in getattr(mrnaFeature, structure[0]) ]
        
        return mrnaFeature.contig, coordinates
    
    def variants(self, contig, coordinates):
        '''
        Returns variants overlapping the specified coordinates.
        
        Parameters:
            contig -- a string indicating the contig on which the gene is located
            coordinates -- a list of tuples indicating the start and end positions
                           of the relevant structural components within the gene;
                           as obtained from the gene() method.
        Returns:
            geneVariants -- a list of variants overlapping the specified coordinates.
        '''
        geneVariants = []
        for start, end in coordinates:
            geneVariants.extend(self.vcf.query(contig, (start, end)))
        return geneVariants
    
    def get_x_y(self, geneID):
        '''
        Returns the x list and y array needed for broken_barh plotting of variant/window
        statistics. Method of operation is dictated by self.statistic.
        Returned objects have a length equal to the number of windows needed to cover
        the gene length, where the y array values indicate the individual (self.windowSize==1)
        or averaged (self.windowSize>1) variant statistic
        
        Parameters:
            geneID -- a string indicating the gene ID
        Returns:
            x -- a list containing tuples of (x, 1) where x is the position index
                 and 1 is the length of the position as needed by broken_barh
            y -- a np.array containing the statFunction transformed values for
                 each variant/window
        '''
        # Get gene coordinates
        contig, coordinates = self.gene(geneID)
        
        # Lay out the y values in a numpy array
        y = []
        for coordinate in coordinates:
            coordsArray = np.zeros(coordinate[1] - coordinate[0] + 1, dtype=self.dtype)
            variants = self.variants(contig, [coordinate])
            for variant in variants:
                position = variant.POS - coordinate[0]
                coordsArray[position] = self.variantFunction(variant)
            y.append(coordsArray)
        y = np.concatenate(y)
        
        # Bin the y values to the specified window size
        if self.windowSize > 1:
            y = Plot.bin_values(y, self.windowSize, summariseMethod=self.windowFunction)
        
        # Create the x array
        x = [ (i, 1) for i in np.arange(0, len(y)) ]
        
        return x, y
    
    def plot(self, outputFileName, idsToPlot=None, typesToPlot=["gene"]):
        '''
        Plots the specified .statistic for all genes in the GFF3 file (if 
        idsToPlot is None) or just for the indicated idsToPlot.
        
        Parameters:
            outputFileName -- a string indicating the location to write plot output to;
                              assumes file name is pre-validated to ensure no overwriting
                              and to have a valid filename suffix to control the file format
            idsToPlot -- a list of strings, corresponding to GFF3 features that should
                         be plotted; OR None to plot every feature in self.gff3.ftypes
            typesToPlot -- a list of strings, corresponding to feature types indexed
                           in self.gff3.ftypes, which are to be presented herein.
        '''
        SPACING = 0.1
        
        # Get gene ids for plotting
        if idsToPlot == None:
            idsToPlot = [ f.ID for ftype in typesToPlot for f in self.gff3.ftypes[ftype] ]
        
        # Obtain data for plotting
        geneArrays = [] # contains [[x1, y1], [x2, y2], ...]
        for geneID in idsToPlot:
            geneArrays.append(self.get_x_y(geneID))
        
        # Sort data from longest to shortest
        geneArrays.sort(key = lambda x: len(x[0])) # [x1,y1] each have same length
        
        # Obtain the minimum and maximum values being plotted
        minValue, maxValue, maxLen = np.inf, -np.inf, -np.inf
        for x, y in geneArrays:
            geneMin, geneMax = np.min(y), np.max(y)
            if geneMin < minValue:
                minValue = geneMin
            if geneMax > maxValue:
                maxValue = geneMax
            if len(y) > maxLen:
                maxLen = len(y)
        
        # Establish colour map
        norm = matplotlib.colors.Normalize(vmin=minValue, vmax=maxValue)
        
        # Configure plot
        fig = plt.figure(figsize=(self.width, self.height), tight_layout=True)
        ax = plt.axes()
        
        # Configure x-axis labels
        if self.windowSize == 1:
            ax.set_xlabel(f"Length (bp)", fontweight="bold")
        else:
            divisor, unit = Plot.INTERVALS[10 ** (len(str(self.windowSize)) - 1)]
            ax.set_xlabel(f"Window number ({self.windowSize} {unit} step and width)", fontweight="bold")
        ax.set_xlim(0, maxLen) # maxLen implicitly adds +1 which includes the last bin
        
        # Plot each gene
        ongoingCount = 0.5 # this centers the contig label
        for x, y in geneArrays:
            ax.broken_barh(x, (ongoingCount+SPACING, 1-(SPACING*2)), facecolors=self.cmap(norm(y)))
            ongoingCount += 1
        ax.set_ylim(0.5+SPACING, ongoingCount-SPACING)
        
        # Indicate y-axis labels
        ax.set_ylabel("") # no y-axis title label
        ax.set_yticks(range(1, len(idsToPlot)+1))
        ax.set_yticklabels(idsToPlot)
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical",
                     label="SNP number" if self.statistic == "snpnumber" else \
                           "Minor Allele Count" if self.statistic == "mac" else \
                           "Minor Allele Frequency" if self.statistic == "maf" else \
                           "Genotype Call Rate" if self.statistic == "callrate" else \
                           "Heterozygosity" if self.statistic == "het" else \
                           "labeltypenothandled")
        
        # Save output file
        plt.savefig(outputFileName)
        plt.close()

class ChromosomesPlot(Plot):
    def __init__(self, statistic, feature, windowSize, vcf, gff3=None, genomeFile=None,
                 width=None, height=None):
        super().__init__(statistic, feature, windowSize, vcf, gff3, genomeFile, width, height)
    
    def chromosome(self, contig):
        '''
        Obtains the start and end coordinates of the specified chromosome/contig.
        
        Parameters:
            contig -- a string indicating the chromosome/contig ID
        Returns:
            coordinates -- a list containing a single tuple of the start and end positions
        '''
        contigLength = self.genomeLengths[contig]
        return [(1, contigLength)]
    
    def variants(self, contig, coordinates):
        '''
        Returns variants overlapping the specified chromosome.
        
        Parameters:
            contig -- a string indicating the contig on which the gene is located
            coordinates -- a list of tuples indicating the start and end positions
                           of any regions within the chromosome.
        Returns:
            geneVariants -- a list of variants overlapping the specified coordinates.
        '''
        geneVariants = []
        for start, end in coordinates:
            geneVariants.extend(self.vcf.query(contig, (start, end)))
        return geneVariants
    
    def get_x_y(self, contig):
        '''
        Returns the x list and y array needed for broken_barh plotting of variant/window
        statistics. Method of operation is dictated by self.statistic.
        Returned objects have a length equal to the number of windows needed to cover
        the chromosome/contig length, where the y array values indicate the individual
        (self.windowSize==1) or averaged (self.windowSize>1) variant statistic
        
        Parameters:
            contig -- a string indicating the chromosome/contig ID
        Returns:
            x -- a list containing tuples of (x, 1) where x is the position index
                 and 1 is the length of the position as needed by broken_barh
            y -- a np.array containing the statFunction transformed values for
                 each variant/window
        '''
        # Get contig ene coordinates
        coordinates = self.chromosome(contig)
        
        # Lay out the y values in a numpy array
        y = []
        for coordinate in coordinates:
            coordsArray = np.zeros(coordinate[1] - coordinate[0] + 1, dtype=self.dtype)
            variants = self.variants(contig, [coordinate])
            for variant in variants:
                position = variant.POS - coordinate[0]
                coordsArray[position] = self.variantFunction(variant)
            y.append(coordsArray)
        y = np.concatenate(y)
        
        # Bin the y values to the specified window size
        if self.windowSize > 1:
            y = Plot.bin_values(y, self.windowSize, summariseMethod=self.windowFunction)
        
        # Create the x array
        x = [ (i, 1) for i in np.arange(0, len(y)) ]
        
        return x, y
    
    def plot(self, outputFileName, idsToPlot=None):
        '''
        Plots the specified .statistic for all contigs in the genome file (if 
        idsToPlot is None) or just for the indicated idsToPlot.
        
        Parameters:
            outputFileName -- a string indicating the location to write plot output to;
                              assumes file name is pre-validated to ensure no overwriting
                              and to have a valid filename suffix to control the file format
            idsToPlot -- a list of strings, corresponding to chromosomes/contigs that should
                         be plotted; OR None to plot every contig indexed by self.genomeLengths
        '''
        SPACING = 0.1
        
        # Get contig ids for plotting
        if idsToPlot == None:
            idsToPlot = [ contig for contig in self.genomeLengths.keys() ]
        
        # Obtain data for plotting
        contigArrays = [] # contains [[x1, y1], [x2, y2], ...]
        for contigID in idsToPlot:
            contigArrays.append(self.get_x_y(contigID))
        
        # Sort data from longest to shortest
        contigArrays.sort(key = lambda x: len(x[0])) # [x1,y1] each have same length
        
        # Obtain the minimum and maximum values being plotted
        minValue, maxValue, maxLen = np.inf, -np.inf, -np.inf
        for x, y in contigArrays:
            geneMin, geneMax = np.min(y), np.max(y)
            if geneMin < minValue:
                minValue = geneMin
            if geneMax > maxValue:
                maxValue = geneMax
            if len(y) > maxLen:
                maxLen = len(y)
        
        # Establish colour map
        norm = matplotlib.colors.Normalize(vmin=minValue, vmax=maxValue)
        
        # Configure plot
        fig = plt.figure(figsize=(self.width, self.height), tight_layout=True)
        ax = plt.axes()
        
        # Configure x-axis labels
        if self.windowSize == 1:
            ax.set_xlabel(f"Length (bp)", fontweight="bold")
        else:
            divisor, unit = Plot.INTERVALS[10 ** (len(str(self.windowSize)) - 1)]
            ax.set_xlabel(f"Window number ({self.windowSize} {unit} step and width)", fontweight="bold")
        ax.set_xlim(0, maxLen) # maxLen implicitly adds +1 which includes the last bin
        
        # Plot each chromosome
        ongoingCount = 0.5 # this centers the contig label
        for (x, y), contigID in zip(contigArrays, idsToPlot):
            ax.broken_barh(x, (ongoingCount+SPACING, 1-(SPACING*2)), facecolors=self.cmap(norm(y)))
            ongoingCount += 1
            
            # Plot gene locations if applicable
            if self.gff3 != None:
                geneFeatures = self.gff3.ncls_finder(0, self.genomeLengths[contigID], "contig", contigID)
                mrnaFeatures = [
                    GFF3Topia.longest_isoform(geneFeature)
                    for geneFeature in geneFeatures
                    if hasattr(geneFeature, "mRNA")
                ]
                mrnaCentres = [
                    int((mrnaFeature.start + mrnaFeature.end) / 2) / self.windowSize # convert to fractional window position
                    for mrnaFeature in mrnaFeatures
                ]
                ax.scatter(mrnaCentres, [ongoingCount+SPACING]*len(mrnaCentres), marker="v", color="#EA8527")
            
                ongoingCount += 1
        ax.set_ylim(0.5+SPACING, ongoingCount-SPACING)
        
        # Indicate y-axis labels
        ax.set_ylabel("") # no y-axis title label
        ax.set_yticks(range(1, len(idsToPlot)+1))
        ax.set_yticklabels(idsToPlot)
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical",
                     label="SNP number" if self.statistic == "snpnumber" else \
                           "Minor Allele Count" if self.statistic == "mac" else \
                           "Minor Allele Frequency" if self.statistic == "maf" else \
                           "Genotype Call Rate" if self.statistic == "callrate" else \
                           "Heterozygosity" if self.statistic == "het" else \
                           "labeltypenothandled")
        
        # Save output file
        plt.savefig(outputFileName)
        plt.close()

class MSAPlot(Plot):
    def __init__(self, statistic, windowSize, vcf, gff3=None, genomeFile=None,
                 width=None, height=None):
        super().__init__(statistic, "genes", windowSize, vcf, gff3, genomeFile, width, height)
    
    def plot(self):
        '''
        This method should be implemented to create a plot of the specified statistic
        for the chromosomes in the genome file.
        '''
        raise NotImplementedError("plot() not yet implemented in MSAPlot subclass")
