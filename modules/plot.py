import os, gzip
from Bio import SeqIO

#from .vcf import VCFTopia
#from .gff3 import GFF3Topia

class Plot:
    Plot.STANDARD_DIMENSION = 5
    
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
        
        if value.nclsObj == {}:
            value.create_ncls_index()
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
    
    def scatter(self, contigID, start, end):
        '''
        Returns data suited for scatter plotting of WindowedNCLS values.
        
        Parameters:
            contigID -- a string indicating the contig ID
            start -- an integer indicating the start position of the region
            end -- an integer indicating the end position of the region
        Returns:
            x -- a numpy array of the x values (positions)
            y -- a numpy array of the y values (statistical values)
        '''
        regionValues = windowedNCLS.find_overlap(contigID, start, end)
        x, y = [], []
        for pos, _, ed in regionValues:
            if clipped and pos < start: # can occur if windowSize > 1; this is treated as our first value
                pos = start
            if clipped and pos > end: # this should not happen
                continue
            x.append(pos)
            y.append(ed)
        x = np.array(x)
        y = np.array(y)
        
        return x, y
    
    
    
    def plot(self):
        raise NotImplementedError("plot() must be implemented in subclasses")

class GenesPlot(Plot):
    def __init__(self, statistic, feature, windowSize, vcf, gff3, genomeFile=None,
                 width=None, height=None):
        super().__init__(statistic, feature, windowSize, vcf, gff3, genomeFile, width, height)
    
    def variants(self, geneID):
        '''
        Returns variants overlapping the specified gene.
        
        Parameters:
            geneID -- a string indicating the gene ID
        Returns:
            x -- a numpy array of the x values (positions)
            y -- a numpy array of the y values (statistical values)
        '''
        geneFeature = self.gff3.index[geneID]
        geneVariants = self.vcf.query(geneFeature["id"], (geneFeature["start"], geneFeature["end"]))
        
        x, y = [], []
        for pos, _, ed in regionValues:
            if clipped and pos < start: # can occur if windowSize > 1; this is treated as our first value
                pos = start
            if clipped and pos > end: # this should not happen
                continue
            x.append(pos)
            y.append(ed)
        x = np.array(x)
        y = np.array(y)
        
        return x, y
    # "snpnumber", "mac", "maf", "callrate", "het"
    
    def snpnumber(self, geneID):
        '''
        Returns the number of SNPs overlapping the specified gene. As an implicitly
        windowSize=1 statistic (prior to any binning of this function's output),
        this function returns a np.array with length equal to the gene length,
        where 0 == no SNP and 1 == SNP presence at that position.
        
        Parameters:
            geneID -- a string indicating the gene ID
        Returns:
            x -- a numpy array of the x values (positions)
            y -- a numpy array of the y values (number of SNPs)
        '''
        geneFeature = self.gff3.index[geneID]
        snpCount = self.vcf.snp_count(geneFeature["id"], (geneFeature["start"], geneFeature["end"]))
        
        x = np.array([geneFeature["start"], geneFeature["end"]])
        y = np.array([snpCount, snpCount])
        
        return x, y
    
    def plot(self):
        '''
        This method should be implemented to create a plot of the specified statistic
        for the genes in the GFF3 file.
        '''
        raise NotImplementedError("plot() must be implemented in GenesPlot subclass")

class ChromosomesPlot(Plot):
    def __init__(self, statistic, feature, windowSize, vcf, gff3=None, genomeFile=None,
                 width=None, height=None):
        super().__init__(statistic, feature, windowSize, vcf, gff3, genomeFile, width, height)

    def plot(self):
        '''
        This method should be implemented to create a plot of the specified statistic
        for the chromosomes in the genome file.
        '''
        raise NotImplementedError("plot() must be implemented in ChromosomesPlot subclass")
