import os, gzip
from Bio import SeqIO

from .vcf import VCFTopia

class Plot:
    def __init__(self, statistic, feature, windowSize, genomeFile=None, gff3File=None,
                 width=None, height=None):
        '''
        Parameters:
            statistic -- a string indicating the statistic to plot; currently supported
                         values are "snpnumber", "mac", "maf", "callrate", and "het"
            feature -- a string indicating the feature to plot; currently supported
                       values are "genes" and "chromosomes"
            windowSize -- a positive integer indicating the size of the window to use
                          for plotting
            gff3File -- a string pointing to a GFF3 file for use in plotting
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
        self.genomeFile = genomeFile
        self.gff3File = gff3File
        
        # Aesthetic parameters
        self.width = width
        self.height = height
        
        # Figure-related parameters (not to be set by user)
        self.fig = None
        self.axs = None
        self.rowNum = None
    
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
            raise ValueError("Genome not loaded. Please set genomeFile before accessing genomeLengths.")
        
        if self._genomeLengths != None:
            return self._genomeLengths
        else:
            genomeRecords = SeqIO.parse(open(self.genomeFile, "r"), "fasta")
            self._genomeLengths = { record.id:len(record) for record in genomeRecords }
            return self._genomeLengths
    
    @property
    def gff3File(self):
        return self._gff3File
    
    @gff3File.setter
    def gff3File(self, value):
        if value is None:
            self._gff3File = None
        elif os.path.isfile(value):
            self._gff3File = value
            self.genomeLengths = value
        else:
            raise FileNotFoundError(f"GFF3 file '{value}' does not exist.")
    
    @property
    def gff3(self):
        if self.genomeFile is None:
            raise ValueError("Genome not loaded. Please set genomeFile before accessing genomeLengths.")
        
        if self._genomeLengths != None:
            return self._genomeLengths
        else:
            genomeRecords = SeqIO.parse(open(self.genomeFile, "r"), "fasta")
            self._genomeLengths = { record.id:len(record) for record in genomeRecords }
            return self._genomeLengths
    
    def plot(self):
        raise NotImplementedError("plot() must be implemented in subclasses")
