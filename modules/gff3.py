#! python3

import re, os, gzip, codecs
from pandas import Series
from contextlib import contextmanager
from ncls import NCLS

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        f.close()
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            f.close()
            return "utf-16"
        except UnicodeDecodeError:
            print(f"Can't tell what codec '{fileName}' is!!")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

class GFF3:
    @staticmethod
    def clean_attributes(attributes):
        '''
        Cleans the attributes string from a GFF3 file by removing unnecessary
        characters and ensuring it is properly formatted.
        
        Parameters:
            attributes -- a string containing the attributes from a GFF3 line.
        Returns:
            cleanedAttributes -- a cleaned string with unnecessary characters removed.
        '''
        return attributes.strip("\r\n\t;'\" ")
    
    def __init__(self, fileLocation, featureTypes=["gene"]):
        self.fileLocation = fileLocation
        self.featureTypes = [ ft.lower() for ft in featureTypes ]
        
        self.features = {}
        
        self.idRegex = re.compile(r"(^|;)ID=(.+?)(;|$)")
        
        self.nclsObj = {}
        self.nclsValues = {}
        
        self.parse_gff3(self.fileLocation)
    
    @property
    def fileLocation(self):
        return self._fileLocation
    
    @fileLocation.setter
    def fileLocation(self, value):
        if not isinstance(value, str):
            raise ValueError("GFF3 file must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"GFF3 file '{value}' is not a file")
        
        self._fileLocation = value
    
    def parse_gff3(self, gff3File):
        '''
        Parses a GFF3 file and populates the graph with features.
        
        Parameters:
            gff3 -- a GFF3 object to parse and populate this graph with.
        '''
        # Reset the object
        self.fileLocation = gff3File
        self.features = {}
        
        # Parse the GFF3 file into a graph structure
        lineCount = 0
        with read_gz_file(self.fileLocation) as fileIn:
            for line in fileIn:
                lineCount += 1
                sl = line.strip("\r\n\t;'\" ").split("\t")
                
                # Skip filler and comment lines
                if line.startswith("#") or len(sl) != 9:
                    continue
                
                # Extract information from this line
                contig, source, ftype, start, end, \
                    score, strand, frame, attributes = sl
                
                # Skip features that are not in the specified feature types
                if ftype.lower() not in self.featureTypes:
                    continue
                
                # Get the ID attribute
                featureID = [ x[1] for x in self.idRegex.findall(GFF3.clean_attributes(attributes)) ] # x == [startCharacter, ID, endCharacter]
                if len(featureID) == 1:
                    featureID = featureID[0]
                elif len(featureID) == 0:
                    featureID = f"{ftype}.{len(self.ftypes[ftype]) + 1}"
                else:
                    raise ValueError(f"GFF3 parsing failed since line #{lineCount} (\"{line}\") has multiple IDs")
                
                # Index the feature
                self.features.setdefault(contig, {})
                self.features[contig][featureID] = {
                    "start": int(start),
                    "end": int(end)
                }
    
    def create_ncls_index(self):
        '''
        Creates an indexed NCLS structure that can be used to find range overlaps
        of indexed features.
        
        Sets:
            self.nclsObj -- a dictionary whre keys are contig names, and values
                            are ncls.NCLS data structures
            self.nclsValues -- a dictionary where keys are contig names, and values
                               are feature IDs in a list with ordering such that
                               returned NCLS integers are the index of the corresponding
                               feature ID.
        '''
        self.nclsObj = {}
        self.nclsValues = {}
        
        for contig, features in self.features.items():
            self.nclsValues[contig] = []
            
            starts, ends, ids = [], [], []
            ongoingCount = 0
            for featureID, valuesDict in features.items():
                starts.append(valuesDict["start"])
                ends.append(valuesDict["end"] + 1) # NCLS indexes 0-based like a range so +1 to make this more logically compliant with gff3 1-based system
                ids.append(ongoingCount)
                self.nclsValues[contig].append(featureID)
                ongoingCount += 1
            
            # Build the NCLS object
            starts = Series(starts)
            ends = Series(ends)
            ids = Series(ids)
            self.nclsObj[contig] = NCLS(starts.values, ends.values, ids.values)
    
    def ncls_finder(self, contig, start, stop):
        '''
        Queries the NCLS structure to find Features that exist within the given
        start->stop range. Specifying the field and value will narrow results
        to only those that have a Feature .field with an equal (==) value.
        
        Parameters:
            contig -- a string indicating the 
            start -- an integer indicating the start position of the feature to check
                     for overlaps
            stop -- an integer indicating the end positon of the feature to check for
                    overlaps; this should be 1-based in GFF3 style e.g., a first
                    position of a feature would be start=1, end=1.
        Returns:
            features -- a list containing Features that overlap the specified range.
                        These Features are NOT deepcopied, so handle them carefully.
        '''
        if not contig in self.nclsObj:
            return []
        
        featureIndices = self.nclsObj[contig].find_overlap(start, stop+1) # Although our ncls is already 1-based, find_overlap acts as a range. We need to +1 to keep everything logically 1-based.
        featureIDs = [ self.nclsValues[i] for i in featureIndices ]
        features = [ self.features[fid] for fid in featureIDs ]
        
        return features
    
    def __iter__(self):
        return iter(self.features.items())
    
    def __contains__(self, value):
        return value in self.features
    
    def keys(self):
        return self.features.keys()
    
    def values(self):
        return self.features.values()
    
    def items(self):
        return self.features.items()
    
    def __repr__(self):
        return "<GFF3 object;file='{0}';num_contigs={1}>".format(
            self.fileLocation,
            len(self.features)
        )
