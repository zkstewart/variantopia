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

class GFF3Topia:
    PARENT_INFERENCE = {
        "CDS": "mRNA",
        "exon": "mRNA",
        "mRNA": "gene",
        "lnc_RNA": "gene",
        "Product": "gene" # Product is a special case, but we treat it as a gene parent
        # "gene": None  # Gene is the top-level feature, no parent should be inferred
    }
    
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
    
    def __init__(self, fileLocation):
        self.fileLocation = fileLocation
        self.ftypes = {}
        self.features = {}
        self.contigs = set()
        
        self.idRegex = re.compile(r"(^|;)ID=(.+?)(;|$)")
        self.parentRegex = re.compile(r"(^|;)Parent=(.+?)(;|$)")
        
        self.ncls = None
        self._nclsType = None
        self._nclsIndex = None
        
        self.parse_gff3(self.fileLocation)
        self.isGFF3Topia = True # flag for easier type checking
    
    @property
    def fileLocation(self):
        return self._fileLocation
    
    @fileLocation.setter
    def fileLocation(self, value):
        if not isinstance(value, str):
            raise ValueError("File location must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"GFF3 file '{value}' is not a file")
        
        self._fileLocation = value
    
    @staticmethod
    def longest_isoform(geneFeature):
        '''
        We pick out the representative gene based on length. If length is identical,
        we'll end up picking the entry listed first in the gff3 file since our > condition
        won't be met. I doubt this will happen much or at all though.
        '''
        if not hasattr(geneFeature, "mRNA"):
            raise ValueError("Longest isoform finding can only occur on features that have mRNA children")
        
        longestMrna = [None, 0]
        for mrnaFeature in geneFeature.mRNA:
            mrnaLen = 0
            
            # Determine the features to use for length calculation
            if hasattr(mrnaFeature, "CDS"):
                features = mrnaFeature.CDS
            elif hasattr(mrnaFeature, "exon"):
                features = mrnaFeature.exon
            else:
                features = []
            
            # Sum the lengths of the CDS (or exon) features
            for subFeature in features:
                mrnaLen += (subFeature.end - subFeature.start + 1)
            
            # Update the longest mRNA if this one is longer
            if mrnaLen > longestMrna[1]:
                longestMrna = [mrnaFeature, mrnaLen]
        return longestMrna[0]
    
    def _get_unique_feature_id(self, inputID):
        ongoingCount = 1
        featureID = inputID
        while featureID in self.features:
            featureID = f"{inputID}.{ongoingCount}"
            if not featureID in self.features:
                break
            ongoingCount += 1
        return featureID
    
    def parse_gff3(self, gff3File):
        '''
        Parses a GFF3 file and populates the graph with features.
        
        Parameters:
            gff3 -- a GFF3 object to parse and populate this graph with.
        '''
        # Reset the graph
        self.fileLocation = gff3File
        self.ftypes = {}
        self.features = {}
        self.contigs = set()
        
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
                start = int(start)
                end = int(end)
                ftype = GFF3Feature.make_ftype_case_appropriate(ftype)
                attributes = GFF3Topia.clean_attributes(attributes) # necessary for regex to work properly
                
                # Establish or populate tracking containers
                self.ftypes.setdefault(ftype, [])
                self.contigs.add(contig)
                
                # Get the ID attribute
                featureID = [ x[1] for x in self.idRegex.findall(attributes) ] # x == [startCharacter, ID, endCharacter]
                if len(featureID) == 1:
                    featureID = featureID[0]
                elif len(featureID) == 0:
                    featureID = f"{ftype}.{len(self.ftypes[ftype]) + 1}"
                else:
                    raise ValueError(f"GFF3 parsing failed since line #{lineCount} (\"{line}\") has multiple IDs")
                
                # Get the parent ID(s)
                parentIDs = [ x[1] for x in self.parentRegex.findall(attributes) ] # x == [startCharacter, ID, endCharacter]
                
                # Create a feature object
                feature = GFF3Feature(ID=featureID, ftype=ftype,
                                      start=start, end=end, strand=strand,
                                      contig=contig, children=[], parents=parentIDs)
                
                # Index the feature if it doesn't already exist
                if not featureID in self.features:
                    self.add(feature)
                # Specifically handle exons or CDS which are allowed to have duplicated IDs
                elif ftype in ["exon", "CDS"]:
                    feature.ID = self._get_unique_feature_id(featureID)
                    self.add(feature)
                # Handle other duplicated feature types
                else:
                    "We assume that the GFF3 is unsorted if we reach this point, so we are detailing an existing feature"
                    feature = self.features[featureID]
                    
                    # Check that inferred details are correct
                    if feature.ftype != ftype:
                        raise ValueError(f"Unsorted GFF3 issue: Feature ID '{featureID}' has a different type '{feature.ftype}' than previously inferred '{ftype}'")
                    if feature.contig != contig:
                        raise ValueError(f"Unsorted GFF3 issue: Feature ID '{featureID}' has a different contig '{feature.contig}' than previously inferred '{contig}'")
                    
                    # Update feature details
                    feature.start = start
                    feature.end = end
                    feature.strand = strand
                    feature.parents.update(parentIDs) # add parents to existing set
                    self.add(feature) # re-add to ensure parents are updated correctly
    
    def add(self, feature):
        # Store feature within the graph
        if not feature.ID in self.features: # only if it doesn't already exist
            self.ftypes.setdefault(feature.ftype, []) # necessary if first occurrence of a subfeature preceeds its parent type
            self.ftypes[feature.ftype].append(feature.ID)
            self.features[feature.ID] = feature
        
        # Update graph features with parent-child relationships
        for parentID in feature.parents:
            # Associate the feature with its existing parents
            if parentID in self.features:
                self.features[parentID].add_child(feature)
            # Create a placeholder for the parent if it doesn't exist
            else:
                if feature.ftype in GFF3Topia.PARENT_INFERENCE:
                    parentFeature = GFF3Feature(parentID, GFF3Topia.PARENT_INFERENCE[feature.ftype],
                                                contig=feature.contig,
                                                children=[feature])
                    self.add(parentFeature)
                else:
                    raise ValueError("Your GFF3 is not sorted in top-down hierarchical order which has caused an error; " +
                                     f"I encountered a {feature.ftype} with ID '{feature.ID}' that has a parent '{parentID}' which has " + 
                                     f"not yet appeared in your GFF3 file. I am unsure what parent type to infer for " +
                                     f"'{feature.ftype}' features, so I cannot continue parsing. Sort your GFF3 file in " +
                                     "conventional top-down hierarchical order before trying again.")
    
    def create_ncls_index(self, typeToIndex=["gene"]):
        '''
        Creates an indexed NCLS structure that can be used to find range overlaps
        for the feature types of interest.
        
        Associates the created index to the .ncls field of this object instance.
        A hidden ._nclsIndex dictionary links the ncls indices to feature objects.
        
        Parameters:
            typeToIndex -- a string (case-sensitive) indicating the entry type
                           to index OR an iterable of strings indicating multiple
                           types to index.
        '''
        if isinstance(typeToIndex, str):
            typeToIndex = [typeToIndex]
        
        for indexType in typeToIndex:
            assert indexType in self.ftypes, \
                "'{0}' not found as a feature type within the parsed GFF3 ('{1}')".format(indexType, self.fileLocation)
        
        nclsIndex = {}
        starts, ends, ids = [], [], []
        
        # Add features of the specified type to our NCLS structure
        ongoingCount = 0
        for indexType in typeToIndex:
            for featureID in self.ftypes[indexType]:
                feature = self.features[featureID]
                starts.append(feature.start)
                ends.append(feature.end + 1) # NCLS indexes 0-based like a range so +1 to make this more logically compliant with gff3 1-based system
                ids.append(ongoingCount)
                nclsIndex[ongoingCount] = feature
                ongoingCount += 1
        
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        
        # Associate it to this instance
        self.ncls = ncls
        self._nclsType = typeToIndex
        self._nclsIndex = nclsIndex
    
    def ncls_finder(self, start, stop, field, value):
        '''
        Queries the NCLS structure to find Features that exist within the given
        start->stop range. Specifying the field and value will narrow results
        to only those that have a Feature .field with an equal (==) value.
        
        Parameters:
            start -- an integer indicating the start position of the feature to check
                     for overlaps
            end -- an integer indicating the end positon of the feature to check for
                   overlaps; this should be 1-based in GFF3 style e.g., a first
                   position of a feature would be start=1, end=1.
            field -- a string (case-sensitive) indicating the field of the Feature
                     object that we want to check. For example, if you want to find
                     features that overlap the given start->stop range on the contig
                     "X", you'd provide "contig" as the field so this function knows
                     to check the Feature.contig field for the value of "X".
            value -- a string (case-sensitive) indicating the value of the Feature
                     field we want to find. As in the above example, if you want to
                     find the value "X" within the .contig field, you'd provide "X" as
                     the value here.
        Returns:
            features -- a list containing Features that overlap the specified range.
                        These Features are NOT deepcopied, so handle them carefully.
        '''
        assert self.ncls != None and self._nclsIndex != None, \
            "Run create_ncls_index before you call this method!"
        
        overlaps = self.ncls.find_overlap(start, stop+1) # Although our ncls is already 1-based, find_overlap acts as a range. We need to +1 to keep everything logically 1-based.
        
        features = []
        for result in overlaps: # result == [start, end, index]
            feature = self._nclsIndex[result[2]]
            if feature.__dict__[field] == value:
                features.append(feature)
        
        # Return list
        return features
    
    def qc(self, typesToCheck=None):
        '''
        Runs a quality control check on the GFF3Topia object to ensure that all
        features are properly linked.
        
        Prints a warning if any features are found that have no parents or children.
        
        Parameters:
            typesToCheck -- an iterable of strings indicating the feature types to check
                            for dangling features. If None, all feature types are checked.
        '''
        danglingFeatures = {}
        for feature in self:
            if typesToCheck == None or feature.ftype in typesToCheck:
                if len(feature.parents) == 0 and len(feature.children) == 0:
                    danglingFeatures.setdefault(feature.ftype, 0)
                    danglingFeatures[feature.ftype] = 1
        
        if len(danglingFeatures) != 0:
            print(f"WARNING: Parsing '{self.fileLocation}' resulted in dangling features with no parents or children, " +
                  "which is likely due to an unsorted or incorrectly formatted GFF3 file. This may cause issues with " +
                  "psQTL's functionality.")
            for ftype, count in danglingFeatures.items():
                print(f"# {count} '{ftype}' feature{'s have' if count > 1 else ' has'} no parents or children")
    
    def __getitem__(self, key):
        return self.features[key]
    
    def __len__(self):
        return len(self.features)
    
    def __iter__(self):
        return iter(self.features.values())
    
    def __contains__(self, item):
        return item.ID in self.features
    
    def has_key(self, key):
        return key in self.features
    
    def keys(self):
        return self.features.keys()
    
    def values(self):
        return self.features.values()
    
    def items(self):
        return self.features.items()
    
    def __repr__(self):
        return "<GFF3Topia object;file='{0}';num_contigs={1};{2}>".format(
            self.fileLocation,
            len(self.contigs),
            ";".join(["num_{0}={1}".format(key, len(self.ftypes[key])) for key in self.ftypes.keys()])
        )


class GFF3Topia:
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
        self.index = {} # allows feature lookup by ID
        self.nclsObj = {}
        self.nclsValues = {}
        
        self.idRegex = re.compile(r"(^|;)ID=(.+?)(;|$)")
        self.parse_gff3(self.fileLocation)
        self.isGFF3Topia = True
    
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
        self.index = {}
        
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
                featureID = [ x[1] for x in self.idRegex.findall(GFF3Topia.clean_attributes(attributes)) ] # x == [startCharacter, ID, endCharacter]
                if len(featureID) == 1:
                    featureID = featureID[0]
                elif len(featureID) == 0:
                    featureID = f"{ftype}.{len(self.ftypes[ftype]) + 1}"
                else:
                    raise ValueError(f"GFF3 parsing failed since line #{lineCount} (\"{line}\") has multiple IDs")
                
                # Index the feature
                self.features.setdefault(contig, {})
                self.features[contig][featureID] = {
                    "id": featureID,
                    "contig": contig,
                    "start": int(start),
                    "end": int(end)
                }
                self.index[featureID] = self.features[contig][featureID]
    
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
        featureIDs = [ self.nclsValues[contig][i] for start, stop, i in featureIndices ] # featureIndices == [[start, end, index], ...]
        features = [ self.features[contig][fid] for fid in featureIDs ]
        
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
        return "<GFF3Topia object;file='{0}';num_contigs={1}>".format(
            self.fileLocation,
            len(self.features)
        )
