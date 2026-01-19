#! python3

import gzip, codecs
from contextlib import contextmanager

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
            print(f"'{fileName}' is neither utf-8 nor utf-16 encoded; please convert to one of these formats.")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

class GzCapableWriter:
    def __init__(self, filename):
        self.filename = filename
        self.file = None
    
    def __enter__(self):
        if self.filename is None:
            return None
        else:
            if self.filename.endswith(".gz"):
                self.file = gzip.open(self.filename, "wt")
            else:
                self.file = open(self.filename, "w")
            return self.file
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file:
            self.file.close()
        if exc_type is not None:
            raise exc_type(exc_val).with_traceback(exc_tb)

def parse_metadata_groups(metadataGroupsFile):
    ACCEPTED_GROUPS = ["0", "1", "2"]
    metadataGroups = {}
    with read_gz_file(metadataGroupsFile) as fileIn:
        for line in fileIn:
            l = line.strip()
            if l != "":
                if "\t" in l:
                    delimiter = "\t"
                elif "," in l:
                    delimiter = ","
                else:
                    raise ValueError(f"Cannot parse metadata groups since I expect a tab or comma delimiter in line '{l}'")
                
                species, group = l.split("\t")
                if not group in ACCEPTED_GROUPS:
                    raise ValueError(f"Group should be in {ACCEPTED_GROUPS}; do not recognise '{group}'")
                
                metadataGroups[species] = group
    return metadataGroups

def parse_2col_tsv_as_dict(tsvFile, keyIsLeft=True):
    '''
    General purpose parser for 2-column TSV to render a dictionary pairing
    one column (key) to the other (value) according to the keyIsLeft behavioural boolean.
    
    Parameters:
        tsvFile -- a string indicating the location of the file to parse
        keyIsLeft -- (OPTIONAL) a boolean controlling whether the left column should be
                     the dictionary key (default; ==True) or if this should be inversed
                     (==False)
    Returns:
        tsvDict -- a dictionary pairing keys and values corresponding to TSV file column contents.
    '''
    tsvDict = {}
    with read_gz_file(tsvFile) as fileIn:
        for line in fileIn:
            delim = "\t" if "\t" in line else ","
            sl = line.strip().split(delim)
            if len(sl) != 2:
                raise ValueError(f"TSV file line should have two columns, but has {len(sl)}; offending line is '{line.rstrip()}'")
            
            if keyIsLeft:
                left, right = sl
            else:
                right, left = sl
            
            tsvDict[left] = right
    return tsvDict
