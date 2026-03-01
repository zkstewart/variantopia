#! python3

import gzip, codecs, math
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

def get_chunking_points(numberToChunk, chunks, isNumOfChunks=True):
    '''
    This is a general purpose function to take in a number of "things"
    that you want to chunk, and find out how to chunk them evenly.
    
    The resulting list should be interpreted as the 0-based indices where
    a new chunk should form. You should check for this index at the start
    of a loop, and form a new file if your index == the value in this list.
    
    Also, this uses "allocated chunking" such that it will try to keep
    the number of things per chunk approximately equal. Even if you specify
    X number of things per chunk, it might be more optimal to have X-1 in each
    chunk so as to make sure the last chunk doesn't contain a single thing.
    This might not be what you want, but usually, allocated chunking leads to
    more optimal code (e.g., a major use of this function could be for
    parallel processing of the chunks).
    
    Params:
        numberToChunk -- an integer value, possibly derived from a list length as example.
        chunks -- an integer value for the desired number of chunks OR the number of
                  sequences to contain within each chunk, determined by
        isNumOfChunks -- a boolean indicating whether you want the number to be the number
                         of chunks (True), or the number of sequences within each chunk (False)
    '''
    assert isinstance(numberToChunk, int)
    assert isinstance(chunks, int)
    if numberToChunk < chunks:
        raise Exception(f"Chunking only valid if chunkSize <= chunks i.e., {chunks} <= {numberToChunk}")
    
    # Derive how many chunks we want to split the file into
    if isNumOfChunks:
        numChunks = chunks
    else:
        numChunks = math.ceil(numberToChunk / chunks)
    
    rawNum = numberToChunk / numChunks # This line is more relevant in the multithreading code I took this from, but it's okay to just leave it.
    numRoundedUp = round((rawNum % 1) * numChunks, 0) # By taking the decimal place and multiplying it by the num of chunks, we can figure out how many chunks need to be rounded up
    
    # Store positions at which to start a new chunk
    chunkPoints = []
    ongoingCount = 0
    for i in range(numChunks):
        
        # Determine where chunks begin in 0-based indexing
        if i < numRoundedUp: # decide if the number of sequences in this chunk should be rounded up
            point = math.ceil(rawNum) + ongoingCount # Round up the rawNum, and also add our ongoingCount which corresponds to the number of things already put into a chunk
            
            # Prevent chunking beyond the last index where a chunk should start
            if point >= numberToChunk: # Without this check, if we have more chunks than things to chunk, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
                break  # This doesn't actually affect program function, but for aesthetic reasons and for clarity of how this function works, I prevent this from occurring.
            
            chunkPoints.append(point)
            ongoingCount += math.ceil(rawNum)
        else:
            point = math.floor(rawNum) + ongoingCount # Round down the rawNum since we've already accounted for any extra uneven numbers
            
            if point >= numberToChunk:
                break
            
            chunkPoints.append(point)
            ongoingCount += math.floor(rawNum)
    
    return chunkPoints
