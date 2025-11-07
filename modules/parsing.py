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

class WriteGzFile:
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
