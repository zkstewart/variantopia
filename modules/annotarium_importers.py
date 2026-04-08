#! python3

import os, sys

def import_annotarium_gff3(annotariumDir):
    try:
        sys.path.append(os.path.dirname(annotariumDir))
        from annotarium import GFF3Feature, GFF3Tarium
        return GFF3Feature, GFF3Tarium
    except ModuleNotFoundError:
        raise ModuleNotFoundError(f"Could not import GFF3 classes from '{args.annotariumDir}'")

def import_annotarium_domains(annotariumDir):
    try:
        sys.path.append(os.path.dirname(annotariumDir))
        from annotarium import Domains, OverlapResolver
        return Domains, OverlapResolver
    except ModuleNotFoundError:
        raise ModuleNotFoundError(f"Could not import Domains classes from '{args.annotariumDir}'")
