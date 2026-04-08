import os, sys

# Make classes accessible to external callers
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.vcftopia import VCFTopia
from modules.vcfplot import GenesPlot, ChromosomesPlot
