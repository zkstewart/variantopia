import os, sys
import numpy as np
import plotly.graph_objects as go
from typing import Union

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from vcf import VCFTopia

def isOverlapping(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return end1 >= start2 and end2 >= start1

def simplify_depth(xy: np.ndarray, tolerance: Union[int, float]) -> np.ndarray:
    '''
    Simplify depth data for plotting by removing data points
    that are redundant or minimally informative.
    
    Algorithm works by obtaining the first and last points, in addition to
    any points that mark a departure from the previously stored point by more
    than the specified tolerance.
    
    Parameters:
        xy -- a 2D numpy array where the first column is the x-axis position
              and the second column is the y-axis value
        tolerance -- an integer or float value indicating the maximum increase
                     or decrease in y-axis value before a change is detected
                     and the point is retained
    Returns:
        simplifiedXY -- a 2D numpy array of the input data with some points removed
                        where deemed to be redundant or minimally informative
    '''
    def process_plateau(simplifiedXY, xPlateau, yPlateau):
        # Process plateaus with only one point
        if len(xPlateau) == 1:
            simplifiedXY.append([xPlateau[0], yPlateau[0]])
        
        # Process plateaus with multiple points
        else:
            # Obtain and store the median value of the plateau
            medianY = np.median(yPlateau)
            
            # Store the plateau's median value
            simplifiedXY.append([xPlateau[0], medianY])
            simplifiedXY.append([x-1, medianY])
    
    # Skip simplification if there are too few points
    if len(xy) < 3:
        return xy
    
    # Note first data points
    simplifiedXY = []
    xPlateau = [xy[0][0]]
    yPlateau = [xy[0][1]]
    
    # Begin iteration through remaining data points
    for i in range(1, len(xy)):
        x, y = xy[i]
        
        # Check for departure from this plateau's tolerance threshold
        if y > yPlateau[0] + tolerance or y < yPlateau[0] - tolerance:
            process_plateau(simplifiedXY, xPlateau, yPlateau)
            
            # Begin a new plateau at the point of departure
            xPlateau = [x]
            yPlateau = [y]
        
        # Otherwise, continue to build the plateau
        else:
            xPlateau.append(x)
            yPlateau.append(y)
    
    # Process the last plateau
    process_plateau(simplifiedXY, xPlateau, yPlateau)
    
    # Return the simplified data as a numpy array
    return np.array(simplifiedXY, dtype=float)

def copynum_plot(args):
    vcf = VCFTopia(args.vcfFile)
    
    # Store variant information if a sample deviates from the majority
    sampleArrays = { s:[[], []] for s in vcf.samples }
    for variant in vcf:
        homHet = VCFTopia.count_homhet(variant, countAbsence=True)
        sampleIntegers = VCFTopia.sample_integers(vcf.samples, variant)
                
        if homHet["het"] > homHet["hom"]: # homozygotes are interesting here
            for sample, integers in sampleIntegers.items():
                intSet = set(integers)
                if len(intSet) == 1: # is a homozygote or deletion
                    value = 0 if -1 in integers else 1
                    sampleArrays[sample][0].append(variant.POS)
                    sampleArrays[sample][1].append(value)
        else: # heterozygotes are interesting here            
            for sample, integers in sampleIntegers.items():
                intSet = set(integers)
                if set(integers) != {0}: # is a heterozygote or deletion
                    value = 0 if -1 in integers else 2
                    sampleArrays[sample][0].append(variant.POS)
                    sampleArrays[sample][1].append(value)
    
    # Build plot object
    fig = go.Figure()
    for sample, (x, y) in sampleArrays.items():
        # Initial simplification to remove data redundancy
        xy = np.column_stack((x, y))
        xy = simplify_depth(xy, 0)
        smoothedX, smoothedY = xy[:, 0], xy[:, 1]
        
        # Smooth over genomic window size
        if args.windowSize > 1:
            flatX, flatY = [], []
            for xstart in range(int(smoothedX[0]), int(smoothedX[-1])+1, args.windowSize):
                xend = xstart + args.windowSize
                regionXY = np.array([
                    (_x, _y)
                    for _x, _y in xy
                    if isOverlapping(_x, _x, xstart, xend)
                ])
                if len(regionXY) == 0:
                    continue
                
                # Calculate the median
                _, medianY = np.median(regionXY, axis=0)
                medianY = round(medianY) # naturally moves value to the "poles" of 2 or 0 if inbetween
                
                # Store values
                flatX.extend([xstart, xend])
                flatY.extend([medianY, medianY])
            
            # Re-simplify after transform
            xy = np.column_stack((flatX, flatY))
            xy = simplify_depth(xy, 0)
            smoothedX, smoothedY = xy[:, 0], xy[:, 1]
        
        # Add sample data to interactive plot
        fig.add_trace(go.Scatter(x=smoothedX, y=smoothedY,
                                 mode="lines",
                                 name=sample,
                                 opacity=0.5,
                                 showlegend=True))
    
    # Write output
    fig.write_html(args.outputFileName)

if __name__ == "__main__":
    pass
