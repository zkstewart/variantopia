# Copyright (C) 2026 Zachary Kenneth Stewart

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os, shutil, subprocess

def run_tabix(vcfFile, tabixExe=None):
    '''
    Create an index for the VCF file using the tabix executable.
    
    Parameters:
        vcfFile -- a string pointing to a VCF file
        tabixExe -- path to the tabix executable
    '''
    if tabixExe == None:
        tabixExe = shutil.which("tabix")
    if ( tabixExe == None ) or ( not os.path.exists(tabixExe) ):
        errorMsg = f"tabix executable not found at '{tabixExe}'" if tabixExe != None \
                   else "tabix executable not locateable through system PATH variable"
        raise FileNotFoundError(errorMsg)
    
    # Run tabix to index the VCF file
    cmd = [tabixExe, "-f", "-p", "vcf", vcfFile]
    run_tabix = subprocess.Popen(" ".join(cmd), shell = True,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE)
    tabixout, tabixerr = run_tabix.communicate()
    
    # Check that tabix ran successfully
    if tabixout.decode("utf-8") != "":
        raise Exception(("tabix stdout is not empty as expected, indicating a probable error. " +
                        f'Please check the stdout ({tabixout.decode("utf-8")}) and stderr ' + 
                        f'({tabixerr.decode("utf-8")}) to make sense of this.'))
    elif tabixerr.decode("utf-8") != "":
        raise Exception(("tabix encountered an error; have a look " +
                        f'at the stdout ({tabixout.decode("utf-8")}) and stderr ' + 
                        f'({tabixerr.decode("utf-8")}) to make sense of this.'))

def run_bcftools_index(vcfFile, bcftoolsExe=None):
    '''
    Create an index for the VCF file using the bcftools executable.
    
    Parameters:
        vcfFile -- a string pointing to a VCF file
        bcftoolsExe -- path to the bcftools executable
    '''
    if bcftoolsExe == None:
        bcftoolsExe = shutil.which("bcftools")
    if ( bcftoolsExe == None ) or ( not os.path.exists(bcftoolsExe) ):
        errorMsg = f"bcftools executable not found at '{bcftoolsExe}'" if bcftoolsExe != None \
                   else "bcftools executable not locateable through system PATH variable"
        raise FileNotFoundError(errorMsg)
    
    # Run tabix to index the VCF file
    cmd = [bcftoolsExe, "index", vcfFile]
    run_index = subprocess.Popen(" ".join(cmd), shell = True,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE)
    indexout, indexerr = run_index.communicate()
    
    # Check that bcftools index ran successfully
    if indexout.decode("utf-8") != "":
        raise Exception(("bcftools index stdout is not empty as expected, indicating a probable error. " +
                        f'Please check the stdout ({indexout.decode("utf-8")}) and stderr ' + 
                        f'({indexerr.decode("utf-8")}) to make sense of this.'))
    elif indexerr.decode("utf-8") != "":
        raise Exception(("bcftools index encountered an error; have a look " +
                        f'at the stdout ({indexout.decode("utf-8")}) and stderr ' + 
                        f'({indexerr.decode("utf-8")}) to make sense of this.'))

UTILITIES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "utilities")

def check_r_packages(packages, rscriptExe):
    '''
    Checks that R packages are installed for the given Rscript executable.

    Parameters:
        packages -- a list of strings of R package names
        rscriptExe -- path to the Rscript executable
    Returns:
        missing -- a list of the packages which are not installed
    '''
    rCode = "for (p in commandArgs(TRUE)) if (!requireNamespace(p, quietly=TRUE)) cat(p, '\\n')"
    result = subprocess.run([rscriptExe, "-e", rCode] + packages,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise Exception(f"Rscript failed to check for R packages; stderr was:\n{result.stderr.decode('utf-8')}")
    return result.stdout.decode("utf-8").split()

def run_rscript(scriptName, args, rscriptExe, logFile=None):
    '''
    Runs an R script from the variantopia utilities directory.

    Parameters:
        scriptName -- a string of the script's file name within utilities
        args -- a list of strings of arguments for the script
        rscriptExe -- path to the Rscript executable
        logFile -- (OPTIONAL) a string of a location to write the script's
                   stdout and stderr to
    '''
    scriptFile = os.path.join(UTILITIES_DIR, scriptName)
    if not os.path.isfile(scriptFile):
        raise FileNotFoundError(f"R script '{scriptFile}' not found")

    cmd = [rscriptExe, scriptFile] + args
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = result.stdout.decode("utf-8")
    if logFile is not None:
        with open(logFile, "w") as fileOut:
            fileOut.write(output)
    if result.returncode != 0:
        raise Exception(f"{scriptName} failed with exit code {result.returncode}; command was " +
                        f"'{' '.join(cmd)}'. Output tail:\n" + "\n".join(output.splitlines()[-15:]))
