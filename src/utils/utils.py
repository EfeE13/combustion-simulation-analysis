import subprocess
from datetime import datetime
from typing import Optional
from scipy.special import erfcinv
import numpy as np

def run(command: str) -> str:
    result = subprocess.run(command, shell = True, capture_output = True, text = True)
    if str(result.stderr) != "":
        raise Exception("Error in " + command + ": " + str(result.stderr))
    return result.stdout

def fixDirFormat(dirNameFromStart:str) -> str:
    """Forces directory name to start with / and removes / from end"""

    if dirNameFromStart[0] != "/":
        dirNameFromStart = "/" + dirNameFromStart
    if dirNameFromStart[-1] == "/":
        dirNameFromStart = dirNameFromStart[:-1]
    
    return dirNameFromStart

def makeDir(dirNameFromStart:str) -> None:
    """Makes directory (including any required, but not present, parent folders)"""

    dirNameFromStart = fixDirFormat(dirNameFromStart)
    
    directoryToMake = ""
    for i in range(1, len(dirNameFromStart.split("/"))):
        directoryToMake += "/" + dirNameFromStart.split("/")[i]
        subprocess.run("mkdir " + directoryToMake, shell = True, capture_output = True, text = True)

def outputTime(message:Optional[str] = "") -> str:
    return str(message) + ", " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S")

def dayAndTime() -> str:
    return datetime.now().strftime("%Y-%m-%dAt%H-%M-%S").strip()

def er_func(Z:float) -> float:
    return float(np.exp(-2 * erfcinv(2*Z)**2))

def log_modulus(x_list):
    return np.sign(x_list) * np.log1p(np.abs(x_list)) / np.log(10)