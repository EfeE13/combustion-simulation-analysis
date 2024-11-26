from src.pdrs_analysis.PDRsRunsND import PDRsRunsND
from src.utils.plotting_utils.ColorMapService import ColorMapService
import numpy as np

class PDRsRuns2D(PDRsRunsND):
    """The data for a number of PDRs Runs where only two of the PDRs inputs (specifically, nickname1 and nickname2) are varied. Note it extends PDRsRuns object

    Attributes
    PDRsRunDict: Dict[str, PDRsRun] with key:value pairs as PDRsOutputFileName:PDRsRunObject
    fileInputsDict: Dict[str, Dict[str, float]] with key:value pairs as PDRsOutputFileName:{"SDR":30.0, "FMX":0.1}
    nickname1: str representing name of the first PDRs input variable varied
    nickname2: str representing name of the second PDRs input variable varied
    """

    def __init__(self, outputFilenameList:list[str], nickname1:str, nickname2:str) -> None:
        """Only two of the PDRs inputs (specifically, nickname1 and nickname2) are varied.

        :param list[str] outputFilenameList: list of outputs one can get from PDRsRunnable.py
        :param str nickname1: nickname for 1st input vaired
        :param str nickname2: nickname for 2nd input vaired
        """

        super().__init__(outputFilenameList)
        self.nickname1 = nickname1
        self.nickname2 = nickname2

        self.xVals = []
        self.yVals = []
        for filename in list(self.fileInputsDict.keys()):
            xVal = float(self.fileInputsDict[filename][self.nickname1])
            if not (xVal in self.xVals):
                self.xVals.append(xVal)
            yVal = float(self.fileInputsDict[filename][self.nickname2])
            if not (yVal in self.yVals):
                self.yVals.append(yVal)

    def makeColorPlot(self, axObj, columnHeaderForZVals:str, typeToPlot:str, xAxisLabel:str, yAxisLabel:str, titleLabel:str, colorbarLabel:str, zMinVal:float = float("-inf"), zMaxVal:float = float("inf"), plotLogX:bool = False, plotLogY:bool = False, plotLogZ:bool = False) -> None:
        """Uses ColorMapService from plotting_utils to make color plot colored by some function of columnHeader's values specified by typeToPlot

        :param _type_ axObj: axis object to plot on to
        :param str columnHeaderForZVals: representing the PDRs output quantity that will be used for the coloring
        :param str typeToPlot: if "max" or "average", the maximum or average value of columnHeader column will be used for the coloring. Different values of type will require changing this function (closed for extension, open for modification!)
        :param str xAxisLabel: x-axis string
        :param str yAxisLabel: y-axis string
        :param str titleLabel: title string
        :param str colorbarLabel: colorbar string
        :param float zMinVal: minimum value for color map, defaults to float("-inf")
        :param float zMaxVal: maximum value for color map, defaults to float("inf")
        :param bool plotLogX: to plot x-axis logarithmically, defaults to False
        :param bool plotLogY: to plot y-axis logarithmically, defaults to False
        :param bool plotLogZ: to plot color logarithmically, defaults to False
        """

        def logBase10Helper(val, takeLog):
            if takeLog:
                return float(np.log10(val))
            else:
                return val

        triplesToColorPlot = [] # List of sublists of length 3 [x, y, z] where location (x, y) will be colored by z
        for filename in list(self.fileInputsDict.keys()):
            tripleToColorPlot = [logBase10Helper(self.fileInputsDict[filename][self.nickname1], plotLogX), logBase10Helper(self.fileInputsDict[filename][self.nickname2], plotLogY)]
            colValues = self.PDRsRunsDict[filename].getCol(columnHeaderForZVals)
            if typeToPlot.lower() == "max":
                zVal = max(colValues)
            elif typeToPlot.lower() in ["avg", "average"]:
                zVal = sum(colValues) / len(colValues)
            else:
                raise Exception("Unimplemented type")
            
            zVal = max(zVal, zMinVal)
            zVal = min(zVal, zMaxVal)
            tripleToColorPlot.append(logBase10Helper(zVal, plotLogZ))
            triplesToColorPlot.append(tripleToColorPlot)
        
        ColorMapMaker = ColorMapService(triplesToColorPlot)
        ColorMapMaker.makeColorPlotContinuous(axObj, xAxisLabel, yAxisLabel, titleLabel, colorbarLabel)