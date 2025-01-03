from src.pdrs_analysis.PDRsRun import PDRsRun
from src.utils.plotting_utils.RegularPlotService import RegularPlotService
import numpy as np

class PDRsRuns():
    """
    The composition of a number of PDRsRun objects

    Attributes
    PDRsRunDict: Dictionary with key:value pairs as PDRsOutputFileName:PDRsRunObject
    fileInputsDict: Dictionary with key:value pairs as PDRsOutputFileName:{"SDR":30.0, "FMX":0.1}
    """
    
    def __init__(self, outputFilenameList:list[str]) -> None:
        """Makes a basic PDRs run object that the PDRsRuns1D and PDRsRuns2D inherits from. n refers to the number of parameters vaired when running PDRs

        :param list[str] outputFilenameList: represent the PDRs output filenames. They should be something like SDR_30.0__FMX_0.1.Y (or other extension)
        """
        
        self.PDRsRunsDict = {}
        for outputFilename in outputFilenameList:
            self.PDRsRunsDict[outputFilename] = PDRsRun(outputFilename)

    def makeLinePlot(self, axObj, columnHeader:str, typeToPlot:str, yAxisLabel:str, titleLabel:str, yMinVal:float = float("-inf"), yMaxVal:float = float("inf"), plotLogY:bool = False) -> None:
        """Makes a regular line plot of 1D data

        :param axObj: the axis object to be drawn on
        :param str columnHeader: representing the PDRs output quantity that will be used for the plotting
        :param str typeToPlot: if "max" or "average", the maximum or average value of columnHeader column will be used for the coloring. Different values of type will require changing this function (closed for extension, open for modification!)
        :param str yAxisLabel: x-axis string
        :param str titleLabel: y-axis string
        :param float yMinVal: minimum value for plotting, defaults to float("-inf")
        :param float yMaxVal: maximum value for plotting, defaults to float("inf")
        :param bool plotLogY: to plot logarithmically, defaults to False
        """

        def logBase10Helper(val, takeLog):
            if takeLog:
                return float(np.log10(val))
            else:
                return val

        pointsToPlot = []
        for filename in list(self.PDRsRunsDict.keys()):
            colValues = self.PDRsRunsDict[filename].getCol(columnHeader)
            if typeToPlot.lower() == "max":
                yVal = max(colValues)
            elif typeToPlot.lower() in ["avg", "average"]:
                yVal = sum(colValues) / len(colValues)
            else:
                raise Exception("Unimplemented type")
            
            yVal = max(yVal, yMinVal)
            yVal = min(yVal, yMaxVal)
            pointsToPlot.append(logBase10Helper(yVal, plotLogY))
        
        RegularPlotMaker = RegularPlotService([i for i in range(1, len(self.PDRsRunsDict.keys()) + 1)], [pointsToPlot])
        RegularPlotMaker.makeRegularPlot(axObj, "File", yAxisLabel, titleLabel)