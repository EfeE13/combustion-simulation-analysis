import os
import numpy as np
from scipy.special import erfcinv
from src.utils.plotting_utils.DataPlot import DataPlot

def getInputFiles(inputFolder:str) -> list[str]:
    """Gets all files in inputFolder; unwanted files will be weeded out in another function

    :param str inputFolder: input folder
    :return list[str]: list of files in inputFolder
    """
    inputFileLocs = []
    for root, _, files in os.walk(inputFolder):
        for file in files:
            inputFileLocs.append(os.path.join(root, file))
    return inputFileLocs

def findChiZAtZSt(ZVals:list[float], chiZAtZ:list[float], ZStVals:list[float]) -> list[float]:
    """Finds dissipation rate at stoichiometric conditions using erfc forula

    :param list[float] ZVals: the Z values for which the dissipation rates are known
    :param list[float] chiZAtZ: the dissipation rates that are known
    :param list[float] ZStVals: the Z values corresponding to stoichiometric conditions
    :return list[float]: the dissipation rate at stoichiometric conditions
    """
    ZVals = np.array(ZVals); chiZAtZ = np.array(chiZAtZ); ZStVals = np.array(ZStVals)
    assert ZVals.shape == chiZAtZ.shape
    assert ZVals.shape == ZStVals.shape

    chiZAtZSt = np.zeros(ZVals.shape)
    for i in range(ZVals.shape[0]):
        for j in range(ZVals.shape[1]):
            chiZAtZSt[i][j] = chiZAtZ[i][j] / np.exp(-2 * (erfcinv(2 * ZVals[i][j])) ** 2) * np.exp(-2 * (erfcinv(2 * ZStVals[i][j])) ** 2)
        
    return chiZAtZSt

def makeAndSaveFigure(timeStep, xVar:str, yVar:str, xAxisLabel:str, yAxisLabel:str, colorbarLabel:str, figureFileLoc:str, zVar:str = None, zList:list[float] = None,
                      titleLabel:str = None, zMinVal:float = float("-inf"), zMaxVal:float = float("inf"), plotLogX:bool = False, 
                      plotLogY:bool = False, plotLogZ:bool = False, increasedDotSize:float = None, diffShapes = False, diffShapesLegend:list[str] = None) -> None:
    """Making and saving figure

    :param _type_ timeStep: TimeStep object
    :param str xVar: variable for x-axis of color plot
    :param str yVar: variable for y-axis of color plot
    :param str xAxisLabel: x-axis label
    :param str yAxisLabel: y-axis label
    :param str colorbarLabel: colorbar label
    :param str figureFileLoc: where to save the files
    :param str zVar: variable for coloring the dots of the color plot (or one can supply zList), defaults to None
    :param list[float] zList: for coloring the dots of the color plot (or one can supply zVar), defaults to None
    :param str titleLabel: title label, defaults to None
    :param float zMinVal: lower cap for coloring, defaults to float("-inf")
    :param float zMaxVal: upper cap for coloring, defaults to float("inf")
    :param bool plotLogX: to plot x-axis logarithmically, defaults to False
    :param bool plotLogY: to plot xy-axis logarithmically, defaults to False
    :param bool plotLogZ: to plot color logarithmically, defaults to False
    :param float increasedDotSize: size to increase dots, defaults to None
    :param bool diffShapes: whether or not to have different dots be different shapes, defaults to False
    :param list[str] diffShapesLegend: for making the legend of the different shapes, defaults to None
    """
    if titleLabel == None:
        titleLabel = colorbarLabel

    dataVisObj = DataPlot(1, 1)
    timeStep.makeHeatMap(dataVisObj.getFreeSubplot(), xVar, yVar, xAxisLabel, yAxisLabel, titleLabel, colorbarLabel, zVar = zVar, zList = zList, zMinVal = zMinVal, zMaxVal = zMaxVal, plotLogX = plotLogX, plotLogY = plotLogY, plotLogZ = plotLogZ, increasedDotSize = increasedDotSize, diffShapes = diffShapes, diffShapesLegend = diffShapesLegend)
    dataVisObj.saveFigure(figureFileLoc[:figureFileLoc.rindex(".")] + "Title" + figureFileLoc[figureFileLoc.rindex("."):], 500)

    dataVisObj = DataPlot(1, 1)
    timeStep.makeHeatMap(dataVisObj.getFreeSubplot(), xVar, yVar, xAxisLabel, yAxisLabel, "", colorbarLabel, zVar = zVar, zList = zList, zMinVal = zMinVal, zMaxVal = zMaxVal, plotLogX = plotLogX, plotLogY = plotLogY, plotLogZ = plotLogZ, increasedDotSize = increasedDotSize, diffShapes = diffShapes, diffShapesLegend = diffShapesLegend)
    dataVisObj.saveFigure(figureFileLoc[:figureFileLoc.rindex(".")] + "NoTitle" + figureFileLoc[figureFileLoc.rindex("."):], 500)
