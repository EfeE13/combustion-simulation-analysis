import os
import struct
import random
import numpy as np
import subprocess
import copy
import math
from typing import Dict

from src.pdrs_analysis.PDRsRunnable import PDRsRunnable
from src.pdrs_analysis.PDRsRun import PDRsRun
from src.utils.plotting_utils.ColorMapService import ColorMapService
from src.nga_analysis import utils as NGAUtils
from src.utils import utils
from src.nga_analysis import constants

class TimeStep:
    """
    TimeStep represents the data of some (not necessarily proper) subset of the cells at one timestep

    Attributes
    data: A dictionary like {"ZMIX":[0.3849, 0.5943, 0.138294], "FMIX":[0.4389754, 0.2839, 0.5403], ...}. Here, the ith element of every list corresponds to the ith cell.
    quantities: A list of strings that are the keys of self.data
    timestepStr: A str, something like "000124" that gives the file extension of the time step to use
    numCells: An int, the length of each of the self.data lists
    """

    def __init__(self, inputFileLocs:list[str], quantities:list[str], timestepStr:str, binary:bool = None, data:Dict[str, list[float]] = None) -> None:
        """
        :param list[str] inputFileLocs: each string in list represents a file location
        :param list[str] quantities: each string in list represents a key that you want in self.data keys.
        :param str timestepStr: represents the file extension of the inputFileLocs you want to keep (e.g., "000124")
        :param bool binary: represents if the input files are in binary format (e.g., ensight) (True) or not (False). If False, a particular text format is needed, defaults to None
        :param Dict[str, list[float]] data: if you already know the data, a very short constructor, defaults to None
        """
        if data != None:
            self.data = copy.deepcopy(data)
            self.timestepStr = timestepStr
        else:
            assert binary != None
            self.timestepStr = timestepStr
            self.data = {quantity:[] for quantity in quantities}
            
            for inputFileLoc in inputFileLocs:
                if ("." in inputFileLoc) and (inputFileLoc[inputFileLoc.rindex(".") + 1:] == timestepStr):
                    if "/" in inputFileLoc:
                        inputFileName = inputFileLoc[inputFileLoc.rindex("/") + 1:]
                    else:
                        inputFileName = inputFileLoc
                    
                    for quantity in quantities:
                        if quantity == inputFileName[:inputFileName.rindex(".")]:
                            if binary:
                                # Reading the binary data files in ensight folder
                                fileSize = os.path.getsize(inputFileLoc)
                                numFloats = fileSize // 4
                                inputFile = open(inputFileLoc, "rb")
                                dataVals = inputFile.read()
                                inputFile.close()
                                assert len(dataVals) == fileSize
                                numbers = list(struct.unpack(f'{numFloats}f', dataVals))
                                self.data[quantity] += numbers
                            else:
                                # Reading text files stored by TimeStep's storeData function on a previous run of this Python file
                                # Advantage: If you do a lot of processing on the raw ensight data using TimeStep's functions, you 
                                # can save the processed data to not have to re-process in later runs of the main Python code
                                inputFile = open(inputFileLoc, "r")
                                inputFileContents = inputFile.read() + "\n"
                                inputFile.close()
                                
                                lines = inputFileContents.split("\n---\n") # This \n---\n separates the data for each cell
                                for line in lines:
                                    if line.strip() != "":
                                        if line.strip() in ["True", "False"]:
                                            valueToAdd = bool(line.strip())
                                        else:
                                            try:
                                                if "BestQ2DFModel" in inputFileLoc:
                                                    assert False
                                                valueToAdd = float(line.strip())
                                            except:
                                                valueToAdd = line.strip()
                                        self.data[quantity].append(valueToAdd)

        self.__assertNumCellsCorrect()

    def storeData(self, storageLocation:str) -> None:
        """If you do a lot of processing on the raw ensight data using TimeStep's functions, you can save the processed data to not have to re-process in later runs of this Python file

        :param str storageLocation: the location to store all of the data of the TimeStep object (all previous files in storageLocation are removed)
        """
        storageLocation = utils.fixDirFormat(storageLocation)
        utils.makeDir(storageLocation)
        # Clears the data you had in storageLocation previously:
        subprocess.run("rm " + storageLocation + "/*", shell = True)
        for quantity in self.getQuantities():
            outputFileContents = ""
            for i in range(self.getNumCells()):
                outputFileContents += str(self.getData(quantity, i)) + "\n---\n"
            
            outputFileLoc = storageLocation + "/" + quantity + "." + self.timestepStr
            utils.run("touch " + outputFileLoc)
            outputFile = open(outputFileLoc, "w")
            outputFile.write(outputFileContents.strip())
            outputFile.close()

    def addQuantity(self, quantity:str) -> None:
        """A large part of the TimeStep object's functionality. Is used to compute new quantities from the existing ones

        :param str quantity: the quantity to be added (the formula for which must be given inside this function)
        """

        if quantity in self.getQuantities():
            raise ValueError(quantity + " already in the data.")

        # Some quantities commonly used in the below formulas that may be useful to have on hand:
        ZSTAR = np.array(self.getDataList("ZSTAR"))
        ZMIX = np.array(self.getDataList("ZMIX"))
        C_ZST = np.array(self.getDataList("C_ZST"))
        C_ZZST = np.array(self.getDataList("C_ZZST"))
        CHI = np.array(self.getDataList("CHI"))
        if "Q2DF" in quantity:
            modelNum = quantity[quantity.index("Q2DF") + 4]   # A str that is one character long
            if "Q2DF" + modelNum + "Eta" in self.getQuantities():
                ETA = np.array(self.getDataList("Q2DF" + modelNum + "Eta"))

        if quantity in ["FMIX"]:
            self.data[quantity] = list(ZSTAR / ZMIX)
        elif quantity in ["Q2DF1ChiEta"]:
            self.data[quantity] = list((CHI-2*(1-ZMIX)/(1-ZSTAR)*C_ZZST+((1-ZMIX)/(1-ZSTAR))**2*C_ZST)/((1-ZSTAR)**2))
        elif quantity in ["Q2DF2ChiEta"]:
            self.data[quantity] = list((C_ZST-2*ZSTAR/ZMIX*C_ZZST+(ZSTAR/ZMIX)**2*CHI)/(ZMIX**2))
        elif quantity in ["Q2DF3ChiEta"]:
            self.data[quantity] = list(C_ZST-2*C_ZZST+CHI)
        elif quantity in ["Q2DF1ChiXi"]:
            self.data[quantity] = list(C_ZST)
        elif quantity in ["Q2DF2ChiXi"]:
            self.data[quantity] = list(CHI)
        elif quantity in ["Q2DF3ChiXi"]:
            self.data[quantity] = list((((1-ZMIX)/(1+ZSTAR-ZMIX))**2*C_ZST+2*ZSTAR*(1-ZMIX)/((1+ZSTAR-ZMIX)**2)*C_ZZST+(ZSTAR/(1+ZSTAR-ZMIX))**2*CHI)/((1+ZSTAR-ZMIX)**2))
        elif quantity in ["Q2DF1Ratio", "Q2DF2Ratio", "Q2DF3Ratio"]:
            # The ratio, chi_eta/chi_xi, that we want to minimize when choosing the model to use
            self.data[quantity] = list(np.array(self.getDataList("Q2DF" + str(modelNum) + "ChiEta")) / np.array(self.getDataList("Q2DF" + str(modelNum) + "ChiXi")))
        elif quantity in ["Q2DF1RatioNorm", "Q2DF2RatioNorm", "Q2DF3RatioNorm"]:
            self.data[quantity] = list(np.array(self.getDataList("Q2DF" + str(modelNum) + "Ratio")) / (np.array(self.getDataList("Q2DF" + str(modelNum) + "Ratio")) + 1))
        elif quantity in ["Ratio12Norm", "Ratio23Norm", "Ratio31Norm"]:
            firstModel = quantity[quantity.index("Norm") - 2]
            secondModel = quantity[quantity.index("Norm") - 1]
            comparisonRatio = np.array(self.getDataList("Q2DF" + firstModel + "Ratio")) / np.array(self.getDataList("Q2DF" + secondModel + "Ratio"))
            self.data[quantity] = list(comparisonRatio / (comparisonRatio + 1))
        elif quantity in ["BestQ2DFModel", "BestQ2DFModel12", "BestQ2DFModel13", "BestQ2DFModel23"]:
            # The data here is always stored under the tag "BestQ2DFModel" so be careful if you passed in something like "BestQ2DFModel12" as the argument to this function
            # BestQ2DFModel12 gives the best model of models 1 and 2, excluding the performance of model 3
            if quantity == "BestQ2DFModel":
                excludedNums = []
            elif quantity == "BestQ2DFModel12":
                excludedNum = [3]
            elif quantity == "BestQ2DFModel13":
                excludedNums = [2]
            elif quantity == "BestQ2DFModel23":
                excludedNums = [1]
            self.data["BestQ2DFModel"] = []
            for i in range(self.getNumCells()):
                ratioValues = [self.getData("Q2DF1Ratio", i), self.getData("Q2DF2Ratio", i), self.getData("Q2DF3Ratio", i)]
                for excludedNum in excludedNums:
                    ratioValues[excludedNum - 1] = float("inf")
                for ratioValue in ratioValues:
                    assert ratioValue > 0
                self.data["BestQ2DFModel"].append(str(ratioValues.index(min(ratioValues)) + 1))
            print()
            self.__assertBestQ2DFModelCorrect()
            
        elif quantity in ["RatioCompetitor"]:
            # How many times better the best model is than the second best, as expressed by the ratio of chi_eta/chi_xi for the two models
            self.data[quantity] = []
            for i in range(self.getNumCells()):
                ratioValues = [self.getData("Q2DF1Ratio", i), self.getData("Q2DF2Ratio", i), self.getData("Q2DF3Ratio", i)]
                for ratioValue in ratioValues:
                    assert ratioValue > 0
                minRatio = min(ratioValues)
                ratioValues.remove(minRatio)
                assert len(ratioValues) == 2
                secondLowestRatio = min(ratioValues)
                self.data[quantity].append(secondLowestRatio / minRatio)
        elif quantity in ["Q2DF1Eta"]:
            self.data[quantity] = list((ZMIX - ZSTAR) / (1 - ZSTAR))
        elif quantity in ["Q2DF2Eta"]:
            self.data[quantity] = list((ZMIX - ZSTAR) / ZMIX)
        elif quantity in ["Q2DF3Eta"]:
            self.data[quantity] = list(ZMIX - ZSTAR)
        elif quantity in ["Q2DF1Xi"]:
            self.data[quantity] = list(ZSTAR)
        elif quantity in ["Q2DF2Xi"]:
            self.data[quantity] = list(ZMIX)
        elif quantity in ["Q2DF3Xi"]:
            self.data[quantity] = list(ZSTAR / (1 + ZSTAR - ZMIX))
        elif quantity in ["BestMVal"]:
            self.data[quantity] = list(((C_ZST - CHI) - np.sqrt((C_ZST-CHI)*(C_ZST-CHI) + 4*C_ZZST*C_ZZST)) / (2*C_ZZST))
            #self.data[quantity] = list(-C_ZZST / C_ZST)
        elif quantity in ["BestXVal"]:
            bestMVal = np.array(self.getDataList("BestMVal"))
            self.data[quantity] = list((1 - (ZSTAR - (1-ZMIX)/bestMVal)) / (1 - (1 - ZMIX - ZSTAR * bestMVal)))
        elif quantity in ["BestXValNormalized"]:
            bestXVal = np.array(self.getDataList("BestXVal"))
            self.data[quantity] = list(bestXVal / (bestXVal + 1))
        else:
            # The place for defining new quantities you want to add is here, in this if-elif-elif-... section, open/closed principle :(
            raise ValueError("Variable " + quantity + " not implemented")

        self.__assertNumCellsCorrect()

    def removeCellsWithValue(self, comparisonQuantity:str, comparisonValue, comparisonType:str) -> None:
        """Calls __removeCellsAtIndices helper function to remove cells based on certain values of comparisonQuantity. Example: Calling removeCellsWithValue("ZMIX", 0.3, "gt") would remove cells with ZMIX > 0.3

        :param str comparisonQuantity: the variable that the reference value refers tos
        :param _type_ comparisonValue: the reference value
        :param str comparisonType: how to compare cells to this reference value
        """
        
        cellsToRemove = []
        for i in range(self.getNumCells()):
            assert comparisonType in ["eqexact", "eq", "lt", "gt", "lte", "gte", "ne", "inrangeinclusive", "inrangeexclusive", "validm", "cosbetweenneg1andpos1"]
            removeEqExact = (comparisonType == "eqexact") and (self.data[comparisonQuantity][i] == comparisonValue)
            removeEq = (comparisonType == "eq") and (abs(self.data[comparisonQuantity][i] - comparisonValue) < 0.00001) # Almost equal
            removeLt = (comparisonType == "lt") and (self.data[comparisonQuantity][i] < comparisonValue)
            removeGt = (comparisonType == "gt") and (self.data[comparisonQuantity][i] > comparisonValue)
            removeLte = (comparisonType == "lte") and (self.data[comparisonQuantity][i] <= comparisonValue)
            removeGte = (comparisonType == "gte") and (self.data[comparisonQuantity][i] >= comparisonValue)
            removeNe = (comparisonType == "ne") and (self.data[comparisonQuantity][i] != comparisonValue)
            removeInRangeInc = (comparisonType == "inrangeinclusive") and (self.data[comparisonQuantity][i] >= comparisonValue[0]) and (self.data[comparisonQuantity][i] <= comparisonValue[1])
            removeInRangeExc = (comparisonType == "inrangeexclusive") and (self.data[comparisonQuantity][i] > comparisonValue[0]) and (self.data[comparisonQuantity][i] < comparisonValue[1])
            
            if comparisonType == "validm":
                z1Val = self.data["ZSTAR"][i]
                z2Val = 1.0 - self.data["ZMIX"][i]
                minM = -(1.0 - z2Val) / z1Val
                maxM = -z2Val / (1.0 - z1Val)
                mVal = self.data["BestMVal"][i]
                removeM = (minM >= mVal) or (mVal >= maxM)
            else:
                removeM = False

            if comparisonType == "cosbetweenneg1andpos1":
                chiVal = self.data["CHI"][i]
                chiZstVal = self.data["C_ZST"][i]
                chiZZstVal = self.data["C_ZZST"][i]
                removeCosInRange = abs(chiZZstVal) >= math.sqrt(chiVal * chiZstVal)
            else:
                removeCosInRange = False

            if (removeEqExact or removeEq or removeLt or removeGt or removeLte or removeGte or removeNe or removeInRangeInc or removeInRangeExc or removeM or removeCosInRange):
                cellsToRemove.append(i)
        
        self.__removeCellsAtIndices(cellsToRemove)
    
    def reduceNumCells(self, **kwargs) -> None:
        """Calls __removeCellsAtIndices helper function to reduced the number of cells
        Examples:
        reduceNumCells(num = 200): Reduces cells to a randomly chosen subset of 200 cells
        reduceNumCells(factor = 3): Makes the number of cells a third of what it was (randomly)
        reduceNumCells(special1 = [["1", "2", "3"], [30, 30, 30]]) which is list[list[str], list[int]]: Takes 30 points from the subset of 
            cells with model 1 best, 30 from the subset with model 2 best, etc. But less than 30 if there are not 30 available.
        """
        key = list(kwargs.keys())[0]
        val = list(kwargs.values())[0]
        oldNumCells = self.getNumCells()

        if key in ["num", "factor"]:
            if key == "num":
                newNumCells = int(val)
            elif key == "factor":
                newNumCells = oldNumCells // val

            self.__removeCellsAtIndices(random.sample(range(self.getNumCells()), oldNumCells - newNumCells))
        elif key == "special1":
            assert len(val[0]) == len(val[1])
            cellsToRemove = []
            print()
            for i in range(self.getNumCells()):
                assert self.getData("BestQ2DFModel", i) in ["1", "2", "3"]

            for modelType in range(len(val[0])):
                model = val[0][modelType]
                numberToKeep = int(val[1][modelType])
                indicesWithModel = [i for i in range(self.getNumCells()) if self.getData("BestQ2DFModel", i) == model]
                addToCellsToRemove = random.sample(indicesWithModel, max(len(indicesWithModel) - numberToKeep, 0))
                cellsToRemove += addToCellsToRemove
                print(model, "|", len(indicesWithModel), "-", numberToKeep, "|", len(indicesWithModel), "getting", len(addToCellsToRemove), "removed")

            self.__removeCellsAtIndices(cellsToRemove)
        else:
            raise ValueError()

    def removeAllButNearestCells(self, coordNames:list[str], coords:list[list[float]]) -> None:
        """This removes all but the nearest cell to (xCoordName = coord[0], yCoordName = coord[1], zCoordName = coord[2], ...) for each coord in coords
        Example: Calling removeAllButNearestCells(["ZMIX", "FMIX"], [[0.2, 0.23], [0.02, 0.1], [0.02, 0.3]]) removes all but three cells: those closest to (ZMIX, FMIX) = (0.2, 0.23), (0.02, 0.01) and (0.02, 0.3)
        It gets a little messy when a point is close to two coords, check implementation for details.

        :param list[str] coordNames: List of strings that represent what each of the elements in a sublist of coords represent
        :param list[list[float]] coords: List of lists representing the coordinates to keep points nearby to
        """

        cellsToKeep = []
        for coord in coords:
            squaredDistances = []
            assert len(coord) == len(coordNames)
            for i in range(self.getNumCells()):
                squaredDistance = 0
                for j in range(len(coordNames)):
                    squaredDistance += (self.getData(coordNames[j], i) - coord[j]) ** 2
                squaredDistances.append(squaredDistance)
            while squaredDistances.index(min(squaredDistances)) in cellsToKeep:
                squaredDistances[squaredDistances.index(min(squaredDistances))] = float("inf")
            cellsToKeep.append(squaredDistances.index(min(squaredDistances)))
        
        cellsToRemove = [i for i in range(self.getNumCells()) if i not in cellsToKeep]
        self.__removeCellsAtIndices(cellsToRemove)
    
    def getQuantities(self) -> list[str]:
        return list(self.data.keys())
    
    def getNumCells(self) -> int:
        return len(self.data[self.getQuantities()[0]])

    def getDataList(self, quantity:str) -> list[float]:
        self.__assertQuantityPresent(quantity)
        return self.data[quantity]
    
    def getData(self, quantity:str, cellNum:int) -> float:
        self.__assertQuantityPresent(quantity)
        return self.getDataList(quantity)[cellNum]

    def getMin(self, quantity:str) -> float:
        self.__assertQuantityPresent(quantity)
        return min(self.getDataList(quantity))
    
    def getMax(self, quantity:str) -> float:
        self.__assertQuantityPresent(quantity)
        return max(self.getDataList(quantity))

    # log x and y functionality not implemented, so I removed those parameters from the function call
    def makeHeatMap(self, axObj, xVar:str, yVar:str, xAxisLabel:str, yAxisLabel:str, titleLabel:str, colorbarLabel:str, 
                    zVar:str = None, zList:list[float] = None, zMinVal:float = float("-inf"), zMaxVal:float = float("inf"), plotLogX:bool = False, 
                    plotLogY:bool = False, plotLogZ:bool = False, dotSize:float = None, diffShapes = False, diffShapesLegend:list[str] = None) -> None:
        """Makes colored scatter plot based on the parameters below

        :param _type_ axObj: axis object to be plotted on.
        :param str xVar: variable for x-axis of color plot
        :param str yVar: variable for y-axis of color plot
        :param str xAxisLabel: x-axis label
        :param str yAxisLabel: y-axis label
        :param str titleLabel: title label
        :param str colorbarLabel: colorbar label
        :param str zVar: variable for coloring the dots of the color plot (or one can supply zList), defaults to None
        :param list[float] zList: for coloring the dots of the color plot (or one can supply zVar), defaults to None
        :param float zMinVal: lower cap for coloring, defaults to float("-inf")
        :param float zMaxVal: upper cap for coloring, defaults to float("inf")
        :param bool plotLogX: to plot x-axis logarithmically, defaults to False
        :param bool plotLogY: to plot y-axis logarithmically, defaults to False
        :param bool plotLogZ: to plot color logarithmically, defaults to False
        :param float dotSize: size of dots, defaults to None
        :param bool diffShapes: whether or not to have different dots be different shapes, defaults to False
        :param list[str] diffShapesLegend: for making the legend of the different shapes, defaults to None
        """
        
        def zBoundToVBound(num:float) -> float:
            if num in [float("-inf"), float("inf")]:
                return None
            else:
                return num

        print("Making heat map " + titleLabel)
        zVarGiven = zVar is not None
        zListGiven = zList is not None
        assert zVarGiven != zListGiven

        xList = np.clip(np.array(self.getDataList(xVar)).astype(float), -constants.CLIPPING_VALUE_COLOR_MAP, constants.CLIPPING_VALUE_COLOR_MAP)
        yList = np.clip(np.array(self.getDataList(yVar)).astype(float), -constants.CLIPPING_VALUE_COLOR_MAP, constants.CLIPPING_VALUE_COLOR_MAP)
        if not zListGiven:
            zList = np.array(self.getDataList(zVar)).astype(float)
            #zList = np.clip(np.array(self.getDataList(zVar)).astype(float), zMinVal, zMaxVal)
        zList = np.clip(zList, zMinVal, zMaxVal)
        
        if plotLogX:
            xList = np.log10(xList)
        if plotLogY:
            xList = np.log10(yList)
        if plotLogZ:
            xList = np.log10(zList)

        if diffShapes:
            model1Best = [i for i in range(self.getNumCells()) if self.getData("BestQ2DFModel", i) == "1"]
            model2Best = [i for i in range(self.getNumCells()) if self.getData("BestQ2DFModel", i) == "2"]
            model3Best = [i for i in range(self.getNumCells()) if self.getData("BestQ2DFModel", i) == "3"]

            assert len(model1Best) + len(model2Best) + len(model3Best) == self.getNumCells()
            diffShapesIndices = [model1Best, model2Best, model3Best]
        else:
            diffShapesIndices = None

        myColorMapService = ColorMapService([[xList[i], yList[i], zList[i]] for  i in range(len(xList))])
        myColorMapService.makeColorPlotNonContinuous(axObj, xAxisLabel, yAxisLabel, titleLabel, colorbarLabel, dotSize = dotSize, diffShapes = diffShapes, diffShapesIndices = diffShapesIndices, diffShapesLegend = diffShapesLegend, vMin = zBoundToVBound(zMinVal), vMax = zBoundToVBound(zMaxVal))
    
    def prepareForQ2DFModels(self, identifier:str, Z1Comp:Dict[str, float], Z2Comp:Dict[str, float], Z3Comp:Dict[str, float], Z1Temp:float, Z2Temp:float, Z3Temp:float, molWeights:Dict[str, float], nuVals:Dict[str, float]):
        """Makes key lists of inputs for the running of PDRs

        :param str identifier: ome nickname for the runs
        :param Dict[str, float] Z1Comp: composition of stream 1
        :param Dict[str, float] Z2Comp: composition of stream 2
        :param Dict[str, float] Z3Comp: composition of stream 3
        :param float Z1Temp: temperature of stream 1 (K)
        :param float Z2Temp: temperature of stream 2 (K)
        :param float Z3Temp: temperature of stream 3 (K)
        :param Dict[str, float] molWeights: molecular weights
        :param Dict[str, float] nuVals: nu values (-1 for O2, i.e., negative for oxidizer)
        """
        def toStr(num:float) -> str:
            return str(round(num, 5)).replace(".", "~")
        
        allChemCompounds = list(set(set(Z1Comp.keys()) | set(Z2Comp.keys()) | set(Z3Comp.keys())))
        for chemCompound in allChemCompounds:
            assert chemCompound in molWeights.keys()
            assert chemCompound in nuVals.keys()
            for dictionary in [Z1Comp, Z2Comp, Z3Comp]:
                if chemCompound not in dictionary.keys():
                    dictionary[chemCompound] = 0
        
        
        for dictionary in Z1Comp, Z2Comp, Z3Comp:
            assert abs(sum(list(dictionary.values())) - 1) == 0

        Z1 = np.array(self.getDataList("ZSTAR"))
        Z2 = 1 - np.array(self.getDataList("ZMIX"))
        Z3 = np.array(self.getDataList("ZMIX")) - np.array(self.getDataList("ZSTAR"))
        
        shape = Z1.shape
        assert len(shape) == 1
        assert shape[0] == self.getNumCells()
        assert Z2.shape == shape
        assert Z3.shape == shape

        xi1 = np.array(self.getDataList("Q2DF1Xi"))
        eta1 = np.array(self.getDataList("Q2DF1Eta"))
        xi2 = np.array(self.getDataList("Q2DF2Xi"))
        eta2 = np.array(self.getDataList("Q2DF2Eta"))
        xi3 = np.array(self.getDataList("Q2DF3Xi"))
        eta3 = np.array(self.getDataList("Q2DF3Eta"))

        oxidZ1 = [np.zeros(shape), np.zeros(shape), np.zeros(shape)]
        oxidZ2 = [1 - eta1, np.ones(shape), 1 - eta3]
        oxidZ3 = [eta1, np.zeros(shape), eta3]
        fuelZ1 = [np.ones(shape), 1 - eta2, 1 - eta3]
        fuelZ2 = [np.zeros(shape), np.zeros(shape), np.zeros(shape)]
        fuelZ3 = [np.zeros(shape), eta2, eta3]

        oxidTemp = [oxidZ1[i] * Z1Temp + oxidZ2[i] * Z2Temp + oxidZ3[i] * Z3Temp for i in [0, 1, 2]]
        fuelTemp = [fuelZ1[i] * Z1Temp + fuelZ2[i] * Z2Temp + fuelZ3[i] * Z3Temp for i in [0, 1, 2]]
        speciesToOxidBC = {}
        speciesToFuelBC = {}
        for chemCompound in allChemCompounds:
            speciesToOxidBC[chemCompound] = np.array([oxidZ1[i] * Z1Comp[chemCompound] + oxidZ2[i] * Z2Comp[chemCompound] + oxidZ3[i] * Z3Comp[chemCompound] for i in [0, 1, 2]])
            speciesToFuelBC[chemCompound] = np.array([fuelZ1[i] * Z1Comp[chemCompound] + fuelZ2[i] * Z2Comp[chemCompound] + fuelZ3[i] * Z3Comp[chemCompound] for i in [0, 1, 2]])

        fuelSpeciesStr = [[0 for j in range(shape[0])] for i in range(3)]
        oxidBC = [[0 for j in range(shape[0])] for i in range(3)]
        fuelBC = [[0 for j in range(shape[0])] for i in range(3)]
        for i in range(3):
            for j in range(shape[0]):
                oxidBCStr = ""
                fuelBCStr = ""
                fuelSpeciesList = []
                for chemCompound in allChemCompounds:
                    if round(speciesToOxidBC[chemCompound][i][j], 12) > 0:
                        oxidBCStr += chemCompound + " " + str(round(speciesToOxidBC[chemCompound][i][j], 12)) + "\n    "
                        if (chemCompound not in fuelSpeciesList) and (chemCompound not in constants.LIST_NONFUELS):
                            fuelSpeciesList.append(chemCompound)
                    if round(speciesToFuelBC[chemCompound][i][j], 12) > 0:
                        fuelBCStr += chemCompound + " " + str(round(speciesToFuelBC[chemCompound][i][j], 12)) + "\n    "
                        if (chemCompound not in fuelSpeciesList) and (chemCompound not in constants.LIST_NONFUELS):
                            fuelSpeciesList.append(chemCompound)
                oxidBC[i][j] = oxidBCStr.strip()
                fuelBC[i][j] = fuelBCStr.strip()
                fuelSpeciesStr[i][j] = " ".join(fuelSpeciesList)

        ZVals = [xi1, xi2, xi3]
        chiZAtZ = [np.array(self.getDataList("Q2DF1ChiXi")), np.array(self.getDataList("Q2DF2ChiXi")), np.array(self.getDataList("Q2DF3ChiXi"))]
        oxidKappa = np.zeros(np.array(ZVals).shape)
        fuelKappa = np.zeros(np.array(ZVals).shape)
        for chemCompound in allChemCompounds:
            oxidKappa -= speciesToOxidBC[chemCompound] / molWeights[chemCompound] * nuVals[chemCompound] * molWeights["O2"]
            fuelKappa -= speciesToFuelBC[chemCompound] / molWeights[chemCompound] * nuVals[chemCompound] * molWeights["O2"]
        ZStVals = oxidKappa / (oxidKappa - fuelKappa)
        chiZAtZSt = NGAUtils.findChiZAtZSt(ZVals, chiZAtZ, ZStVals)

        validZSt = np.logical_and(np.greater_equal(ZStVals, np.zeros(ZStVals.shape)), np.less_equal(ZStVals, np.ones(ZStVals.shape)))
        #validZSt = np.append(validZSt, np.array([[validZSt[int(self.getData("BestQ2DFModel", i)) - 1][i] for i in range(self.getNumCells())]]), axis = 0)

        fileNameZeta = [[0 for j in range(shape[0])] for i in range(3)]
        fileNameQ2DF = [0 for j in range(shape[0])]
        fileNameCont = [0 for j in range(shape[0])]
        for j in range(shape[0]):
            descriptiveText = "ZMIX" + toStr(self.getData("ZMIX", j)) + "FMIX" + toStr(self.getData("FMIX", j))
            for i in range(3):
                fileNameZeta[i][j] = identifier + "ZetaQ2DF" + str(i + 1) + descriptiveText
            fileNameQ2DF[j] = identifier + "ManiQ2DF" + descriptiveText
            fileNameCont[j] = identifier + "CONTManiQ2DF" + descriptiveText
        
        # Push all the data into self.data in preparation for PDRs running
        for m in range(3):
            model = str(m + 1)
            assert model in ["1", "2", "3"]
            self.data[identifier + "ZetaQ2DF" + model + "ZValue"] = list(ZVals[m])
            self.data[identifier + "ZetaQ2DF" + model + "OxidTemp"] = list(oxidTemp[m])
            self.data[identifier + "ZetaQ2DF" + model + "FuelTemp"] = list(fuelTemp[m])
            self.data[identifier + "ZetaQ2DF" + model + "OxidBC"] = list(oxidBC[m])
            self.data[identifier + "ZetaQ2DF" + model + "FuelBC"] = list(fuelBC[m])
            self.data[identifier + "ZetaQ2DF" + model + "FuelSpecies"] = list(fuelSpeciesStr[m])
            self.data[identifier + "ZetaQ2DF" + model + "ScalarDissRate"] = list(chiZAtZSt[m])
            self.data[identifier + "ZetaQ2DF" + model + "FileName"] = list(fileNameZeta[m])
            self.data[identifier + "ZetaQ2DF" + model + "Valid"] = list(validZSt[m])
        
        self.data[identifier + "ManiQ2DFFileName"] = fileNameQ2DF
        self.data[identifier + "ManiQ2DFZValue"] = []
        self.data[identifier + "ManiQ2DFValid"] = []
        # TODO: Continuous filename and z value
        self.data[identifier + "CONTManiQ2DFFileName"] = fileNameCont
        for i in range(self.getNumCells()):
            self.data[identifier + "ManiQ2DFZValue"].append(float(ZVals[int(self.getData("BestQ2DFModel", i)) - 1][i]))
            self.data[identifier + "ManiQ2DFValid"].append(validZSt[int(self.getData("BestQ2DFModel", i)) - 1][i])

        self.__removeBadInputs()
        self.__assertNumCellsCorrect()
            
    def runManiZ(self, identifier:str, model:str, pdrsQuantities:list[str], inputFileSpecs:str, inputFolder:str, outputFolder:str):
        """Runs Q2DF using maniZ. Remember to also specify the same output folder in the PDRs input file you are passing in here.

        :param str identifier: some nickname for the runs
        :param str model: "1", "2", "3" or "Best", representing which Q2DF model will be used
        :param list[str] pdrsQuantities: the PDRs outputs you care about
        :param str inputFileSpecs: the input file, but with different keywords like "ZSTAR" replacing where key numbers/inputs would be
        :param str inputFolder: input folder for PDRs
        :param str outputFolder: output folder for PDRs
        """
        self.__removeBadInputs()

        assert model in ["1", "2", "3", "Best"]
        inputFolder = utils.fixDirFormat(inputFolder)
        outputFolder = utils.fixDirFormat(outputFolder)
        
        for pdrsQuantity in pdrsQuantities:
            self.data[identifier + "ZetaQ2DF" + model + "Result" + pdrsQuantity] = []
        for i in range(self.getNumCells()):
            newInputFileSpecs = copy.deepcopy(inputFileSpecs)
            if model == "Best":
                modelStr = self.getData("BestQ2DFModel", i)
            else:
                modelStr = model
            
            for inputEntry in ["OxidTemp", "FuelTemp", "OxidBC", "FuelBC", "FuelSpecies", "ScalarDissRate"]:
                newInputFileSpecs = newInputFileSpecs.replace(inputEntry + "_INPUT", str(self.getData(identifier + "ZetaQ2DF" + modelStr + inputEntry, i)))
            for commentEntry in ["ZMIX", "ZSTAR", "FMIX", "CHI", "C_ZZST", "C_ZST", "BestQ2DFModel"]:
                newInputFileSpecs = newInputFileSpecs.replace(commentEntry + "_INPUT", str(self.getData(commentEntry, i)))

            filename = self.getData(identifier + "ZetaQ2DF" + modelStr + "FileName", i)
            myPDRsRunnable = PDRsRunnable(newInputFileSpecs, [filename + ".multicomponent"], [filename + ".Y"])
            myPDRsRun = PDRsRun(myPDRsRunnable.runPDRs(inputFolder, outputFolder)[0])

            for pdrsQuantity in pdrsQuantities:
                self.data[identifier + "ZetaQ2DF" + model + "Result" + pdrsQuantity].append(myPDRsRun.interpolate("Z", self.getData(identifier + "ZetaQ2DF" + modelStr + "ZValue", i), pdrsQuantity))

    def runManiQ2DF(self, identifier:str, pdrsQuantities:list[str], inputFileSpecs:str, inputFolder:str, outputFolder:str, continuous:bool = False):
        """Runs maniQ2DF (and maniQ2DFCont I think). Remember to also specify the same output folder in the PDRs input file you are passing in here.

        :param str identifier: some nickname for the runs
        :param list[str] pdrsQuantities: the PDRs outputs you care about
        :param str inputFileSpecs: the input file, but with different keywords like "ZSTAR" replacing where key numbers/inputs would be
        :param str inputFolder: input folder for PDRs
        :param str outputFolder: output folder for PDRs
        :param bool continuous: whether Q2DF is continuous or not, defaults to False
        """
        self.__removeBadInputs()

        if continuous:
            typeLabel = "CONTManiQ2DF"
        else:
            typeLabel = "ManiQ2DF"

        inputFolder = utils.fixDirFormat(inputFolder)
        outputFolder = utils.fixDirFormat(outputFolder)
        for pdrsQuantity in pdrsQuantities:
            self.data[identifier + typeLabel + "Result" + pdrsQuantity] = []
        
        for i in range(self.getNumCells()):
            newInputFileSpecs = copy.deepcopy(inputFileSpecs)
            for inputEntry in ["ZMIX", "ZSTAR", "CHI", "C_ZST", "C_ZZST"]:
                newInputFileSpecs = newInputFileSpecs.replace(inputEntry + "_INPUT", str(self.getData(inputEntry, i)))
            for commentEntry in ["BestQ2DFModel"]:
                newInputFileSpecs = newInputFileSpecs.replace(commentEntry + "_INPUT", str(self.getData(commentEntry, i)))

            filename = self.getData(identifier + typeLabel + "FileName", i)
            myPDRsRunnable = PDRsRunnable(newInputFileSpecs, [filename + ".multicomponent"], [filename + ".Y"])
            myPDRsRun = PDRsRun(myPDRsRunnable.runPDRs(inputFolder, outputFolder)[0])

            for pdrsQuantity in pdrsQuantities:
                self.data[identifier + typeLabel + "Result" + pdrsQuantity].append(myPDRsRun.interpolate("Z", self.getData(identifier + typeLabel + "ZValue", i), pdrsQuantity))

    def __removeBadInputs(self) -> None:
        for badInput in constants.BAD_INPUT_FILES_RUNNING_LIST:
            self.removeCellsWithValue("TolAirHepManiQ2DFFileName", badInput, "eqexact")
            for m in ["1", "2", "3"]:
                self.removeCellsWithValue("TolAirHepZetaQ2DF" + m + "FileName", badInput, "eqexact")

    def __removeCellsAtIndices(self, cellsToRemove:list[int]) -> None:
        """Helper function that removes the cells at indices given by cellsToRemove

        :param list[int] cellsToRemove: indices of cells to remove
        """
        if len(cellsToRemove) > 0:
            assert max(cellsToRemove) <= self.getNumCells() - 1
        self.__assertNumCellsCorrect()
        for quantity in self.getQuantities():
            self.data[quantity] = np.delete(np.array(self.data[quantity]), cellsToRemove).tolist()

        self.__assertNumCellsCorrect()
    
    def __assertNumCellsCorrect(self) -> None:
        for quantity in self.getQuantities():
            if len(self.getDataList(quantity)) != self.getNumCells():
                raise Exception(quantity + " only has " + str(len(self.getDataList(quantity))) + " values, but self.getNumCells() = " + str(self.getNumCells()))
    
    def __assertBestQ2DFModelCorrect(self) -> None:
        if "BestQ2DFModel" in self.getQuantities():
            for i in range(self.getNumCells()):
                if self.data["BestQ2DFModel"][i] not in ["1", "2", "3"]:
                    raise ValueError("BestQ2DFModel given the value of " + str(self.data["BestQ2DFModel"][i]))
    
    def __assertQuantityPresent(self, quantity:str) -> None:
        if quantity not in self.getQuantities():
            raise Exception(quantity + " not in self.getQuantities()")