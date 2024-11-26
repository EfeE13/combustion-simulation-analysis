import copy
import numpy as np

from src.nga_analysis import constants
from src.nga_analysis import utils
from src.nga_analysis.TimeStep import TimeStep
from src.utils.plotting_utils.DataPlot import DataPlot

ALL_QUANTITIES = ['ZSTAR', 'ZMIX', 'Temp', 'C_ZST', 'C_ZZST', 'CHI', 'ADDS', 'DIFF', 'P', 'VISC', 'Y_CO2', 'Y_H2', 'Y_O2', 'ZSTA2', 'CHIREF', 'dRHO', 'KSGS', 'RHO', 'Y_CO', 'Y_F', 'Y_H2O', 'Y_OH', 'ZMIX2', 'Q2DF1ChiEta', 'Q2DF2ChiEta', 'Q2DF3ChiEta', 'Q2DF1ChiXi', 'Q2DF2ChiXi', 'Q2DF3ChiXi', 'Q2DF1Eta', 'Q2DF2Eta', 'Q2DF3Eta', 'FMIX', 'Q2DF1Ratio', 'Q2DF2Ratio', 'Q2DF3Ratio', 'Q2DF1Xi', 'Q2DF2Xi', 'Q2DF3Xi', "Q2DF1RatioNorm", "Q2DF2RatioNorm", "Q2DF3RatioNorm", "Ratio12Norm", "Ratio23Norm", "Ratio31Norm", "BestQ2DFModel"]
quants = copy.deepcopy(ALL_QUANTITIES)
quants += ["TolAirHepManiQ2DFFileName", "TolAirHepManiQ2DFZValue", "TolAirHepManiQ2DFValid"]
for model in ["1", "2", "3"]:
    quants += ["TolAirHepZetaQ2DF" + model + "ZValue", "TolAirHepZetaQ2DF" + model + "OxidTemp", "TolAirHepZetaQ2DF" + model + "FuelTemp", "TolAirHepZetaQ2DF" + model + "OxidBC", "TolAirHepZetaQ2DF" + model + "FuelBC", "TolAirHepZetaQ2DF" + model + "FuelSpecies", "TolAirHepZetaQ2DF" + model + "ScalarDissRate", "TolAirHepZetaQ2DF" + model + "FileName", "TolAirHepZetaQ2DF" + model + "Valid"]

def initialDataSetup():
    allQuantities = []
    myTimeStep = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/ensight_input/summer24"), ["ZSTAR", "ZMIX", "Temp", "C_ZST", "C_ZZST", "CHI", 'ADDS', 'DIFF', 'P', 'VISC', 'Y_CO2', 'Y_H2', 'Y_O2', 'ZSTA2', 'CHIREF', 'dRHO', 'KSGS', 'RHO', 'Y_CO', 'Y_F', 'Y_H2O', 'Y_OH', 'ZMIX2'], "000124", True)
    allQuantities += ["ZSTAR", "ZMIX", "Temp", "C_ZST", "C_ZZST", "CHI", 'ADDS', 'DIFF', 'P', 'VISC', 'Y_CO2', 'Y_H2', 'Y_O2', 'ZSTA2', 'CHIREF', 'dRHO', 'KSGS', 'RHO', 'Y_CO', 'Y_F', 'Y_H2O', 'Y_OH', 'ZMIX2']

    myTimeStep.addQuantity("FMIX")
    allQuantities += ["FMIX"]

    myTimeStep.removeCellsWithValue("ZMIX", 0, "lte")
    myTimeStep.removeCellsWithValue("ZMIX", 1, "gte")
    myTimeStep.removeCellsWithValue("FMIX", 0, "lte")
    myTimeStep.removeCellsWithValue("FMIX", 1, "gte")
    myTimeStep.removeCellsWithValue("CHI", 0, "lt")
    myTimeStep.removeCellsWithValue("C_ZST", 0, "lt")
    myTimeStep.removeCellsWithValue("", -100000.0, "cosbetweenneg1andpos1")

    quantitiesToAdd = ["Q2DF1ChiEta", "Q2DF2ChiEta", "Q2DF3ChiEta", "Q2DF1ChiXi", "Q2DF2ChiXi", "Q2DF3ChiXi"]
    quantitiesToAdd += ["Q2DF1Eta", "Q2DF2Eta", "Q2DF3Eta", "Q2DF1Xi", "Q2DF2Xi", "Q2DF3Xi"]
    quantitiesToAdd += ["Q2DF1Ratio", "Q2DF2Ratio", "Q2DF3Ratio"]
    quantitiesToAdd += ["Q2DF1RatioNorm", "Q2DF2RatioNorm", "Q2DF3RatioNorm", "Ratio12Norm", "Ratio23Norm", "Ratio31Norm", "BestQ2DFModel", "RatioCompetitor"]
    for quantity in quantitiesToAdd:
        myTimeStep.addQuantity(quantity)
    allQuantities += quantitiesToAdd

    if True:
        origNumCells = myTimeStep.getNumCells()
        print(0, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF1ChiEta", 0, "lt")
        print(1, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF1ChiXi", 0, "lt")
        print(2, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF2ChiEta", 0, "lt")
        print(3, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF2ChiXi", 0, "lt")
        print(4, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF3ChiEta", 0, "lt")
        print(5, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF3ChiXi", 0, "lt")
        print(6, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF1Xi", 0, "lte")
        print(9, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF1Xi", 1, "gte")
        print(10, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF1Eta", 0, "lte")
        print(11, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF1Eta", 1, "gte")
        print(12, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF2Xi", 0, "lte")
        print(13, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF2Xi", 1, "gte")
        print(14, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF2Eta", 0, "lte")
        print(15, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF2Eta", 1, "gte")
        print(16, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF3Xi", 0, "lte")
        print(17, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF3Xi", 1, "gte")
        print(18, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF3Eta", 0, "lte")
        print(19, myTimeStep.getNumCells())
        myTimeStep.removeCellsWithValue("Q2DF3Eta", 1, "gte")
        print("Comparison:", origNumCells, myTimeStep.getNumCells())
        assert origNumCells == myTimeStep.getNumCells()

    print("v" * 50 + "\n", allQuantities, "\n" + "^" * 50)

    myTimeStep.storeData("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/FullSquareFullNumPoints")

def trimmingData():
    myTimeStepFull = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/FullSquareFullNumPoints"), ALL_QUANTITIES, "000124", False)
    assert myTimeStepFull.getNumCells() == 46135
    myTimeStepFull.removeCellsWithValue("ZMIX", 0.0001, "lt")
    myTimeStepFull.removeCellsWithValue("ZMIX", 0.9999, "gt")
    print("Full Size:", myTimeStepFull.getNumCells())
    myTimeStepFull.storeData("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/ZMIXExtremesTrimmedFullNumPoints")

    myTimeStepFullToReduced = copy.deepcopy(myTimeStepFull)
    myTimeStepFullToReduced.reduceNumCells(num = 1000)
    myTimeStepFullToReduced.storeData("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/ZMIXExtremesTrimmedFewerPoints")

def paperNonPDRsPlotting(dataLoc:str):
    myTimeStep = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/" + dataLoc), ALL_QUANTITIES, "000124", False)
    assert myTimeStep.getNumCells() == 1000

    folderLoc = "/home/efeeroz/Documents/CombustionModelAnalysis/graphs/paperfigs"
    utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "Temperature (K)", folderLoc + "/Temperature.pdf", zVar = "Temp")
    utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "Preferred Model Based on $\psi$", folderLoc + "/BestQ2DFModel.pdf", zMinVal = 1.0, zMaxVal = 4.0, zVar = "BestQ2DFModel")
    for i in [1, 2, 3]:
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "$\psi_" + str(i) + " / (\psi_" + str(i) + " + 1)$", folderLoc + "/Psi" + str(i) + "Norm.pdf", zVar = "Q2DF" + str(i) + "RatioNorm")
    for compareValues in [[1, 2], [2, 3], [3, 1]]:
        i = str(compareValues[0]); j = str(compareValues[1])
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "$\psi_" + i + "/\psi_" + j + "$ Normalized", folderLoc + "/Ratio" + i + j + "Norm.pdf", zVar = "Ratio" + i + j + "Norm")
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "$\psi_" + i + "/\psi_" + j + "$ Normalized", folderLoc + "/Ratio" + i + j + "NormCutoff.pdf", zVar = "Ratio" + i + j + "Norm", zMinVal = 0.0, zMaxVal = 1.0)
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "$\psi_" + i + "/\psi_" + j + "$ Normalized", folderLoc + "/Ratio" + i + j + "NormDiffShapes.pdf", zVar = "Ratio" + i + j + "Norm", zMinVal = 0.0, zMaxVal = 1.0, dotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"])

    myTimeStepForX = copy.deepcopy(myTimeStep)
    myTimeStepForX.addQuantity("BestMVal")
    myTimeStepForX.removeCellsWithValue("", -10000.0, "validm")
    myTimeStepForX.addQuantity("BestXVal")
    myTimeStepForX.addQuantity("BestXValNormalized")
    utils.makeAndSaveFigure(myTimeStepForX, "ZMIX", "FMIX", "$Z$", "$F$", "$x$", folderLoc + "/xValNorm.pdf", titleLabel = "$x$ Value", zMinVal = 0.0, zMaxVal = 1.0, zVar = "BestXValNormalized")

def PDRsPreparation():
    myTimeStepFull = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/ZMIXExtremesTrimmedFullNumPoints"), ALL_QUANTITIES, "000124", False)
    myTimeStepFull.prepareForQ2DFModels("TolAirHep", {"A1CH3":1}, {"O2":0.23292, "N2":0.76708}, {"NXC7H16":1}, 300, 300, 300, {"A1CH3":92.14, "O2":31.999, "N2":28.01, "NXC7H16":100.21}, {"A1CH3":9, "O2":-1, "N2":0, "NXC7H16":11})

    AllModsValid = copy.deepcopy(myTimeStepFull)
    for m in range(3):
        AllModsValid.removeCellsWithValue("TolAirHepZetaQ2DF" + str(m + 1) + "Valid", False, "eqexact")

    Mod2Mod3Valid = copy.deepcopy(myTimeStepFull)
    for m in [1, 2]:
        Mod2Mod3Valid.removeCellsWithValue("TolAirHepZetaQ2DF" + str(m + 1) + "Valid", False, "eqexact")
    Mod2Mod3Valid.removeCellsWithValue("TolAirHepZetaQ2DF1Valid", True, "eqexact")
    Mod2Mod3Valid.addQuantity("BestQ2DFModel23")

    ManiQ2DFValid = copy.deepcopy(myTimeStepFull)
    ManiQ2DFValid.removeCellsWithValue("TolAirHepManiQ2DFValid", False, "eqexact")
    
    DESIRED_NUM_CELLS = 50
    AllModsValid.reduceNumCells(special1 = [["1", "2", "3"], [DESIRED_NUM_CELLS, DESIRED_NUM_CELLS, DESIRED_NUM_CELLS]])
    Mod2Mod3Valid.reduceNumCells(special1 = [["2", "3"], [DESIRED_NUM_CELLS, DESIRED_NUM_CELLS]])
    ManiQ2DFValid.reduceNumCells(special1 = [["1", "2", "3"], [DESIRED_NUM_CELLS, DESIRED_NUM_CELLS, DESIRED_NUM_CELLS]])
    AllModsValid.storeData("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/AllModsValidReduced")
    Mod2Mod3Valid.storeData("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/Mod2Mod3ValidReduced")
    #ManiQ2DFValid.storeData("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/ManiQ2DFValidReduced")

def PDRsRunning():
    AllModsValid = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/AllModsValidReduced"), quants, "000124", False)
    Mod2Mod3Valid = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/Mod2Mod3ValidReduced"), quants, "000124", False)
    #ManiQ2DFValid = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/ManiQ2DFValidReduced"), quants, "000124", False)
    print("Main Func - Stored data retrieved")

    for model in ["1", "2", "3"]:
        AllModsValid.runManiZ("TolAirHep", model, ["Y-CO"], constants.ZETA_INPUT_ZETA_FOLDER, "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/24_11_03FinalZeta", "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/24_11_03FinalZeta")
    for model in ["2", "3"]:
        Mod2Mod3Valid.runManiZ("TolAirHep", model, ["Y-CO"], constants.ZETA_INPUT_ZETA_FOLDER, "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/24_11_03FinalZeta", "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/24_11_03FinalZeta")
    #ManiQ2DFValid.runManiQ2DF("TolAirHep", ["Y-CO"], constants.MANI_Q2DF_INPUT, "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/24_11_03Q2DF", "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/24_11_03Q2DF")
    #ManiQ2DFValid.runManiZ("TolAirHep", "Best", ["Y-CO"], constants.ZETA_INPUT_Q2DF_FOLDER, "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/24_11_03Q2DF", "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/24_11_03Q2DF")
    print("Main Func - PDRs ran")

    AllModsValidMod1YCO = np.array(AllModsValid.getDataList("TolAirHepZetaQ2DF1ResultY-CO"))
    AllModsValidMod2YCO = np.array(AllModsValid.getDataList("TolAirHepZetaQ2DF2ResultY-CO"))
    AllModsValidMod3YCO = np.array(AllModsValid.getDataList("TolAirHepZetaQ2DF3ResultY-CO"))
    Mod2Mod3ValidMod2YCO = np.array(Mod2Mod3Valid.getDataList("TolAirHepZetaQ2DF2ResultY-CO"))
    Mod2Mod3ValidMod3YCO = np.array(Mod2Mod3Valid.getDataList("TolAirHepZetaQ2DF3ResultY-CO"))

    '''print(1, AllModsValidMod1YCO)
    print(2, AllModsValidMod2YCO)
    print(3, AllModsValidMod3YCO)
    print(4, Mod2Mod3ValidMod2YCO)
    print(5, Mod2Mod3ValidMod3YCO)'''

    #ManiQ2DFValidYCO = np.array(ManiQ2DFValid.getDataList("TolAirHepManiQ2DFResultY-CO"))
    #ManiQ2DFValidWithZetaBestYCO = np.array(ManiQ2DFValid.getDataList("TolAirHepZetaQ2DFBestResultY-CO"))

    folderLoc = "/home/efeeroz/Documents/CombustionModelAnalysis/graphs/paperfigs"
    folderLocQ2DFTesting = "/home/efeeroz/Documents/CombustionModelAnalysis/graphs/q2df_testing"

    SHARED_COLORBAR_MIN = -10
    SHARED_COLORBAR_MAX = 0
    utils.makeAndSaveFigure(AllModsValid, "ZMIX", "FMIX", "$Z$", "$F$", "$\log_{10}|Y_{CO}(Q2DF1) - Y_{CO}(Q2DF2)|$", folderLoc + "/AllModsValidMod12.pdf", zList = np.log10(np.abs(AllModsValidMod1YCO - AllModsValidMod2YCO)), dotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"], zMinVal = SHARED_COLORBAR_MIN, zMaxVal = SHARED_COLORBAR_MAX)
    utils.makeAndSaveFigure(AllModsValid, "ZMIX", "FMIX", "$Z$", "$F$", "$\log_{10}|Y_{CO}(Q2DF2) - Y_{CO}(Q2DF3)|$", folderLoc + "/AllModsValidMod23.pdf", zList = np.log10(np.abs(AllModsValidMod2YCO - AllModsValidMod3YCO)), dotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"], zMinVal = SHARED_COLORBAR_MIN, zMaxVal = SHARED_COLORBAR_MAX)
    utils.makeAndSaveFigure(AllModsValid, "ZMIX", "FMIX", "$Z$", "$F$", "$\log_{10}|Y_{CO}(Q2DF3) - Y_{CO}(Q2DF1)|$", folderLoc + "/AllModsValidMod31.pdf", zList = np.log10(np.abs(AllModsValidMod3YCO - AllModsValidMod1YCO)), dotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"], zMinVal = SHARED_COLORBAR_MIN, zMaxVal = SHARED_COLORBAR_MAX)
    utils.makeAndSaveFigure(Mod2Mod3Valid, "ZMIX", "FMIX", "$Z$", "$F$", "$\log_{10}|Y_{CO}(Q2DF2) - Y_{CO}(Q2DF3)|$", folderLoc + "/Mod2Mod3ValidMod23.pdf", zList = np.log10(np.abs(Mod2Mod3ValidMod2YCO - Mod2Mod3ValidMod3YCO)), dotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"], zMinVal = SHARED_COLORBAR_MIN, zMaxVal = SHARED_COLORBAR_MAX)
    #utils.makeAndSaveFigure(ManiQ2DFValid, "ZMIX", "FMIX", "$Z$", "$F$", "$Y_{CO}(\text{ManiQ2DF}) / Y_{CO}(\text{ManiZBest})$", folderLocQ2DFTesting + "/Q2DFTesting.pdf", zList = np.abs(ManiQ2DFValidYCO / ManiQ2DFValidWithZetaBestYCO), dotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"])
    print("Main Func - PDRs results plotted")
    #print("results,", np.abs(ManiQ2DFValidYCO / ManiQ2DFValidWithZetaBestYCO))
    #print(ManiQ2DFValidYCO, ManiQ2DFValidWithZetaBestYCO)

def codeTest():
    '''
    myTimeStepFull = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/ZMIXExtremesTrimmedFullNumPoints"), ALL_QUANTITIES, "000124", False)
    myTimeStepFull.prepareForQ2DFModels("TolAirHep", {"A1CH3":1}, {"O2":0.23292, "N2":0.76708}, {"NXC7H16":1}, 300, 300, 300, {"A1CH3":92.14, "O2":31.999, "N2":28.01, "NXC7H16":100.21}, {"A1CH3":9, "O2":-1, "N2":0, "NXC7H16":11})
    for m in range(3):
        myTimeStepFull.removeCellsWithValue("TolAirHepZetaQ2DF" + str(m + 1) + "Valid", False, "eqexact")
    myTimeStepFull.reduceNumCells(num = 3)
    for i in [1]:
        for key in myTimeStepFull.data.keys():
            print(key + ":", myTimeStepFull.getData(key, i))
            print("$$")
        print("********************\n"*5)

    myTimeStepFull.storeData("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/Nov3CodeTesting")
    '''
    testingTimeStep = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/Nov3CodeTesting"), quants, "000124", False)

    for model in ["1", "2", "3"]:
        testingTimeStep.runManiZ("TolAirHep", model, ["Y-CO"], constants.ZETA_INPUT_ZETA_FOLDER, "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/24_11_03Test", "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/24_11_03Test")

    Mod1YCO = np.array(testingTimeStep.getDataList("TolAirHepZetaQ2DF1ResultY-CO"))
    Mod2YCO = np.array(testingTimeStep.getDataList("TolAirHepZetaQ2DF2ResultY-CO"))
    Mod3YCO = np.array(testingTimeStep.getDataList("TolAirHepZetaQ2DF3ResultY-CO"))
    print(Mod1YCO, "^^", Mod2YCO, "^^", Mod3YCO)

def maniQ2DFTest():
    # myTimeStepFull = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/ZMIXExtremesTrimmedFullNumPoints"), ALL_QUANTITIES, "000124", False)
    # myTimeStepFull.prepareForQ2DFModels("TolAirHep", {"A1CH3":1}, {"O2":0.23292, "N2":0.76708}, {"NXC7H16":1}, 300, 300, 300, {"A1CH3":92.14, "O2":31.999, "N2":28.01, "NXC7H16":100.21}, {"A1CH3":9, "O2":-1, "N2":0, "NXC7H16":11})
    # print("reached here in testing")
    # myTimeStepFull.removeCellsWithValue("TolAirHepManiQ2DFValid", False, "eqexact")
    # myTimeStepFull.reduceNumCells(num = 1)
    # myTimeStepFull.storeData("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/JunkDELETE")

    myTimeStepFull = TimeStep(utils.getInputFiles("/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/JunkDELETE"), quants, "000124", False)
    myTimeStepFull.runManiZ("TolAirHep", "Best", ["Y-CO"], constants.ZETA_INPUT_Q2DF_FOLDER, "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/24_11_03TestQ2DF", "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/24_11_03TestQ2DF")
    myTimeStepFull.runManiQ2DF("TolAirHep", ["Y-CO"], constants.MANI_Q2DF_INPUT, "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/24_11_03TestQ2DF", "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/24_11_03TestQ2DF")

    ManiQ2DFValidYCO = np.array(myTimeStepFull.getDataList("TolAirHepManiQ2DFResultY-CO"))
    ManiQ2DFValidWithZetaBestYCO = np.array(myTimeStepFull.getDataList("TolAirHepZetaQ2DFBestResultY-CO"))
    print(ManiQ2DFValidYCO, "and", ManiQ2DFValidWithZetaBestYCO)


def maniQ2DFNov15Test():
    testStep = TimeStep(None, None, None, data = {'ZSTAR': [0.016853515426432764, 0.025714044423171264, 0.01780163838773446, 0.031518620906190224, 0.012921340854308912, 0.025992279171368583, 0.030726908718502352, 0.049682904057937904, 0.036254013034948816, 0.014373284789867786, 0.001798478138111176, 0.00757325312230353, 0.009400059007875855, 0.015966328074508717, 0.0022600663744325625], 'ZMIX': [0.03744471628034196, 0.03040525123477511, 0.029791240086982196, 0.040834186097218894, 0.031514249384210435, 0.043505587707020066, 0.03389795035859733, 0.049953121724504834, 0.0470427317357196, 0.03332774916251802, 0.046006254276218406, 0.044496089359940205, 0.04788550103966027, 0.04719102589415225, 0.04033988997915683], 'CHI': [6.604099461721531, 6.971846333307038, 4.292228320266545, 5.00884332174878, 5.396007234318317, 4.8964612655053745, 6.50147616247029, 6.993686864591076, 0.9609581060916272, 0.7381924414978014, 5.893947916072884, 5.348756306287921, 4.321499769912338, 8.385242103270665, 3.32368751331803], 'C_ZZST': [-0.08374260926253775, -1.6390599767977783, 0.571072857967963, 0.3291792555792996, 0.19697819076960613, -1.6355934450546479, -0.6226196153527671, 1.1846472649150923, -0.5842824031803963, -0.5312484482687941, -0.01666435888821575, -0.11867556728665551, -0.17155712608898055, 0.3176691576776195, -0.17239512396888368], 'C_ZST': [0.9767324143217997, 0.44845710338059674, 0.7423809627099739, 0.5364336402693972, 0.18663738932530594, 0.8869888516572825, 0.38810022060538973, 0.2647693536193543, 0.39284226004248246, 0.8460863329277816, 0.010462729994124675, 0.008012078900702102, 0.09624210342974526, 0.08505181569400877, 0.0097479235488932]})
    for parm in ['FMIX', 'Q2DF1Xi', 'Q2DF2Xi', 'Q2DF3Xi', 'Q2DF1Eta', 'Q2DF2Eta', 'Q2DF3Eta', 'Q2DF1ChiXi', 'Q2DF2ChiXi', 'Q2DF3ChiXi', 'Q2DF1ChiEta', 'Q2DF2ChiEta', 'Q2DF3ChiEta', 'Q2DF1Ratio', 'Q2DF2Ratio', 'Q2DF3Ratio', 'BestQ2DFModel']:
        testStep.addQuantity(parm)
    testStep.prepareForQ2DFModels("TolAirHep", {"A1CH3":1}, {"O2":0.23292, "N2":0.76708}, {"NXC7H16":1}, 300, 300, 300, {"A1CH3":92.14, "O2":31.999, "N2":28.01, "NXC7H16":100.21}, {"A1CH3":9, "O2":-1, "N2":0, "NXC7H16":11})
    print(testStep.getNumCells(), "before")
    testStep.removeCellsWithValue("TolAirHepManiQ2DFValid", False, "eqexact")
    print(testStep.getNumCells(), "after")
    
    testStep.runManiQ2DF("TolAirHep", ["Y-CO"], constants.MANI_Q2DF_INPUT, "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/24_11_15TESTQ2DF", "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/24_11_15TESTQ2DF")
    testStep.runManiZ("TolAirHep", "Best", ["Y-CO"], constants.ZETA_INPUT_Q2DF_FOLDER, "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/24_11_15TESTQ2DF", "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/24_11_15TESTQ2DF")

    testStepQ2DFYCO = np.array(testStep.getDataList("TolAirHepManiQ2DFResultY-CO"))
    testStepWithZetaBestYCO = np.array(testStep.getDataList("TolAirHepZetaQ2DFBestResultY-CO"))

    #for quant in ["Q2DF1ChiEta", "Q2DF1ChiXi", "Q2DF1Ratio", "Q2DF2ChiEta", "Q2DF2ChiXi", "Q2DF2Ratio", "Q2DF3ChiEta", "Q2DF3ChiXi", "Q2DF3Ratio", "BestQ2DFModel"]:
    #    testStep.addQuantity(quant)
    error = (testStepQ2DFYCO - testStepWithZetaBestYCO) / testStepWithZetaBestYCO * 100
    for i in range(testStep.getNumCells()):
        print(testStep.getData("BestQ2DFModel", i), "Best:", round(error[i], 7), "|  ", round(testStepWithZetaBestYCO[i], 7), "#", round(testStep.getData("TolAirHepZetaQ2DF" + testStep.getData("BestQ2DFModel", i) + "ZValue", i), 7), round(testStep.getData("TolAirHepManiQ2DFZValue", i), 7))
        if abs(error[i]) > 0.1:
            print("\t", testStepWithZetaBestYCO[i], testStepQ2DFYCO[i])
    
    print("Max:", max(np.abs(error)))
    print("^^^")
    print(testStepQ2DFYCO - testStepWithZetaBestYCO)
    print(max(testStepQ2DFYCO - testStepWithZetaBestYCO))

'''
Make sure to change:
- Change DESIRED_NUM_CELLS
- Change the PDRs input and output folders when runManiZ or runManiQ2DF is called, and in constants.py where output directory is specified in the PDRs input
- Consider changing the (two) folderLoc variables or folderLocQ2DFTesting variable to change where plots are saved to
'''
if __name__ == "__main__":
    pass
    maniQ2DFNov15Test()
    #initialDataSetup()
    #trimmingData()
    #paperNonPDRsPlotting("ZMIXExtremesTrimmedFewerPoints")
    #PDRsPreparation()
    #PDRsRunning()
    #codeTest()
    #maniQ2DFTest()
