import copy
import numpy as np
import random
import textwrap
from scipy.special import erfcinv

from src.nga_analysis import constants
from src.nga_analysis import utils
from src.nga_analysis.TimeStep import TimeStep
from src.utils.plotting_utils.DataPlot import DataPlot
from src.nga_analysis.Q2DFGeneralizer import Q2DFGeneralizer
from src.nga_analysis.OnTheFlyAndQ2DF2Runner import OnTheFlyAndQ2DF2Runner

ENSIGHT_FOLDER = "/home/efeeroz/Documents/CombustionModelAnalysis/ensight_input/after_isat_bugfix"

LABEL  = constants.LABEL
FOLDER_FIGS = "/home/efeeroz/Documents/CombustionModelAnalysis/graphs/" + LABEL + "paperfigs"
FOLDER_Q2DF_TESTING_FIG = "/home/efeeroz/Documents/CombustionModelAnalysis/graphs/" + LABEL + "q2df_testing"

STORE_FULL_SQUARE_FULL_POINTS = "/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/" + LABEL + "FullSquareFullNumPoints"
STORE_ZMIX_TRIM_FULL_POINTS = "/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/" + LABEL + "ZMIXExtremesTrimmedFullNumPoints"
STORE_ZMIX_TRIM_FEWER_POINTS = "/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/" + LABEL + "ZMIXExtremesTrimmedFewerPoints"
STORE_ALL_MODELS_REDUCED = "/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/" + LABEL + "AllModsValidReduced"
STORE_MOD2_MOD3_REDUCED = "/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/" + LABEL + "Mod2Mod3ValidReduced"
STORE_MANIQ2DF_REDUCED = "/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/" + LABEL + "ManiQ2DFValidReduced"
STORE_IDAM_DATA = "/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/" + LABEL + "IDAMData"

INPUT_FOLDER_ZETA = "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/" + LABEL + "Zeta"
OUTPUT_FOLDER_ZETA = "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/" + LABEL + "Zeta"
INPUT_FOLDER_Q2DF = "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/" + LABEL + "Q2DF"
OUTPUT_FOLDER_Q2DF = "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/" + LABEL + "Q2DF"

ALL_QUANTITIES = ['ZSTAR', 'ZMIX', 'Temp', 'C_ZST', 'C_ZZST', 'CHI', 'ADDS', 'DIFF', 'P', 'VISC', 'Y_CO2', 'Y_H2', 'Y_O2', 'ZSTA2', 'CHIREF', 'dRHO', 'KSGS', 'RHO', 'Y_CO', 'Y_F', 'Y_H2O', 'Y_OH', 'ZMIX2', 'Q2DF1ChiEta', 'Q2DF2ChiEta', 'Q2DF3ChiEta', 'Q2DF1ChiXi', 'Q2DF2ChiXi', 'Q2DF3ChiXi', 'Q2DF1Eta', 'Q2DF2Eta', 'Q2DF3Eta', 'FMIX', 'Q2DF1Ratio', 'Q2DF2Ratio', 'Q2DF3Ratio', 'Q2DF1Xi', 'Q2DF2Xi', 'Q2DF3Xi', "Q2DF1RatioNorm", "Q2DF2RatioNorm", "Q2DF3RatioNorm", "Ratio12Norm", "Ratio23Norm", "Ratio31Norm", "BestQ2DFModel"]
quants = copy.deepcopy(ALL_QUANTITIES)
quants += ["TolAirHepManiQ2DFFileName", "TolAirHepManiQ2DFZValue", "TolAirHepManiQ2DFValid"]
for model in ["1", "2", "3"]:
    quants += ["TolAirHepZetaQ2DF" + model + "ZValue", "TolAirHepZetaQ2DF" + model + "OxidTemp", "TolAirHepZetaQ2DF" + model + "FuelTemp", "TolAirHepZetaQ2DF" + model + "OxidBC", "TolAirHepZetaQ2DF" + model + "FuelBC", "TolAirHepZetaQ2DF" + model + "FuelSpecies", "TolAirHepZetaQ2DF" + model + "ScalarDissRate", "TolAirHepZetaQ2DF" + model + "FileName", "TolAirHepZetaQ2DF" + model + "Valid"]

def initialDataSetup():
    allQuantities = []
    myTimeStep = TimeStep(utils.getInputFiles(ENSIGHT_FOLDER), ["ZSTAR", "ZMIX", "Temp", "C_ZST", "C_ZZST", "CHI", 'ADDS', 'DIFF', 'P', 'VISC', 'Y_CO2', 'Y_H2', 'Y_O2', 'ZSTA2', 'CHIREF', 'dRHO', 'KSGS', 'RHO', 'Y_CO', 'Y_F', 'Y_H2O', 'Y_OH', 'ZMIX2'], "000123", True)
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

    myTimeStep.storeData(STORE_FULL_SQUARE_FULL_POINTS)

def trimmingData():
    myTimeStepFull = TimeStep(utils.getInputFiles(STORE_FULL_SQUARE_FULL_POINTS), ALL_QUANTITIES, "000123", False)
    assert myTimeStepFull.getNumCells() == 45206   # 46135
    myTimeStepFull.removeCellsWithValue("ZMIX", 0.0001, "lt")
    myTimeStepFull.removeCellsWithValue("ZMIX", 0.9999, "gt")
    print("Full Size:", myTimeStepFull.getNumCells())
    myTimeStepFull.storeData(STORE_ZMIX_TRIM_FULL_POINTS)

    myTimeStepFullToReduced = copy.deepcopy(myTimeStepFull)
    myTimeStepFullToReduced.reduceNumCells(num = 1000)
    myTimeStepFullToReduced.storeData(STORE_ZMIX_TRIM_FEWER_POINTS)

def paperNonPDRsPlotting():
    myTimeStep = TimeStep(utils.getInputFiles(STORE_ZMIX_TRIM_FEWER_POINTS), ALL_QUANTITIES, "000123", False)
    assert myTimeStep.getNumCells() == 1000

    folderLoc = FOLDER_FIGS
    utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "Temperature (K)", folderLoc + "/Temperature.pdf", zVar = "Temp")
    utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "Preferred Model Based on $\psi$", folderLoc + "/BestQ2DFModel.pdf", zMinVal = 1.0, zMaxVal = 3.0, zVar = "BestQ2DFModel")
    for i in [1, 2, 3]:
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "$\psi_" + str(i) + " / (\psi_" + str(i) + " + 1)$", folderLoc + "/Psi" + str(i) + "Norm.pdf", zVar = "Q2DF" + str(i) + "RatioNorm")
    for compareValues in [[1, 2], [2, 3], [3, 1]]:
        i = str(compareValues[0]); j = str(compareValues[1])
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "$\psi_" + i + "/\psi_" + j + "$ Normalized", folderLoc + "/Ratio" + i + j + "Norm.pdf", zVar = "Ratio" + i + j + "Norm")
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "$\psi_" + i + "/\psi_" + j + "$ Normalized", folderLoc + "/Ratio" + i + j + "NormCutoff.pdf", zVar = "Ratio" + i + j + "Norm", zMinVal = 0.0, zMaxVal = 1.0)
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "$\psi_" + i + "/\psi_" + j + "$ Normalized", folderLoc + "/Ratio" + i + j + "NormDiffShapes.pdf", zVar = "Ratio" + i + j + "Norm", zMinVal = 0.0, zMaxVal = 1.0, increasedDotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"])
    
    myTimeStep.addQuantity("xPrime")
    utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", "$x'$", folderLoc + "/xPrimeVal.pdf", titleLabel = "$x'$ Value", zMinVal = 0.0, zMaxVal = 6.0, zVar = "xPrime")

    myTimeStepForX = copy.deepcopy(myTimeStep)
    myTimeStepForX.addQuantity("BestMVal")
    myTimeStepForX.removeCellsWithValue("", -10000.0, "validm")
    myTimeStepForX.addQuantity("BestXVal")
    myTimeStepForX.addQuantity("BestXValNormalized")
    utils.makeAndSaveFigure(myTimeStepForX, "ZMIX", "FMIX", "$Z$", "$F$", "$x$", folderLoc + "/xValNorm.pdf", titleLabel = "$x$ Value", zMinVal = 0.0, zMaxVal = 1.0, zVar = "BestXValNormalized")

def PDRsPreparation():
    myTimeStepFull = TimeStep(utils.getInputFiles(STORE_ZMIX_TRIM_FULL_POINTS), ALL_QUANTITIES, "000123", False)
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
    AllModsValid.storeData(STORE_ALL_MODELS_REDUCED)
    Mod2Mod3Valid.storeData(STORE_MOD2_MOD3_REDUCED)
    ManiQ2DFValid.storeData(STORE_MANIQ2DF_REDUCED)

def PDRsRunning():
    AllModsValid = TimeStep(utils.getInputFiles(STORE_ALL_MODELS_REDUCED), quants, "000123", False)
    Mod2Mod3Valid = TimeStep(utils.getInputFiles(STORE_MOD2_MOD3_REDUCED), quants, "000123", False)
    ManiQ2DFValid = TimeStep(utils.getInputFiles(STORE_MANIQ2DF_REDUCED), quants, "000123", False)
    print("Main Func - Stored data retrieved")

    for model in ["1", "2", "3"]:
        AllModsValid.runManiZ("TolAirHep", model, ["Y-CO"], constants.ZETA_INPUT_ZETA_FOLDER, INPUT_FOLDER_ZETA, OUTPUT_FOLDER_ZETA)
    for model in ["2", "3"]:
        Mod2Mod3Valid.runManiZ("TolAirHep", model, ["Y-CO"], constants.ZETA_INPUT_ZETA_FOLDER, INPUT_FOLDER_ZETA, OUTPUT_FOLDER_ZETA)
    #ManiQ2DFValid.runManiQ2DF("TolAirHep", ["Y-CO"], constants.MANI_Q2DF_INPUT, INPUT_FOLDER_Q2DF, OUTPUT_FOLDER_Q2DF)
    #ManiQ2DFValid.runManiZ("TolAirHep", "Best", ["Y-CO"], constants.ZETA_INPUT_Q2DF_FOLDER, INPUT_FOLDER_Q2DF, OUTPUT_FOLDER_Q2DF)
    print("Main Func - PDRs ran")

    AllModsValidMod1YCO = np.array(AllModsValid.getDataList("TolAirHepZetaQ2DF1ResultY-CO"))
    AllModsValidMod2YCO = np.array(AllModsValid.getDataList("TolAirHepZetaQ2DF2ResultY-CO"))
    AllModsValidMod3YCO = np.array(AllModsValid.getDataList("TolAirHepZetaQ2DF3ResultY-CO"))
    Mod2Mod3ValidMod2YCO = np.array(Mod2Mod3Valid.getDataList("TolAirHepZetaQ2DF2ResultY-CO"))
    Mod2Mod3ValidMod3YCO = np.array(Mod2Mod3Valid.getDataList("TolAirHepZetaQ2DF3ResultY-CO"))

    #ManiQ2DFValidYCO = np.array(ManiQ2DFValid.getDataList("TolAirHepManiQ2DFResultY-CO"))
    #ManiQ2DFValidWithZetaBestYCO = np.array(ManiQ2DFValid.getDataList("TolAirHepZetaQ2DFBestResultY-CO"))

    folderLoc = FOLDER_FIGS

    SHARED_COLORBAR_MIN = -10
    SHARED_COLORBAR_MAX = 0
    utils.makeAndSaveFigure(AllModsValid, "ZMIX", "FMIX", "$Z$", "$F$", "$\log_{10}|Y_{CO}(Q2DF1) - Y_{CO}(Q2DF2)|$", folderLoc + "/AllModsValidMod12.pdf", zList = np.log10(np.abs(AllModsValidMod1YCO - AllModsValidMod2YCO)), increasedDotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"], zMinVal = SHARED_COLORBAR_MIN, zMaxVal = SHARED_COLORBAR_MAX)
    utils.makeAndSaveFigure(AllModsValid, "ZMIX", "FMIX", "$Z$", "$F$", "$\log_{10}|Y_{CO}(Q2DF2) - Y_{CO}(Q2DF3)|$", folderLoc + "/AllModsValidMod23.pdf", zList = np.log10(np.abs(AllModsValidMod2YCO - AllModsValidMod3YCO)), increasedDotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"], zMinVal = SHARED_COLORBAR_MIN, zMaxVal = SHARED_COLORBAR_MAX)
    utils.makeAndSaveFigure(AllModsValid, "ZMIX", "FMIX", "$Z$", "$F$", "$\log_{10}|Y_{CO}(Q2DF3) - Y_{CO}(Q2DF1)|$", folderLoc + "/AllModsValidMod31.pdf", zList = np.log10(np.abs(AllModsValidMod3YCO - AllModsValidMod1YCO)), increasedDotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"], zMinVal = SHARED_COLORBAR_MIN, zMaxVal = SHARED_COLORBAR_MAX)
    utils.makeAndSaveFigure(Mod2Mod3Valid, "ZMIX", "FMIX", "$Z$", "$F$", "$\log_{10}|Y_{CO}(Q2DF2) - Y_{CO}(Q2DF3)|$ ", folderLoc + "/Mod2Mod3ValidMod23.pdf", zList = np.log10(np.abs(Mod2Mod3ValidMod2YCO - Mod2Mod3ValidMod3YCO)), increasedDotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"], zMinVal = SHARED_COLORBAR_MIN, zMaxVal = SHARED_COLORBAR_MAX)
    #utils.makeAndSaveFigure(ManiQ2DFValid, "ZMIX", "FMIX", "$Z$", "$F$", "$Y_{CO}(\text{ManiQ2DF}) / Y_{CO}(\text{ManiZBest})$", FOLDER_Q2DF_TESTING_FIG + "/Q2DFTesting.pdf", zList = np.abs(ManiQ2DFValidYCO / ManiQ2DFValidWithZetaBestYCO), increasedDotSize = 4, diffShapes = True, diffShapesLegend = ["Q2DF1 Best", "Q2DF2 Best", "Q2DF3 Best"])
    print("Main Func - PDRs results plotted")
    #print("results,", np.abs(ManiQ2DFValidYCO / ManiQ2DFValidWithZetaBestYCO))
    #print(ManiQ2DFValidYCO, ManiQ2DFValidWithZetaBestYCO)

def RawMixingDataIDAM():
    myTimeStep = TimeStep(utils.getInputFiles(STORE_ZMIX_TRIM_FEWER_POINTS), ALL_QUANTITIES, "000123", False)
    NUM_CELLS = 1000
    assert myTimeStep.getNumCells() == NUM_CELLS

    cell_order = list(range(NUM_CELLS))
    random.seed(13)
    random.shuffle(cell_order)

    raw_mixing_data = []
    for cell_num in cell_order:
        Z1 = myTimeStep.getData("ZSTAR", cell_num)
        Z2 = 1 - myTimeStep.getData("ZMIX", cell_num)
        chi11 = myTimeStep.getData("C_ZST", cell_num)
        chi12 = -myTimeStep.getData("C_ZZST", cell_num)
        chi22 = myTimeStep.getData("CHI", cell_num)
        raw_mixing_data.append([Z1, Z2, chi11, chi12, chi22])

    raw_mixing_data_1D = []
    for i in range(5):
        raw_mixing_data_1D += [raw_mixing_data[j][i] for j in range(NUM_CELLS)]
    raw_mixing_data_1D = [str(elm) + "_WP" for elm in raw_mixing_data_1D]

    str_to_print = "raw_mixing_data = reshape([" + ", ".join(raw_mixing_data_1D) + "], shape(raw_mixing_data))"
    chunks = textwrap.wrap(str_to_print, width = 80, break_long_words = False)
    str_to_print = " &\n\t".join(chunks)
    with open("/home/efeeroz/Other/IDAM_raw_data/IDAM_raw_data.txt", "w") as raw_data_file:
        raw_data_file.write(str_to_print)
    
    cell_shuffling_dict = {}
    for elm in range(1, 1000+1):
        cell_shuffling_dict[elm] = cell_order[elm - 1]
    print(cell_shuffling_dict)

def plotIDAMOutput():
    myTimeStep = TimeStep(utils.getInputFiles(STORE_ZMIX_TRIM_FEWER_POINTS), ALL_QUANTITIES, "000123", False)
    assert myTimeStep.getNumCells() == 1000
    for quantity in ['IDAM_RHO', 'IDAM_DIFF', 'IDAM_VISC', 'IDAM_Temp', 'IDAM_Y_CO2', 'IDAM_Y_O2', 'IDAM_Y_CO', 'IDAM_Y_H2O', 'IDAM_Y_H2', 'IDAM_Y_OH']:
        myTimeStep.addQuantity(quantity)
    
    myTimeStep.removeCellsWithValue("IDAM_DIFF", "NA", "eqexact")
    assert myTimeStep.getNumCells() == 727

    for quantity in ['IDAM_RHO_diff', 'IDAM_DIFF_diff', 'IDAM_VISC_diff', 'IDAM_Temp_diff', 'IDAM_Y_CO2_diff', 'IDAM_Y_O2_diff', 'IDAM_Y_CO_diff', 'IDAM_Y_H2O_diff', 'IDAM_Y_H2_diff', 'IDAM_Y_OH_diff']:
        myTimeStep.addQuantity(quantity)

    for quantity in ['IDAM_RHO_percent', 'IDAM_DIFF_percent', 'IDAM_VISC_percent', 'IDAM_Temp_percent', 'IDAM_Y_CO2_percent', 'IDAM_Y_O2_percent', 'IDAM_Y_CO_percent', 'IDAM_Y_H2O_percent', 'IDAM_Y_H2_percent', 'IDAM_Y_OH_percent']:
        myTimeStep.addQuantity(quantity)
    
    myTimeStep.storeData(STORE_IDAM_DATA)

    folderLoc = FOLDER_FIGS

    labelHelper = {'IDAM_RHO':"rho", 'IDAM_DIFF':"diff", 'IDAM_VISC':"visc", 'IDAM_Temp':"$T$", 'IDAM_Y_CO2':"$Y_{CO2}$", 'IDAM_Y_O2':"$Y_{O2}$", 'IDAM_Y_CO':"$Y_{CO}$", 'IDAM_Y_H2O':"$Y_{H2O}$", 'IDAM_Y_H2':"$Y_{H2}$", 'IDAM_Y_OH':"$Y_{OH}$"}
    for quantity in ['IDAM_RHO', 'IDAM_DIFF', 'IDAM_VISC', 'IDAM_Temp', 'IDAM_Y_CO2', 'IDAM_Y_O2', 'IDAM_Y_CO', 'IDAM_Y_H2O', 'IDAM_Y_H2', 'IDAM_Y_OH']:
        figureLabel = "IDAM " + labelHelper[quantity]
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", figureLabel, folderLoc + "/IDAM_outputs/" + quantity + ".pdf", zVar = quantity)

    labelHelper = {'IDAM_RHO_diff':"rho", 'IDAM_DIFF_diff':"diff", 'IDAM_VISC_diff':"visc", 'IDAM_Temp_diff':"T", 'IDAM_Y_CO2_diff':"Y_{CO2}", 'IDAM_Y_O2_diff':"Y_{O2}", 'IDAM_Y_CO_diff':"Y_{CO}", 'IDAM_Y_H2O_diff':"Y_{H2O}", 'IDAM_Y_H2_diff':"Y_{H2}", 'IDAM_Y_OH_diff':"Y_{OH}"}
    for quantity in ['IDAM_RHO_diff', 'IDAM_DIFF_diff', 'IDAM_VISC_diff', 'IDAM_Temp_diff', 'IDAM_Y_CO2_diff', 'IDAM_Y_O2_diff', 'IDAM_Y_CO_diff', 'IDAM_Y_H2O_diff', 'IDAM_Y_H2_diff', 'IDAM_Y_OH_diff']:
        figureLabel = "$\log|" + labelHelper[quantity] + "($IDAM$) - " + labelHelper[quantity] + "($Q2DF2$)|$"
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", figureLabel, folderLoc + "/IDAM_outputs/log_diff_" + quantity + ".pdf", zVar = quantity)
    
    labelHelper = {'IDAM_RHO_percent':"rho", 'IDAM_DIFF_percent':"diff", 'IDAM_VISC_percent':"visc", 'IDAM_Temp_percent':"T", 'IDAM_Y_CO2_percent':"Y_{CO2}", 'IDAM_Y_O2_percent':"Y_{O2}", 'IDAM_Y_CO_percent':"Y_{CO}", 'IDAM_Y_H2O_percent':"Y_{H2O}", 'IDAM_Y_H2_percent':"Y_{H2}", 'IDAM_Y_OH_percent':"Y_{OH}"}
    for quantity in ['IDAM_RHO_percent', 'IDAM_DIFF_percent', 'IDAM_VISC_percent', 'IDAM_Temp_percent', 'IDAM_Y_CO2_percent', 'IDAM_Y_O2_percent', 'IDAM_Y_CO_percent', 'IDAM_Y_H2O_percent', 'IDAM_Y_H2_percent', 'IDAM_Y_OH_percent']:
        figureLabel = "$\log$ of % diff for $" + labelHelper[quantity] + "$"
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", figureLabel, folderLoc + "/IDAM_outputs/percent_" + quantity + ".pdf", zVar = quantity)

def plotNonIDAMOutput():
    myTimeStep = TimeStep(utils.getInputFiles(STORE_ZMIX_TRIM_FEWER_POINTS), ALL_QUANTITIES, "000123", False)
    assert myTimeStep.getNumCells() == 1000
    folderLoc = FOLDER_FIGS
    labelHelper = {'RHO':"rho", 'DIFF':"diff", 'VISC':"visc", 'Temp':"$T$", 'Y_CO2':"$Y_{CO2}$", 'Y_O2':"$Y_{O2}$", 'Y_CO':"$Y_{CO}$", 'Y_H2O':"$Y_{H2O}$", 'Y_H2':"$Y_{H2}$", 'Y_OH':"$Y_{OH}$"}
    for quantity in ['RHO', 'DIFF', 'VISC', 'Temp', 'Y_CO2', 'Y_O2', 'Y_CO', 'Y_H2O', 'Y_H2', 'Y_OH']:
        figureLabel = "Q2DF Model 2 " + labelHelper[quantity]
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", figureLabel, folderLoc + "/IDAM_outputs/" + "Q2DF2_" + quantity + ".pdf", zVar = quantity)

def plotQ2DF2AndIDAMChiRatios():
    new_storage_location = "/home/efeeroz/Documents/CombustionModelAnalysis/data_storage/" + LABEL + "ZMIXExtremesTrimmedFewerPointsWithXPrime"
    if True:
        myTimeStep = TimeStep(utils.getInputFiles(STORE_ZMIX_TRIM_FEWER_POINTS), ALL_QUANTITIES, "000123", False)
        assert myTimeStep.getNumCells() == 1000
        myTimeStep.addQuantity("xPrime")
        #myTimeStep.storeData(new_storage_location)
    else:
        myTimeStep = TimeStep(utils.getInputFiles(new_storage_location), ALL_QUANTITIES + ["xPrime", "xPrimeChiRatioNorm"], "000123", False)
        assert myTimeStep.getNumCells() == 1000

    for elm in range(1000):
        if not (myTimeStep.getData("Q2DF2RatioNorm", elm) >= myTimeStep.getData("xPrimeChiRatioNorm", elm)):
            print(elm, myTimeStep.getData("Q2DF2RatioNorm", elm), myTimeStep.getData("xPrimeChiRatioNorm", elm))

    folderLoc = FOLDER_FIGS
    labelHelper = {"Q2DF2RatioNorm":"Q2DF Model 2", "xPrimeChiRatioNorm":"On-the-Fly"}
    for quantity in ["Q2DF2RatioNorm", "xPrimeChiRatioNorm"]:
        figureLabel = "$\psi / (\psi + 1)$ for " + labelHelper[quantity]
        utils.makeAndSaveFigure(myTimeStep, "ZMIX", "FMIX", "$Z$", "$F$", figureLabel, folderLoc + "/Chi_Ratios/" + quantity + ".pdf", zVar = quantity, zMinVal = 0, zMaxVal = 1)

def runOnFlyAndQ2DF2():
    myTimeStep = TimeStep(utils.getInputFiles(STORE_ZMIX_TRIM_FEWER_POINTS), ALL_QUANTITIES, "000123", False)
    assert myTimeStep.getNumCells() == 1000
    myTimeStep.addQuantity("5BasicVariables")

    myTimeStep.shuffleData(constants.shuffling_timestep)
    
    my_runner = OnTheFlyAndQ2DF2Runner(myTimeStep, "toluene", "air", "n-heptane")
    on_fly_data_dict, q2df_data_dict = my_runner.run_pdrs("/home/efeeroz/Documents/OnFlyVersusQ2DF/inputs_pdrs", "/home/efeeroz/Documents/OnFlyVersusQ2DF/outputs_pdrs", "/home/efeeroz/Documents/OnFlyVersusQ2DF/log.txt")

'''
Make sure to change:
- Change DESIRED_NUM_CELLS
- Change the folders at the top of the program, and in constants.py where output directory is specified in the PDRs input
- Make sure unneeded function calls in main are commented out to prevent overriding of randomly-selected data
'''
if __name__ == "__main__":
    #initialDataSetup()
    #trimmingData()
    #paperNonPDRsPlotting()
    #PDRsPreparation()
    #PDRsRunning()

    #RawMixingDataIDAM()
    #plotIDAMOutput()

    #plotIDAMOutput()
    #plotNonIDAMOutput()

    #getManualOnTheFlyTestingPoints()
    runOnFlyAndQ2DF2()