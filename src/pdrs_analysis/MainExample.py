from src.pdrs_analysis.PDRsRuns2D import PDRsRuns2D
from src.utils.plotting_utils.DataPlot import DataPlot
from src.pdrs_analysis.PDRsRunnable import PDRsRunnable
from src.pdrs_analysis import constants

def main():
    inputFileSpecs = constants.VARY_SDR_FMX
    inputFolderLoc = "/home/efeeroz/Documents/CombustionModelAnalysis/inputs_pdrs/pdrs_cartesian_product/24-10-15"
    outputFolderLoc = "/home/efeeroz/Documents/CombustionModelAnalysis/outputs_pdrs/pdrs_cartesian_product/24-10-15"

    myPDRsRunnable = PDRsRunnable(inputFileSpecs)
    outputFileLocs = myPDRsRunnable.runPDRs(inputFolderLoc, outputFolderLoc)
    myPDRsRuns2D = PDRsRuns2D(outputFileLocs, "SDR", "FMX")

    dataVis = DataPlot(3, 6)
    for quantityName in "Y-SXCH2 Y-CO Y-H T[K]".strip().split():
        axObj = dataVis.getFreeSubplot()
        myPDRsRuns2D.makeColorPlot(axObj, quantityName, "max", "Scalar Dissipation Rate (1/s))", "FMIX", quantityName + " Color Map", quantityName, plotLogX = True)
    dataVis.saveFigure("/home/efeeroz/Documents/CombustionModelAnalysis/graphs/misc/24-10-15TestYCOColorMap", 500)

if __name__ == "__main__":
    main()