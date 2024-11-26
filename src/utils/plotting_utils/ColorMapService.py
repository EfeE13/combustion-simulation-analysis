import numpy as np
from typing import Optional
import matplotlib.tri as tri
import matplotlib.colors as mcolors

DECIMAL_PLACES_PRINTING = 4
PRINT_CONTINUOUS_COLOR_MAP = True
CMAP_CONTINUOUS_COLOR_MAP = "RdBu_r"
CMAP_NON_CONTINUOUS_COLOR_MAP = "RdBu_r"
DEFAULT_DOT_SIZE = 2
CMAP_SHAPE_ORDER = ['^', 's', '*', 'X', 'v', 'D']

class ColorMapService:
    def __init__(self, triplesToColorPlot:list[list[float]]) -> None:
        """Holds data for continuous or noncontinuous (i.e., with dots) color plotting

        :param list[list[float]] triplesToColorPlot: a list with sublists of length 3, where: 1) The 1st element is for x-axis values (e.g., SDR, scalar dissipation rate, used for PDRs run), 2) The 2nd element is for y-axis values (e.g., FMX used for PDRs run), 3) The 3rd element is what to color by (e.g., the max Y-CO attained in the PDRs run with specified SDR and FMX)
        """
        self.triplesToColorPlot = triplesToColorPlot
        self.xVals = []; self.yVals = []; self.zVals = []
        for tripleList in self.triplesToColorPlot:
            self.xVals.append(tripleList[0])
            self.yVals.append(tripleList[1])
            self.zVals.append(tripleList[2])

    def makeColorPlotContinuous(self, axObj, xAxisLabel:str, yAxisLabel:str, titleLabel:str, colorbarLabel:str) -> None:
        """Makes a color plot where the entire plot is colored (i.e., not colored dots)

        :param _type_ axObj: the axis object to be drawn on
        :param str xAxisLabel: x-axis label
        :param str yAxisLabel: y-axis label
        :param str titleLabel: title label
        :param str colorbarLabel: colorbar label
        """
        assert len(set(self.xVals)) * len(set(self.yVals)) == len(self.zVals)
        uniqXVals = list(sorted(list(set(self.xVals))))
        uniqYVals = list(sorted(list(set(self.yVals))))
        
        X = np.array([uniqXVals] * len(uniqYVals))
        Y = np.transpose(np.array([uniqYVals] * len(uniqXVals)))
        assert X.shape == Y.shape
        Z = np.zeros(X.shape)

        ZHelperDict = {(self.xVals[i], self.yVals[i]):self.zVals[i] for i in range(len(self.xVals))}
        for i in range(len(uniqYVals)):
            for j in range(len(uniqXVals)):
                Z[i][j] = ZHelperDict[(uniqXVals[j], uniqYVals[i])]

        if PRINT_CONTINUOUS_COLOR_MAP:
            print("\n" + "- " * 40)
            for i in range(len(uniqYVals)):
                print("\n" + str(round(uniqYVals[len(uniqYVals) - 1 - i], DECIMAL_PLACES_PRINTING)), end = "")
                for j in range(len(uniqXVals)):
                    print("\t" + str(round(Z[len(uniqYVals) - 1 - i][j], DECIMAL_PLACES_PRINTING)), end = "")

            print()
            for xVal in uniqXVals:
                print("\t" + str(round(xVal, DECIMAL_PLACES_PRINTING)), end = "")
            print("\n" + "- " * 40)

        c = axObj.pcolormesh(X, Y, Z, shading = "nearest", cmap = CMAP_CONTINUOUS_COLOR_MAP)
        pos = axObj.get_position()
        cax = axObj.get_figure().add_axes([pos.x1 + 0.01, pos.y0, 0.01, pos.height])
        cbar = axObj.get_figure().colorbar(c, cax=cax)
        if colorbarLabel != None:
            cbar.set_label(colorbarLabel)

        axObj.set_xlabel(xAxisLabel)
        axObj.set_ylabel(yAxisLabel)
        if titleLabel != None:
            axObj.set_title(titleLabel)

    def makeColorPlotNonContinuous(self, axObj, xAxisLabel:str, yAxisLabel:str, titleLabel:str, colorbarLabel:str, dotSize:float = None, 
                                   diffShapes:bool = False, diffShapesIndices:list[list[int]] = None, diffShapesLegend:list[str] = None, 
                                   vMin:float = None, vMax:float = None) -> None:
        """Makes plot of colored dots

        :param _type_ axObj: the axis object to be drawn on
        :param str xAxisLabel: x-axis label
        :param str yAxisLabel: y-axis label
        :param str titleLabel: title label
        :param str colorbarLabel: colorbar label
        :param float dotSize: size of dot, defaults to None
        :param bool diffShapes: whether or not to have different dots be different shapes, defaults to False
        :param list[list[int]] diffShapesIndices: each of the sublists contains the indices to be assignned a certain shape I believe, defaults to None
        :param list[str] diffShapesLegend: for making the legend of the different shapes, defaults to None
        :param float vMin: lower cap for coloring, defaults to None
        :param float vMax: upper cap for coloring, defaults to None
        """
        
        diffShapesIndicesProvided = diffShapesIndices is not None
        diffShapesLegendProvided = diffShapesLegend is not None
        assert diffShapes == diffShapesIndicesProvided
        assert diffShapes == diffShapesLegendProvided

        xList = [tripleToColorPlot[0] for tripleToColorPlot in self.triplesToColorPlot]
        yList = [tripleToColorPlot[1] for tripleToColorPlot in self.triplesToColorPlot]
        zList = [tripleToColorPlot[2] for tripleToColorPlot in self.triplesToColorPlot]

        if dotSize == None:
            dotSize = DEFAULT_DOT_SIZE
        
        if not diffShapes:
            c = axObj.scatter(xList, yList, c = zList, s = dotSize, cmap = CMAP_NON_CONTINUOUS_COLOR_MAP, vmin = vMin, vmax = vMax)
        else:
            if vMin != None:
                vminVal = vMin
            else:
                vminVal = min(list(zList))
            if vMax != None:
                vmaxVal = vMax
            else:
                vmaxVal = max(list(zList))

            shapeCounter = 0
            #norm = mcolors.Normalize(vmin = min(list(zList)), vmax = max(list(zList)))
            for shapeIndices in diffShapesIndices:
                xListForShape = [xList[i] for i in shapeIndices]
                yListForShape = [yList[i] for i in shapeIndices]
                zListForShape = [zList[i] for i in shapeIndices]

                c = axObj.scatter(xListForShape, yListForShape, c = zListForShape, marker = CMAP_SHAPE_ORDER[shapeCounter], label = diffShapesLegend[shapeCounter], s = dotSize, cmap = CMAP_NON_CONTINUOUS_COLOR_MAP, vmin = vminVal, vmax = vmaxVal)
                shapeCounter += 1
            axObj.legend()

        pos = axObj.get_position()
        cax = axObj.get_figure().add_axes([pos.x1 + 0.01, pos.y0, 0.01, pos.height])
        cbar = axObj.get_figure().colorbar(c, cax=cax)
        if colorbarLabel != None:
            cbar.set_label(colorbarLabel)

        axObj.set_xlabel(xAxisLabel)
        axObj.set_ylabel(yAxisLabel)
        if titleLabel != None:
            axObj.set_title(titleLabel)
