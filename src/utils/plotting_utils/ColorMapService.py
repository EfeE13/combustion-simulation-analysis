import numpy as np
from typing import Optional
import matplotlib.tri as tri
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib import cm

DECIMAL_PLACES_PRINTING = 4
PRINT_CONTINUOUS_COLOR_MAP = True
CMAP_CONTINUOUS_COLOR_MAP = "RdBu_r"
CMAP_NON_CONTINUOUS_COLOR_MAP = "RdBu_r"
DEFAULT_DOT_SIZE = 15
DOT_BORDER_WIDTH = 1.01
CMAP_SHAPE_ORDER = ['^', 's', '*', 'X', 'v', 'D']

FONT_SIZE = 19

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

    def makeColorPlotNonContinuous(self, axObj, xAxisLabel:str, yAxisLabel:str, titleLabel:str, colorbarLabel:str, increasedDotSize:float = None, 
                                   diffShapes:bool = False, diffShapesIndices:list[list[int]] = None, diffShapesLegend:list[str] = None, 
                                   vMin:float = None, vMax:float = None) -> None:
        """Makes plot of colored dots

        :param _type_ axObj: the axis object to be drawn on
        :param str xAxisLabel: x-axis label
        :param str yAxisLabel: y-axis label
        :param str titleLabel: title label
        :param str colorbarLabel: colorbar label
        :param float increasedDotSize: increased size of dot, defaults to None
        :param bool diffShapes: whether or not to have different dots be different shapes, defaults to False
        :param list[list[int]] diffShapesIndices: each of the sublists contains the indices to be assignned a certain shape I believe, defaults to None
        :param list[str] diffShapesLegend: for making the legend of the different shapes, defaults to None
        :param float vMin: lower cap for coloring, defaults to None
        :param float vMax: upper cap for coloring, defaults to None
        """
        plt.rcParams.update({'font.size': FONT_SIZE})
        
        diffShapesIndicesProvided = diffShapesIndices is not None
        diffShapesLegendProvided = diffShapesLegend is not None
        assert diffShapes == diffShapesIndicesProvided
        assert diffShapes == diffShapesLegendProvided

        xList = [tripleToColorPlot[0] for tripleToColorPlot in self.triplesToColorPlot]
        yList = [tripleToColorPlot[1] for tripleToColorPlot in self.triplesToColorPlot]
        zList = [tripleToColorPlot[2] for tripleToColorPlot in self.triplesToColorPlot]

        if increasedDotSize == None:
            increasedDotSize = 0
        
        if not diffShapes:
            # TODO: remove +10 from dotSize in line below
            norm = plt.Normalize(vmin=vMin, vmax=vMax)  # Normalize z values
            edge_colors = cm.RdBu_r(norm(zList))
            if colorbarLabel == "$x'$":
                edge_colors = cm.twilight_shifted(norm(zList))
            elif colorbarLabel == "Preferred Model Based on $\psi$":
                edge_colors = cm.coolwarm(norm(zList))
            c = axObj.scatter(xList, yList, s = DEFAULT_DOT_SIZE + increasedDotSize, facecolors = "none", edgecolors = edge_colors, linewidths = DOT_BORDER_WIDTH)
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

                norm = plt.Normalize(vmin=vminVal, vmax=vmaxVal)  # Normalize z values
                edge_colors = cm.RdBu_r(norm(zListForShape))
                c = axObj.scatter(xListForShape, yListForShape, marker = CMAP_SHAPE_ORDER[shapeCounter], label = diffShapesLegend[shapeCounter], s = DEFAULT_DOT_SIZE + increasedDotSize, facecolors = "none", edgecolors = edge_colors, linewidths = DOT_BORDER_WIDTH)
                shapeCounter += 1
            axObj.legend()

        pos = axObj.get_position()
        cax = axObj.get_figure().add_axes([pos.x1 + 0.01, pos.y0, 0.01, pos.height])
        
        cmapToUse = CMAP_NON_CONTINUOUS_COLOR_MAP
        if colorbarLabel == "$x'$":
            cmapToUse = "twilight_shifted"
        cbar = axObj.get_figure().colorbar(cm.ScalarMappable(norm = norm, cmap = cmapToUse), cax = cax)
        if colorbarLabel != None:
            cbar.set_label(colorbarLabel)

        axObj.set_xlabel(xAxisLabel)
        axObj.set_ylabel(yAxisLabel)
        if titleLabel != None:
            axObj.set_title(titleLabel)

        # Fixes overlapping x-tick labels for the base-10 log of abs. difference in Q2DF2 and Q2DF3
        # predictions where both models are valid but Q2DF1 is not:
        if colorbarLabel == "$\log_{10}|Y_{CO}(Q2DF2) - Y_{CO}(Q2DF3)|$ ":
            xTickLabels = ["0.07", "0.075", "0.08", "0.085"]
            axObj.set_xticks([float(elm) for elm in xTickLabels])
            axObj.set_xticklabels(xTickLabels)