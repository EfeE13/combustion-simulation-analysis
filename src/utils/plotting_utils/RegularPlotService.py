import numpy as np
import matplotlib.pyplot as plt

DECIMAL_PLACES_PRINTING = 4
DEFAULT_DOT_SIZE = 15
DOT_BORDER_WIDTH = 1.01
SHAPE_ORDER = ['^', 's', '*', 'X', 'v', 'D']
COLOR_ORDER = ["red", "orange", "yellow", "green", "blue"]

FONT_SIZE = 19

class RegularPlotService:
    def __init__(self, xVals:list[float], yVals:list[list[float]]) -> None:
        """Holds data for plotting 1D data, including maybe multiple lines

        :param list[float] xVals: a list of the x values
        :param list[list[float]] yVals: a list with sublists where each sublist represents a line (or otherwise 1D data) to be plotted
        """
        self.xVals = xVals
        self.yVals = yVals

    def makeRegularPlot(self, axObj, xAxisLabel:str, yAxisLabel:str, titleLabel:str) -> None:
        """Makes regular plot of 1D data, with possibly multiple lines

        :param axObj: the axis object to be drawn on
        :param str xAxisLabel: x-axis label
        :param str yAxisLabel: y-axis label
        :param str titleLabel: title label
        """
        for lineNum in range(len(self.yVals)):
            axObj.plot(self.xVals, self.yVals[lineNum])

        axObj.set_xlabel(xAxisLabel)
        axObj.set_ylabel(yAxisLabel)
        if titleLabel != None:
            axObj.set_title(titleLabel)