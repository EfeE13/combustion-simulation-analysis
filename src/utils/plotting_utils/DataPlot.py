import matplotlib.pyplot as plt
from src.utils import utils
import time

class DataPlot:
    """
    Object containing some number of independent subplots, for ease of display (not all subplots need to be used before saving)

    Attributes
    rowCount: int
    colCount: int
    fig: plt-based object
    axes: plt-based objects
    freeSubplots: Starts off as something like [0, 1, 2, ..., 17] representing the positions where there are subplots that haven't been used
    """

    def __init__(self, rowCount: int, colCount: int) -> None:
        """An object for having greater than or equal to one plot per window

        :param int rowCount: number of rows
        :param int colCount: number of columns
        """

        self.rowCount = rowCount
        self.colCount = colCount
        self.fig, self.axes = plt.subplots(rowCount, colCount)
        self.freeSubplots = list(range(self.rowCount * self.colCount))
        if rowCount == 3 and colCount == 6:
            # Specifications give decent apperance for 3 rows, 6 columns:
            self.fig.set_size_inches(30.0, 11.0)
            self.fig.subplots_adjust(left = 0.03, right = 0.97, top = 0.9, bottom = 0.07, wspace = 0.9, hspace = 0.5)
        elif rowCount == 1 and colCount == 1:
            pass
        else:
            print("The figure size can be adjusted in the DataPlot.py module.")
    
    def getFigure(self):
        return self.fig

    def getFreeSubplot(self, pos:int = None):
        """Returns a subplot that's free

        :param int pos: None if you want it to just use the first free subplot location, int if you have a preferred location, defaults to None
        :return: plt-based object
        """
        if self.rowCount == 1 and self.colCount == 1:
            return self.axes
        else:
            if pos == None:
                assert len(self.freeSubplots) != 0
                pos = min(self.freeSubplots)
            assert pos in self.freeSubplots

            xCor = pos % self.colCount
            yCor = int((pos - xCor) / self.colCount)
            self.freeSubplots.remove(pos)
            return self.axes[yCor][xCor]
    
    def saveFigure(self, saveLocation: str, res: int):
        """Will override other figures in the same location with the same filename

        :param str saveLocation: location you want the figure saved to (include extension like .pdf or .png, if no extension it will save as both .pdf and .png)
        :param int res: represents the resolution (500 seems more than sufficient)
        """
        saveDir = utils.fixDirFormat(saveLocation[:saveLocation.rindex("/")])
        utils.makeDir(saveDir)
        
        print("- Saving to " + saveLocation)
        startTime = time.time()
        if "." in saveLocation:
            self.fig.savefig(saveLocation, dpi = res, bbox_inches = "tight")
        else:
            self.fig.savefig(saveLocation + ".png", dpi = res, bbox_inches = "tight")
            self.fig.savefig(saveLocation + ".pdf", dpi = res, bbox_inches = "tight")
        print("\t- Saved to " + saveLocation + " in", round(time.time() - startTime, 2), "sec.")
        plt.close()