class PDRsRun:
    """
    Given outputFileName, generates a dictionary from column headers to data as a list
    
    Attributes
    headerToDataDict: A dictionary like {"Z":[0.00, 0.023, ..., 1.00], "Y-CO":[0.00, 0.34789, ..., 0.00], ...}
    """

    def __init__(self, outputFilename:str) -> None:
        """Populates headerToDataDict with the data from the output file passed

        :param str outputFilename: PDRs output file (e.g., .Y file) the PDRsRun object will correspond to
        """

        def strToFloat(strNum: str) -> float:
            """For exponents with three digits like 2.819755818043E-105, PDRs may write 2.819755818043-105 without the E. This function tries to convert all string numbers, even if the formatting is like above, to the corresponding float.

            :param str strNum: string representing number
            :return float: the number itself
            """
            if ('-' in strNum[1:] or '+' in strNum) and (not ('e' in strNum.lower())):
                posMinusOrPlus = max(strNum[1:].find('-') + 1, strNum.find('+'))   # The 1's are because there could be a negative sign at very front of number
                decimalFormatStrNum = strNum[:posMinusOrPlus] + 'E' + strNum[posMinusOrPlus:]
            else:
                decimalFormatStrNum = strNum
            return float(decimalFormatStrNum)

        outputFile = open(outputFilename, "r")
        headersList = outputFile.readline().strip().split()
        assert len(headersList) > 1

        self.headerToDataDict = {}
        for header in headersList:
            self.headerToDataDict[header] = []

        for line in outputFile:
            if line.strip() != "":
                splitLineList = line.strip().split()
                assert len(headersList) == len(splitLineList)
                for i in range(len(headersList)):
                    self.headerToDataDict[headersList[i]].append(strToFloat(splitLineList[i]))

        outputFile.close()

    def getCol(self, columnHeader:str) -> list[float]:
        assert columnHeader in self.headerToDataDict.keys()
        return self.headerToDataDict[columnHeader]

    def interpolate(self, inputHeader:str, inputVal:float, outputHeader:str) -> float:
        """There's a list of data corresponding to inputHeader and for outputHeader.
        If inputVal is 30% of the way from one value in inputHeader list to another, return 30% of the way between the corresponding outputHeader list values
        Assumes inputHeader's data is everywhere non-decreasing. If inputVal < min(inputHeader list) or > max(outputHeader list), it just returns the 
        first or last element of outputHeader list.

        :param str inputHeader: the header for the variable with the value that will be used for interpolation
        :param float inputVal: the value that will be used for interpolation
        :param str outputHeader: the output header (the header representing what will actually be plotted)
        :return float: representing the interpolated value as explained above
        """
        inputCol = self.getCol(inputHeader)
        outputCol = self.getCol(outputHeader)
        if inputVal <= inputCol[0]:
            return outputCol[0]

        if inputVal >= inputCol[len(inputCol) - 1]:
            return outputCol[len(inputCol) - 1]
        
        for i in range(len(inputCol) - 1):
            if inputCol[i] <= inputVal and inputVal <= inputCol[i + 1]:
                return outputCol[i] + (outputCol[i + 1] - outputCol[i]) * (inputVal - inputCol[i]) / (inputCol[i + 1] - inputCol[i])
        
        raise Exception(str(inputVal) + " could not be interpolated in " + str(inputCol))