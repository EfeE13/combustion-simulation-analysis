from src.pdrs_analysis.PDRsRun import PDRsRun

class PDRsRunsND:
    """
    The composition of a number of PDRsRun objects

    Attributes
    PDRsRunDict: Dictionary with key:value pairs as PDRsOutputFileName:PDRsRunObject
    fileInputsDict: Dictionary with key:value pairs as PDRsOutputFileName:{"SDR":30.0, "FMX":0.1}
    """
    
    def __init__(self, outputFilenameList:list[str]) -> None:
        """Makes a basic PDRs run object that the PDRsRuns1D and PDRsRuns2D inherits from. n refers to the number of parameters vaired when running PDRs

        :param list[str] outputFilenameList: represent the PDRs output filenames. They should be something like SDR_30.0__FMX_0.1.Y (or other extension)
        """
        
        self.PDRsRunsDict = {}
        self.fileInputsDict = {}
        for outputFilename in outputFilenameList:
            self.PDRsRunsDict[outputFilename] = PDRsRun(outputFilename)
            
            self.fileInputsDict[outputFilename] = {}
            fileInputsList = outputFilename[outputFilename.rindex("/") + 1:outputFilename.rindex(".")].split("__")
            for fileInput in fileInputsList:
                fileInputSplit = fileInput.split("_")
                assert len(fileInputSplit) == 2
                nickname = fileInputSplit[0]; value = fileInputSplit[1]
                self.fileInputsDict[outputFilename][nickname] = float(value.replace(",", "."))