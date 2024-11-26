from src.utils import utils
import subprocess
import itertools
import time

class PDRsRunnable:
    def __init__(self, inputFileSpecs:str, inputFileNames:list[str] = None, outputFileNames:list[str] = None) -> None:
        """Can run a Cartesian product of inputs (or just one)

        :param str inputFileSpecs: the input file with values (commas allow for Cartesian product of inputs to be run)
        :param list[str] inputFileNames: the input filenames, defaults to None
        :param list[str] outputFileNames: the output filenames, defaults to None
        """
        newInputFileSpecs = ""
        commentedLines = ""
        # Preprocessing to remove commented lines:
        for line in inputFileSpecs.split("\n"):
            if line.strip() != "":
                if line.strip()[0] == "!":
                    commentedLines += line + "\n"
                else:
                    newInputFileSpecs += line + "\n"
        
        inputFileSpecs = newInputFileSpecs.strip()
        commentedLines = commentedLines.strip()

        paddedInputFileSpecs = "\n" + inputFileSpecs.strip() + "\n"
        splitInputFileSpecs = paddedInputFileSpecs.split(":")
        
        inputDescs = []
        inputNicknames = []
        inputVals = []
        self.numPDRsRuns = 1
        dimPDRsRuns = 0
        for i in range(0, len(splitInputFileSpecs) - 1):
            assert "\n" in splitInputFileSpecs[i]
            inputDescAndNickname = splitInputFileSpecs[i][splitInputFileSpecs[i].rindex("\n") + 1:].strip()
            inputValsSublist = splitInputFileSpecs[i + 1][:splitInputFileSpecs[i + 1].rindex("\n")].strip().split("|")

            if "{" in inputDescAndNickname:
                assert "}" in inputDescAndNickname
                inputDesc = inputDescAndNickname[:inputDescAndNickname.index("{")].strip()
                inputNickname = inputDescAndNickname[inputDescAndNickname.index("{") + 1:inputDescAndNickname.index("}")].strip()
            else:
                assert "}" not in inputDescAndNickname
                inputDesc = inputDescAndNickname
                inputNickname = ""
            
            inputDescs.append(inputDesc)
            inputNicknames.append(inputNickname)
            inputVals.append(inputValsSublist)
            if len(inputValsSublist) != 1:
                self.numPDRsRuns *= len(inputValsSublist)
                dimPDRsRuns += 1
        
        cartesianProductInputVals = list(itertools.product(*inputVals))
        assert len(cartesianProductInputVals) == self.numPDRsRuns

        self.inputFileContents = []
        if inputFileNames != None:
            assert len(inputFileNames) == self.numPDRsRuns
            self.inputFileNames = inputFileNames
        else:
            self.inputFileNames = []
        if outputFileNames != None:
            assert len(outputFileNames) == self.numPDRsRuns
            self.outputFileNames = outputFileNames
        else:
            self.outputFileNames = []
        
        for run in range(self.numPDRsRuns):
            inputFileContent = commentedLines + "\n"
            fileName = ""
            for inputNum in range(len(inputDescs)):
                inputVal = cartesianProductInputVals[run][inputNum]
                inputFileContent += inputDescs[inputNum] + " : " + inputVal + "\n"
                if inputNicknames[inputNum] != "":
                    fileName += inputNicknames[inputNum] + "_" + f"{float(inputVal):.5g}".replace(".", ",") + "__"
        
            self.inputFileContents.append(inputFileContent.strip())
            if inputFileNames == None:
                assert len(fileName) > 0
                self.inputFileNames.append(fileName[:-2] + ".multicomponent")
            if outputFileNames == None:
                assert len(fileName) > 0
                self.outputFileNames.append(fileName[:-2] + ".Y")
    
    def runPDRs(self, inputFolderFromStart:str, outputFolderFromStart:str, changeOutputNames:bool = True) -> list[str]:
        """Runs PDRs, changing output filenames for convenience

        :param str inputFolderFromStart: input folder name
        :param str outputFolderFromStart: output folder name
        :param bool changeOutputNames: whether or not to change output names, defaults to True
        :return list[str]: A list of output file locations that can then be used to construct another one of this package's objects for PDRs analysis
        """
        inputFolderFromStart = utils.fixDirFormat(inputFolderFromStart)
        outputFolderFromStart = utils.fixDirFormat(outputFolderFromStart)
        utils.makeDir(inputFolderFromStart)
        utils.makeDir(outputFolderFromStart)

        outputFileLocs = []
        for i in range(self.numPDRsRuns):
            inputFileLoc = inputFolderFromStart + "/" + self.inputFileNames[i]
            outputFileLoc = outputFolderFromStart + "/" + self.outputFileNames[i]
            print("- " + str(i + 1) + ". Running PDRs on", inputFileLoc, "(run", i + 1, "of", str(self.numPDRsRuns) + ")")
            if subprocess.run("du " + outputFileLoc,  shell = True, capture_output = True, text = True).stderr != "":
                startTime = time.time()
                inputFile = open(inputFileLoc, "w")
                inputFile.write(self.inputFileContents[i])
                inputFile.close()

                pdrsCommand = "pdrs " + inputFileLoc.replace("/home/efeeroz/Documents/CombustionModelAnalysis/", "") # The .replace is a temporary fix
                #subprocess.run(pdrsCommand, shell = True) # Lots of output/verbose
                utils.run(pdrsCommand) # Output suppressed
                if changeOutputNames:
                    utils.run("mv " + outputFolderFromStart + "/" + utils.run("ls -Art " + outputFolderFromStart + " | tail -n 1").strip() + " " + outputFileLoc)
                    finalOutputFileLoc = outputFileLoc
                else:
                    finalOutputFileLoc = utils.run("ls -Art " + outputFolderFromStart + " | tail -n 1").strip()
                outputFileLocs.append(finalOutputFileLoc)
                print("\t- Finished in", round(time.time() - startTime, 2), "sec., output at", finalOutputFileLoc)
            else:
                outputFileLocs.append(outputFileLoc)
                print("\t- Already ran previously, output at", outputFileLoc)
        
        return outputFileLocs