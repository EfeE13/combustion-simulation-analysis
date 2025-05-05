import subprocess
import os
from src.utils import utils
from src.pdrs_analysis.PDRsRun import PDRsRun

class PDRsRunnableNonCartesian:   # 1D maniZ-type manifold calculations
    def __init__(self, input_format:str, inputs, input_file_loc:str, output_file_loc:str, Z_val:float):
        """Generates the input file contents based on the format and input values

        :param str input_format: The input format with the input variables replaced by keywords like "OMIX1"
        :param _type_ inputs: Dictionary from the aforementioned keywords like "OMIX1" to the corresponding values
        :param str input_file_loc: The location of the input file (includes input file name and abs. path)
        :param str output_file_loc: The location of the output file (includes output file name and abs. path)
        :param float Z_val: The Z value to extract from the profile
        """
        self.input = input_format.strip()
        for keyword in inputs.keys():
            self.input = self.input.replace(keyword, str(inputs[keyword]))
        
        self.input_file_loc = utils.fixDirFormat(input_file_loc)
        self.output_file_loc = utils.fixDirFormat(output_file_loc)
        self.input = self.input.replace("OUTPUT_DIR", self.output_file_loc[:self.output_file_loc.rindex("/")])
        # Adds Output file name field to the input file, requiring a change to the PDRs src code
        self.input = self.input.replace("OUTPUT_FILE_NAME", self.output_file_loc[self.output_file_loc.rindex("/") + 1:])

        self.Z_val = Z_val
    
    def make_input(self):
        with open(self.input_file_loc, "w") as file:
            file.write(self.input)

    def run_pdrs(self):
        """Runs PDRs and returns the below dictionary

        :return: Dictionary from the output variables like Y_CO2 to the values at self.Z_val
        """
        error = False
        
        try:
            print("Running PDRs with", self.input_file_loc[self.input_file_loc.rindex("/") + 1:])
            subprocess.run("pdrs " + self.input_file_loc[self.input_file_loc.rindex("/") + 1:], timeout = 250, shell = True, cwd = self.input_file_loc[:self.input_file_loc.rindex("/")])
        except subprocess.TimeoutExpired:
            print("\tTimed out")
            error = True

        if not error:
            return self.get_output()
        else:
            return "NA"
    
    def get_output(self):
        if os.path.exists(self.output_file_loc):
            pdrs_run = PDRsRun(self.output_file_loc)
            pdrs_headers = pdrs_run.getHeaders()
            assert "Z" in pdrs_headers
            pdrs_headers.remove("Z")

            pdrs_data = {}
            for header in pdrs_headers:
                pdrs_data[header] = float(pdrs_run.interpolate("Z", self.Z_val, header))

            return pdrs_data
        else:
            return "NA"