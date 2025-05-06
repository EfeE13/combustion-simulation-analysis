from src.nga_analysis.TimeStep import TimeStep
from src.nga_analysis.Q2DFGeneralizer import Q2DFGeneralizer
from src.utils import utils
from src.pdrs_analysis.PDRsRunnableNonCartesian import PDRsRunnableNonCartesian
from src.nga_analysis import constants

class OnTheFlyAndQ2DF2Runner:
    def __init__(self, timestep, stream1:str, stream2:str, stream3:str):
        """
        :param timestep: A TimeStep object
        :param str stream1: e.g., "toluene"
        :param str stream2: e.g., "air"
        :param str stream3: e.g., "n-heptane"
        """
        self.timestep = timestep
        self.stream1 = stream1
        self.stream2 = stream2
        self.stream3 = stream3
    
    def run_pdrs(self, input_folder:str, output_folder:str, log_file:str):
        input_folder = utils.fixDirFormat(input_folder)
        output_folder = utils.fixDirFormat(output_folder)

        on_fly_data_dict = {}
        q2df_data_dict = {}

        for i in range(self.timestep.getNumCells()):
            Z1 = self.timestep.getData("Z1", i)
            Z2 = self.timestep.getData("Z2", i)
            chi11 = self.timestep.getData("chi11", i)
            chi12 = self.timestep.getData("chi12", i)
            chi22 = self.timestep.getData("chi22", i)

            my_q2df_gen = Q2DFGeneralizer(Z1, Z2, chi11, chi12, chi22, self.stream1, self.stream2, self.stream3)
            OMIX, FMIX, chi, chi_eta, x_prime, eta, Z_opt, Z_stoic = my_q2df_gen.get_optimal_mapping()
            chi_ref = utils.er_func(0.5) * chi / utils.er_func(Z_opt)
            inputs = {"CHI_REF":chi_ref, "OMIX1":OMIX[0], "OMIX2":OMIX[1], "OMIX3":OMIX[2], "FMIX1":FMIX[0], "FMIX2":FMIX[1], "FMIX3":FMIX[2]}
            name = "OnFly_" + "0"*(4-len(str(i))) + str(i) + "_"
            my_runnable = PDRsRunnableNonCartesian(constants.MANI_OF_INPUT, inputs, input_folder + "/" +  name + "Input", output_folder + "/" + name + "Output.Y", Z_opt)
            pdrs_run_data_on_fly = my_runnable.get_output()

            Z_opt = 1 - Z2
            chi = chi22
            chi_ref = utils.er_func(0.5) * chi / utils.er_func(Z_opt)
            Z3 = 1-Z1-Z2
            FMIX = [Z1/(Z1+Z3), 0, Z3/(Z1+Z3)]
            inputs = {"CHI_REF":chi_ref, "OMIX1":"0.0", "OMIX2":"1.0", "OMIX3":"0.0", "FMIX1":FMIX[0], "FMIX2":FMIX[1], "FMIX3":FMIX[2]}
            name = "Q2DF2_" + "0"*(4-len(str(i))) + str(i) + "_"
            my_runnable = PDRsRunnableNonCartesian(constants.MANI_OF_INPUT, inputs, input_folder + "/" +  name + "Input", output_folder + "/" + name + "Output.Y", Z_opt)
            pdrs_run_data_q2df2 = my_runnable.get_output()
    
            add_dicts(on_fly_data_dict, pdrs_run_data_on_fly, i)
            add_dicts(q2df_data_dict, pdrs_run_data_q2df2, i)
            if len(list(on_fly_data_dict.keys())) != 0:
                assert len(on_fly_data_dict[list(on_fly_data_dict.keys())[0]]) == i + 1
            if len(list(q2df_data_dict.keys())) != 0:
                assert len(q2df_data_dict[list(q2df_data_dict.keys())[0]]) == i + 1

            # with open(log_file, "a") as file:
            #     file.write("\npdrs_run_data_on_the_fly[" + str(i) + "] = " + str(pdrs_run_data_on_fly) + "\npdrs_run_data_q2df2[" + str(i) + "] = " + str(pdrs_run_data_q2df2) + "\n# " + "-" * 50)
        
        # with open(log_file, "a") as file:
        #     file.write("\nTWO ENTIRE DICTIONARIES BELOW (On The Fly, then Q2DF2) " + "-"*30 + "\n" + str(on_fly_data_dict) + "\n" + str(q2df_data_dict))
        # print(on_fly_data_dict, q2df_data_dict)
        return on_fly_data_dict, q2df_data_dict

def add_dicts(dict1, dict2, i:int):
    # Adds the data (one "row" of data) from dict2 to dict1
    if dict2 != "NA":
        for key in dict2.keys():
            if key not in dict1.keys():
                dict1[key] = ["NA"] * i
            dict1[key].append(dict2[key])
    else:
        for key in dict1.keys():
            dict1[key].append("NA")