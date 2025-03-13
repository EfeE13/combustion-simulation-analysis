import numpy as np
from scipy.optimize import fsolve

# STREAMS is a dictionary from types of streams to the grams of oxygen that would 
# be "leftover" if 1 gram of the stream underwent complete combustion. As
# expected, the "leftover" oxygen values are negative for the fuels toluene
# and n-heptane.
LEFTOVER_OX = {"toluene":-1/92.14*9*31.999, "air":1*0.23292, "n-heptane":-1/100.21*11*31.999}
TOL = 0.00000001

class Q2DFGeneralizer:
    def __init__(self, Z1_CFD:float, Z2_CFD:float, chi11_CFD:float, chi12_CFD:float, chi22_CFD:float, stream1:str, stream2:str, stream3:str) -> None:
        self.Z_CFD = [Z1_CFD, Z2_CFD, 1 - Z1_CFD - Z2_CFD]

        chi13_CFD = -chi11_CFD - chi12_CFD
        chi23_CFD = -chi22_CFD - chi12_CFD
        chi33_CFD = chi11_CFD + 2*chi12_CFD + chi22_CFD
        self.chi_CFD = [[chi11_CFD, chi12_CFD, chi13_CFD], 
                [chi12_CFD, chi22_CFD, chi23_CFD],
                [chi13_CFD, chi23_CFD, chi33_CFD]]

        self.streams = [stream1, stream2, stream3]
        self.leftover_ox_CFD = [LEFTOVER_OX[self.streams[0]], LEFTOVER_OX[self.streams[1]], LEFTOVER_OX[self.streams[2]]]

    def get_optimal_mapping(self):
        optimal_mapping_candidates = {}   # chi ratio: [i, j, x] where i and j are the indices of the compositions that will be labeled stream 1 and stream 2, respectively
        for i in [0, 1, 2]:
            j_options = [0, 1, 2]
            j_options.remove(i)
            for j in j_options:
                Z1 = self.Z_CFD[i]
                Z2 = self.Z_CFD[j]
                chi11 = self.chi_CFD[i][i]
                chi12 = self.chi_CFD[i][j]
                chi22 = self.chi_CFD[j][j]
                leftover_ox = [self.leftover_ox_CFD[i], self.leftover_ox_CFD[j], self.leftover_ox_CFD[0+1+2-i-j]]

                extrema = self.get_extrema_chi_ratio(Z1, Z2, chi11, chi12, chi22)
                valid_x_ranges = self.get_valid_x_ranges(Z1, Z2, leftover_ox)
                #print(">>For " + self.streams[i] + ", " + self.streams[j] + ", " + self.streams[0+1+2-i-j] + ", the valid ranges are:", valid_x_ranges)

                optimal_x_candidates = []
                for extremum in extrema:
                    for valid_x_range in valid_x_ranges:
                        if valid_x_range[0] < extremum and extremum < valid_x_range[1]:
                            optimal_x_candidates.append(extremum)
                            break
                for valid_x_range in valid_x_ranges:
                    optimal_x_candidates.append(valid_x_range[0])
                    optimal_x_candidates.append(valid_x_range[1])
                #print("^^^^^^^", optimal_x_candidates)
        
                for x_candidate in optimal_x_candidates:
                    chi_eta = self.get_chi_eta(Z1, Z2, chi11, chi12, chi22, x_candidate)
                    chi_xi = self.get_chi_xi(Z1, Z2, chi11, chi12, chi22, x_candidate)
                    chi_ratio = chi_eta / chi_xi

                    if not np.isnan(chi_ratio):
                        optimal_mapping_candidates[float(chi_ratio)] = [i, j, x_candidate]
        
        #print("*********", optimal_mapping_candidates, "and", min(list(optimal_mapping_candidates.keys())))
        
        optimal_mapping = optimal_mapping_candidates[min(list(optimal_mapping_candidates.keys()))]
        i, j, x = optimal_mapping[0], optimal_mapping[1], optimal_mapping[2]

        return self.get_outputs(i, j, x)

    def get_outputs(self, i:int, j:int, x:float):
        Z1 = self.Z_CFD[i]
        Z2 = self.Z_CFD[j]
        chi11 = self.chi_CFD[i][i]
        chi12 = self.chi_CFD[i][j]
        chi22 = self.chi_CFD[j][j]
        leftover_ox = [self.leftover_ox_CFD[i], self.leftover_ox_CFD[j], self.leftover_ox_CFD[0+1+2-i-j]]

        b = float(self.get_b(Z1, Z2, x))
        a = float(x * b)

        leftover_ox_left = a * leftover_ox[0] + (1-a) * leftover_ox[2]
        leftover_ox_right = b * leftover_ox[0] + (1-b) * leftover_ox[1]
        if leftover_ox_left >= -TOL and leftover_ox_right <= TOL:
            OMIX_shuffled = [a, 0, 1-a]
            FMIX_shuffled = [b, 1-b, 0]
        elif leftover_ox_left <= TOL and leftover_ox_right >= -TOL:
            OMIX_shuffled = [b, 1-b, 0]
            FMIX_shuffled = [a, 0, 1-a]
        else:
            raise Exception("Leftover oxygen is " + str(leftover_ox_left) + " and " + str(leftover_ox_right))
        
        OMIX = [-1, -1, -1]
        FMIX = [-1, -1, -1]
        OMIX[i] = OMIX_shuffled[0]
        FMIX[i] = FMIX_shuffled[0]
        OMIX[j] = OMIX_shuffled[1]
        FMIX[j] = FMIX_shuffled[1]
        OMIX[0+1+2-i-j] = OMIX_shuffled[2]
        FMIX[0+1+2-i-j] = FMIX_shuffled[2]

        chi = self.get_chi_xi(Z1, Z2, chi11, chi12, chi22, x)
        #print("chi ratio is", self.get_chi_ratio(Z1, Z2, chi11, chi12, chi22, x))

        x_prime_dict = {"12":[3,4], "13":[5,4], "21":[3,2], "23":[1,2], "31":[5,6], "32":[1,0]}
        x_prime_left = x_prime_dict[str(i + 1) + str(j + 1)][0]
        x_prime_right = x_prime_dict[str(i + 1) + str(j + 1)][1]
        x_prime = x_prime_left + x * (x_prime_right - x_prime_left)
        return chi, OMIX, FMIX, x_prime

    def get_valid_x_ranges(self, Z1:float, Z2:float, leftover_ox:list[float]) -> list[list[float]]:
        # stream1_stream3 variables denote constant-eta line intersecting stoichiometric stream1-stream3 mixture if it exists
        # # stream1_stream2 variables denote constant-eta line intersecting stoichiometric stream1-stream2 mixture if it exists
        stream1_stream3_a = leftover_ox[2] / (leftover_ox[2] - leftover_ox[0])
        stream1_stream2_b = leftover_ox[1] / (leftover_ox[1] - leftover_ox[0])

        stream1_stream3_b = 1 + Z2 / ((1-Z1-Z2)/(1-stream1_stream3_a) - 1)
        stream1_stream2_a = 1 - (1-Z1-Z2)/(Z2/(stream1_stream2_b-1) + 1)

        if 0 <= stream1_stream3_a and stream1_stream3_a <= 1 and 0 <= stream1_stream3_b and stream1_stream3_b <= 1:
            stream1_stream3_x = stream1_stream3_a / stream1_stream3_b
        else:
            stream1_stream3_x = -1

        if 0 <= stream1_stream2_a and stream1_stream2_a <= 1 and 0 <= stream1_stream2_b and stream1_stream2_b <= 1:
            stream1_stream2_x = stream1_stream2_a / stream1_stream2_b
        else:
            stream1_stream2_x = -1

        x_boundaries = sorted([elm for elm in [0, 1, stream1_stream3_x, stream1_stream2_x] if elm >= 0 and elm <= 1])
        x_boundaries.reverse()
        assert 0 in x_boundaries and 1 in x_boundaries
        #print("\tFor below, x_boundaries =", x_boundaries)
        
        x1_left_leftover_ox = leftover_ox[0] * Z1 + leftover_ox[2] * (1 - Z1)
        x1_right_leftover_ox = leftover_ox[0] * Z1 + leftover_ox[1] * (1 - Z1)
        x1_valid = (x1_left_leftover_ox >= 0 and x1_right_leftover_ox <= 0) or (x1_left_leftover_ox <= 0 and x1_right_leftover_ox >= 0)
        x_boundaries_valid = [x1_valid]
        for i in range(len(x_boundaries) - 2):
            x_boundaries_valid.append(not x_boundaries_valid[-1])

        valid_x_ranges = []
        for i in range(len(x_boundaries_valid)):
            if x_boundaries_valid[i]:
                valid_x_ranges.append([x_boundaries[i + 1], x_boundaries[i]])

        return valid_x_ranges

    def get_chi_ratio(self, Z1:float, Z2:float, chi11:float, chi12:float, chi22:float, x:float) -> float:
        return self.get_chi_eta(Z1, Z2, chi11, chi12, chi22, x) / self.get_chi_xi(Z1, Z2, chi11, chi12, chi22, x)
    
    def get_chi_eta(self, Z1:float, Z2:float, chi11:float, chi12:float, chi22:float, x:float) -> float:
        eta, xi = self.get_eta_xi(Z1, Z2, x)
        numer = (1-eta)**2*chi11 - 2*eta*(1-eta)*(1-x)*chi12 + (eta*(1-x))**2*chi22
        denom = ((1-eta)*(x+(1-x)*xi) + eta*xi*(1-x))**2
        return numer/denom

    def get_chi_xi(self, Z1:float, Z2:float, chi11:float, chi12:float, chi22:float, x:float) -> float:
        eta, xi = self.get_eta_xi(Z1, Z2, x)
        numer = xi**2*chi11 + 2*xi*(x+(1-x)*xi)*chi12 + (x+(1-x)*xi)**2*chi22
        denom = ((1-eta)*(x+(1-x)*xi) + eta*xi*(1-x))**2
        return numer/denom

    def get_b(self, Z1:float, Z2:float, x:float) -> float:
        Z3 = 1 - Z1 - Z2
        b = (Z1+Z2+x*(Z1+Z3) - np.sqrt((Z1+Z2+x*(Z1+Z3))**2 - 4*Z1*x)) / (2*x)
        return b

    def get_eta_xi(self, Z1:float, Z2:float, x:float):
        b = self.get_b(Z1, Z2, x)
        eta = b
        xi = Z2 / (1 - b)
        return eta, xi
    
    def get_deriv_eta_xi(self, Z1:float, Z2:float, x:float):
        Z3 = 1 - Z1 - Z2
        b = self.get_b(Z1, Z2, x)
        db = 1/(4*x**2) * ((Z1+Z3-0.5/np.sqrt((Z1+Z2+x*(Z1+Z3))**2-4*Z1*x)*(2*(Z1+Z2+x*(Z1+Z3))*(Z1+Z3)-4*Z1)) * 2*x
                           - 2 * (Z1+Z2+x*(Z1+Z3) - np.sqrt((Z1+Z2+x*(Z1+Z3))**2 - 4*Z1*x)))
        deta = db
        dxi = Z2/((1-b)**2) * db
        return deta, dxi

    def get_extrema_chi_ratio(self, Z1, Z2, chi11, chi12, chi22):
        def deriv(x):
            eta, xi = self.get_eta_xi(Z1, Z2, x)
            deta, dxi = self.get_deriv_eta_xi(Z1, Z2, x)
            deriv_val = (-2*(1-eta)*deta*chi11 + 2*deta*(2*eta-1)*(1-x)*chi12 + 2*eta*(1-eta)*chi12 + 2*eta*(1-x)*(deta*(1-x)-eta)*chi22) * (xi**2*chi11 + 2*xi*(x+(1-x)*xi)*chi12 + (x+(1-x)*xi)**2*chi22) - (2*xi*dxi*chi11 + 2*dxi*(x+(1-x)*xi)*chi12 + 2*xi*(1+(1-x)*dxi-xi)*chi12 + 2*(x+(1-x)*xi)*(1+(1-x)*dxi-xi)*chi22) * ((1-eta)**2*chi11 - 2*eta*(1-eta)*(1-x)*chi12 + (eta*(1-x))**2*chi22)
            return deriv_val
    
        return fsolve(deriv, 0.5)
    
    def __assert_local_extremum(self, Z1:float, Z2:float, chi11:float, chi12:float, chi22:float, x:float):
        PERTURBATION = 0.000001
        low_chi_ratio = self.get_chi_ratio(Z1, Z2, chi11, chi12, chi22, x - PERTURBATION)
        mid_chi_ratio = self.get_chi_ratio(Z1, Z2, chi11, chi12, chi22, x)
        hig_chi_ratio = self.get_chi_ratio(Z1, Z2, chi11, chi12, chi22, x + PERTURBATION)

        if (mid_chi_ratio - low_chi_ratio) * (mid_chi_ratio - hig_chi_ratio) > 0:
            return True
        else:
            print("NON-EXTREMUM DETECTED")
            return False