import numpy as np
from scipy.optimize import fsolve
import textwrap

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
        #print("------------------------"*2)
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

                #print(">>>>>>>EXTREMA", i+1, j+1, extrema)

                optimal_x_candidates = []
                for extremum in extrema:
                    for valid_x_range in valid_x_ranges:
                        if valid_x_range[0] < extremum and extremum < valid_x_range[1]:
                            optimal_x_candidates.append(extremum)
                            break
                for valid_x_range in valid_x_ranges:
                    optimal_x_candidates.append(valid_x_range[0])
                    optimal_x_candidates.append(valid_x_range[1])
                #print("Candidates", i+1, j+1, optimal_x_candidates)
        
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
            complement_xi = False
        elif leftover_ox_left <= TOL and leftover_ox_right >= -TOL:
            OMIX_shuffled = [b, 1-b, 0]
            FMIX_shuffled = [a, 0, 1-a]
            complement_xi = True
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

        chi_eta = self.get_chi_eta(Z1, Z2, chi11, chi12, chi22, x)
        eta, Z_opt = self.get_eta_xi(Z1, Z2, x)
        Z_stoic = leftover_ox_left / (leftover_ox_left - leftover_ox_right)
        if complement_xi:
            Z_opt = 1 - Z_opt
            Z_stoic = 1 - Z_stoic

        return OMIX, FMIX, chi, chi_eta, x_prime, eta, Z_opt, Z_stoic

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
        assert x_boundaries[0] == 1 and x_boundaries[-1] == 0
        for i in range(len(x_boundaries) - 1):
            assert x_boundaries[i] >= x_boundaries[i + 1]
        #print("\tFor below, x_boundaries =", x_boundaries)
        
        x1_left_leftover_ox = leftover_ox[0] * Z1 + leftover_ox[2] * (1 - Z1)
        x1_right_leftover_ox = leftover_ox[0] * Z1 + leftover_ox[1] * (1 - Z1)
        x1_valid = (x1_left_leftover_ox >= 0 and x1_right_leftover_ox <= 0) or (x1_left_leftover_ox <= 0 and x1_right_leftover_ox >= 0)
        x_boundaries_valid = [x1_valid]
        for i in range(len(x_boundaries) - 2):
            x_boundaries_valid.append(not x_boundaries_valid[-1])

        #print("LLLL", x1_left_leftover_ox, x1_right_leftover_ox, x1_valid)

        valid_x_ranges = []
        for i in range(len(x_boundaries_valid)):
            if x_boundaries_valid[i]:
                valid_x_ranges.append([x_boundaries[i + 1], x_boundaries[i]])

        #print("RUNGES", valid_x_ranges)
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

    def get_deriv_chi_ratio(self, Z1:float, Z2:float, chi11:float, chi12:float, chi22:float, x:float) -> float:
        eta, xi = self.get_eta_xi(Z1, Z2, x)
        deta, dxi = self.get_deriv_eta_xi(Z1, Z2, x)
        deriv_val = (-2*(1-eta)*deta*chi11 + 2*deta*(2*eta-1)*(1-x)*chi12 + 2*eta*(1-eta)*chi12 + 2*eta*(1-x)*(deta*(1-x)-eta)*chi22) * (xi**2*chi11 + 2*xi*(x+(1-x)*xi)*chi12 + (x+(1-x)*xi)**2*chi22) - (2*xi*dxi*chi11 + 2*dxi*(x+(1-x)*xi)*chi12 + 2*xi*(1+(1-x)*dxi-xi)*chi12 + 2*(x+(1-x)*xi)*(1+(1-x)*dxi-xi)*chi22) * ((1-eta)**2*chi11 - 2*eta*(1-eta)*(1-x)*chi12 + (eta*(1-x))**2*chi22)
        return deriv_val

    def get_extrema_chi_ratio(self, Z1:float, Z2:float, chi11:float, chi12:float, chi22:float):
        GRID = 1000
        x_list = np.linspace(0.000001, 1-0.000001, GRID)
        deriv_chi_ratio_list = [self.get_deriv_chi_ratio(Z1, Z2, chi11, chi12, chi22, x) for x in x_list]
        deriv_product_list = [deriv_chi_ratio_list[i] * deriv_chi_ratio_list[i + 1] for i in range(len(deriv_chi_ratio_list) - 1)]

        extrema = []
        for i, product in enumerate(deriv_product_list):
            if product <= 0:
                extrema.append((x_list[i] + x_list[i + 1]) / 2)

        return extrema

    """
    def get_extrema_chi_ratio2(self, Z1, Z2, chi11, chi12, chi22):
        def deriv(x):
            eta, xi = self.get_eta_xi(Z1, Z2, x)
            deta, dxi = self.get_deriv_eta_xi(Z1, Z2, x)
            deriv_val = (-2*(1-eta)*deta*chi11 + 2*deta*(2*eta-1)*(1-x)*chi12 + 2*eta*(1-eta)*chi12 + 2*eta*(1-x)*(deta*(1-x)-eta)*chi22) * (xi**2*chi11 + 2*xi*(x+(1-x)*xi)*chi12 + (x+(1-x)*xi)**2*chi22) - (2*xi*dxi*chi11 + 2*dxi*(x+(1-x)*xi)*chi12 + 2*xi*(1+(1-x)*dxi-xi)*chi12 + 2*(x+(1-x)*xi)*(1+(1-x)*dxi-xi)*chi22) * ((1-eta)**2*chi11 - 2*eta*(1-eta)*(1-x)*chi12 + (eta*(1-x))**2*chi22)
            return deriv_val
    
        return fsolve(deriv, 0.5)
    """
    
    def assert_local_extremum(self, Z1:float, Z2:float, chi11:float, chi12:float, chi22:float, x:float):
        PERTURBATION = 0.01
        low_chi_ratio = self.get_chi_ratio(Z1, Z2, chi11, chi12, chi22, max(x - PERTURBATION, 0))
        mid_chi_ratio = self.get_chi_ratio(Z1, Z2, chi11, chi12, chi22, x)
        hig_chi_ratio = self.get_chi_ratio(Z1, Z2, chi11, chi12, chi22, min(x + PERTURBATION, 1))

        if (mid_chi_ratio - low_chi_ratio) * (mid_chi_ratio - hig_chi_ratio) > 0:
            return True
        else:
            print("NON-EXTREMUM DETECTED")
            return False

def unit_test_fortran():
    NUM_TESTS = 148
    list_inputs = []

    Z1_CFD_list = []; Z2_CFD_list = []; chi11_CFD_list = []; chi12_CFD_list = []; chi22_CFD_list = []
    OMIX_truth = []; FMIX_truth = []; chi_xi_truth = []; chi_eta_truth = []; x_prime_truth = []; eta_opt_truth = []; Z_opt_truth = []; Z_stoic_truth = []

    for i in range(NUM_TESTS):
        Z1, Z2, chi11, chi12, chi22 = get_rand()
        list_inputs.append([Z1, Z2, chi11, chi12, chi22])

    list_inputs = [[0.2199014833190605, 0.3705680937512448, 0.37448988013761686, np.float64(-0.5320037866604121), 0.9767321738539616], [0.33724889867020524, 0.34392899696950824, 1.2080334384442448, np.float64(0.23563586287195093), 0.2640920929278583], [0.4986614756676147, 0.32705817515950003, 0.5504055865908986, np.float64(0.2945280199871081), 0.2299821167604854], [0.5473975725738289, 0.31998335176127896, 1.0925378386908877, np.float64(0.478713175030655), 1.0120684731457443], [0.6518446856741693, 0.046115196195361974, 1.0989881299301363, np.float64(-0.6923906372141617), 1.708256001518653], [0.10011579482263833, 0.4497745859137792, 1.9749569160344527, np.float64(0.2676786392517665), 0.7938850940011903], [0.2540913890650932, 0.4504820151801578, 0.5024357055548725, np.float64(-0.22443454591710185), 0.42306362958535493], [0.22912771734103338, 0.31941769907896295, 0.6823016980747854, np.float64(0.24440230493603893), 0.23568383253036562], [0.5532287089156661, 0.0814056258180373, 0.0273243485975192, np.float64(-0.09347096424153109), 1.5960817452382452], [0.42617771347983435, 0.18621384531728674, 1.8482880046591774, np.float64(0.8661307950603732), 0.7069485633681638], [0.05017197579677181, 0.4739657505431522, 0.8464867171963026, np.float64(0.6624171190445489), 0.5942785990647463], [0.39596483660459614, 0.3064783260060796, 0.729763741342554, np.float64(-0.04887573540833268), 0.09409702609616755], [0.46850268483227786, 0.38796853211730536, 0.32913223002243175, np.float64(-0.057337822486441814), 0.3634071367874687], [0.21152852794938165, 0.24532477722231166, 1.5601136768866235, np.float64(0.30925022269070035), 1.9145055695548616], [0.31978147817863894, 0.24313528433168716, 0.3967115478546368, np.float64(-0.1604091744958509), 1.9394571621156174], [0.08521825767419386, 0.6421625137595265, 0.6252017135087902, np.float64(0.4120984945549112), 0.40815156154820187], [0.05860632441878657, 0.3626880748917589, 0.46485487718235663, np.float64(-0.28127129242758414), 0.2830523779120533], [0.26087143946165775, 0.2717566926176325, 1.5170183562340402, np.float64(1.0458068601154276), 1.3078704865042496], [0.6852275764235453, 0.048905252402952896, 0.8723349175607615, np.float64(0.35438080247724835), 0.315940641705168], [0.3353462937557347, 0.4642142604341024, 1.8883864428777344, np.float64(-0.8725486001797165), 0.6385183027054233], [0.5603224618444705, 0.05617292409388961, 0.1672949127537655, np.float64(0.20552544948667334), 0.5427153842709038], [0.28192800758686465, 0.13177021107482748, 0.834374688562945, np.float64(-0.994208845838252), 1.357977209296319], [0.15240056342074515, 0.10494491750435105, 1.410442075366921, np.float64(-1.3753994835604468), 1.540279139499389], [0.38160982109391983, 0.3870044922330624, 0.25396961381471295, np.float64(-0.13042237092158096), 0.1345394993070026], [0.15648022727515504, 0.46322163807163014, 1.3932056861325763, np.float64(0.007107620815502758), 0.5279132176174415], [0.4127767238789505, 0.04916554287202971, 0.40150820183539215, np.float64(-0.1224633564874176), 0.18249120545801456], [0.19185443818191086, 0.4114918693296534, 1.068154220090624, np.float64(0.440479577838399), 0.6106636155704146], [0.7227385333773537, 0.1398867526165806, 1.7562223826988685, np.float64(1.120208378032193), 0.9904161597982089], [0.08481392176514818, 0.6609007393246171, 1.2726279788886494, np.float64(-0.80663702806795), 1.9810419461654414], [0.3557030916482899, 0.36441405001385996, 0.35558205073943383, np.float64(0.4387220050661861), 0.9163101187783578], [0.44986806407635177, 0.018435103788960176, 1.5897733657280848, np.float64(0.42374144689981924), 0.13140705966415278], [0.5291482372882572, 0.28617649604162654, 0.8628043291847629, np.float64(-0.30303341667467576), 0.2724687795662817], [0.18575446369421877, 0.4142045497766184, 0.09657572886132093, np.float64(-0.2509403468688064), 0.830356065406284], [0.410379911989409, 0.3868160343638349, 0.06437340550152304, np.float64(-0.13207927714366266), 0.8950795907034639], [0.4262601527358208, 0.06911246898010968, 0.9155817521092293, np.float64(-0.1302705180215593), 1.8873531752438442], [0.10301395132197028, 0.6207405109209677, 0.7674240652837501, np.float64(1.0086455296692498), 1.3497359251170455], [0.26732406560547023, 0.2826328581986067, 1.8635666321390718, np.float64(0.06130944592639698), 1.3313246913747685], [0.015577285366665146, 0.5664602650861379, 0.5247805688594531, np.float64(0.26853266846917656), 0.9906787344183767], [0.440566052479085, 0.2819561026635625, 1.3702097619538702, np.float64(-0.8167952437943864), 1.2530021881821847], [0.2182344295966475, 0.5350636983204851, 1.6981280138795412, np.float64(-1.1224999668605493), 1.2753393253474221], [0.1402416982133526, 0.7730029852720429, 0.6404995369312549, np.float64(0.14625478817439308), 0.7930912536617669], [0.5785432965280617, 0.26556542417113727, 0.701071945978126, np.float64(-0.1522014148272574), 0.31046533720711], [0.3403689992201045, 0.36132427507712656, 1.5119059331790106, np.float64(-1.194010427956042), 1.8283370624450053], [0.6238387205390453, 0.2685228765249786, 0.5275943002401089, np.float64(-0.5450983601742062), 0.625673633237817], [0.3248215986019519, 0.41756888548365234, 0.7180401350831824, np.float64(1.1381741716023244), 1.949154596462041], [0.23609184978666217, 0.38940371020191705, 1.1904368800158878, np.float64(0.8185083306916844), 0.5932859052820307], [0.5447778257083787, 0.4517058857701718, 1.7647003400074834, np.float64(0.9029168287906786), 0.7861454939316261], [0.37332424851416274, 0.2925467696602508, 0.2925069867300407, np.float64(0.3205044306725915), 1.02306371341613], [0.3790366126450805, 0.5158959100072742, 1.4652238128030857, np.float64(-1.1663052832946548), 1.7744547268038224], [0.2984766060199765, 0.2064114849399907, 0.4076734744295445, np.float64(-0.321513349589935), 0.9727149037656684], [0.09141904812902701, 0.784052168928812, 1.2258662323573903, np.float64(-0.2311043871134415), 0.7919226533241905], [0.5794287269844878, 0.1856240320522414, 0.7350400795391094, np.float64(-0.7130532011711956), 1.5848924101211843], [0.4702399633396789, 0.37477895221933083, 0.07141439700788532, np.float64(-0.05214470234861428), 0.4079812716625326], [0.38082351796535896, 0.34126810812432, 0.31286561789027023, np.float64(-0.3090089215338304), 0.969344589596074], [0.366452490710751, 0.27506489629970954, 1.9125013592185598, np.float64(-0.12038731386336754), 0.912155691675336], [0.41042525537703367, 0.4338615394967612, 1.4058174561693906, np.float64(1.422446084034911), 1.9565635873338854], [0.32471802425049406, 0.3334863868284931, 0.44441661252136155, np.float64(0.34628940193636926), 0.6810648357324167], [0.39981158709912334, 0.3989873628897666, 0.5222624226660373, np.float64(0.35551191346501765), 0.9247913432983428], [0.6527941297948753, 0.17275641071008677, 1.3951101074366286, np.float64(-0.8087339077372864), 0.9403531793858524], [0.4043399535575942, 0.18672046802426834, 1.6447689543571302, np.float64(0.9170578212528193), 1.1461975820128807], [0.3737644840093439, 0.17866973403018435, 1.4554268516222, np.float64(-0.6663547469505713), 0.4940694348629191], [0.33539588601218445, 0.24402278507822608, 0.3825984765893611, np.float64(-0.21719495544653233), 0.6539174775158023], [0.4693843619898315, 0.23494653535259788, 1.7627177830302416, np.float64(-0.03230668135477632), 0.9955305789590359], [0.3417361669300392, 0.3088041275444134, 0.4217835541174666, np.float64(0.14523454601698482), 1.6012321147979975], [0.42444859120006884, 0.2463194502096872, 1.053895106439828, np.float64(0.03437743783244551), 1.034700391688032], [0.5936092401607916, 0.04882008621365231, 1.7726365914389606, np.float64(-0.08119054340141996), 1.6320217939066162], [0.4046750123424271, 0.43651196458612196, 1.3237436348992115, np.float64(0.3715922386193339), 0.16019390392001243], [0.3995551672280777, 0.0895003733545358, 1.3051290571271017, np.float64(0.3946647601593958), 1.1383713729218887], [0.41840635583207475, 0.37416858022650495, 0.8602153714801837, np.float64(0.6561867102068776), 0.6519252895384495], [0.3207757646093318, 0.47002297238415236, 0.5645751340876433, np.float64(0.6750836374977776), 1.8830039833941365], [0.4645072359392396, 0.10479788953145505, 0.20816495254153233, np.float64(0.21321738972189763), 0.44714725551420287], [0.23103225424586576, 0.4031677857690975, 1.9213858921410247, np.float64(0.05924700990519216), 0.48593463333132636], [0.14031228757077374, 0.28612020597500043, 1.0981524961503244, np.float64(0.39046700208212015), 1.3552940652182421], [0.5679987279397191, 0.008492308474438865, 0.36669991435611315, np.float64(0.45195126837380206), 1.1858084893097764], [0.042928147514109914, 0.5518617637366549, 1.6823613653769602, np.float64(0.9197628809966887), 1.3994261643728207], [0.08382232557046158, 0.4065115450893835, 0.07059586738277801, np.float64(0.02583081965901865), 0.5763485789042675], [0.6448766846548963, 0.1611668592886759, 0.956955922331739, np.float64(0.6323430068550735), 1.196715601680969], [0.28958802854562904, 0.3201205215964496, 0.30525397377881736, np.float64(0.01995865536266095), 1.2022363144875858], [0.2050861958148722, 0.6251392296144408, 1.5117377454715273, np.float64(0.6077609082730555), 1.5734165755739793], [0.10083073983905352, 0.2016286212834643, 1.2623391616877935, np.float64(0.13844779773203397), 1.5896383588150451], [0.28072608914963615, 0.35270984329647836, 1.3278824921452705, np.float64(0.5649853464925103), 0.28582500348348594], [0.3704117950187408, 0.2900585557345272, 1.1992303352251623, np.float64(-0.847124480048162), 0.6001123225065945], [0.27863002074775045, 0.462715639482829, 0.8635962753437723, np.float64(-0.012177059803286427), 1.363753615741378], [0.586435512639438, 0.13559452113671347, 1.793024279328282, np.float64(-1.599243963243535), 1.4498873280402937], [0.8626431115515417, 0.0895590273985805, 0.712881958640577, np.float64(-0.18156026271803466), 0.19231815963567866], [0.207390174164105, 0.5911628581781774, 0.24396886160691245, np.float64(-0.34582851903870293), 1.7997687608046742], [0.3232471329590064, 0.29806084528983595, 0.6068370527734206, np.float64(0.4947595669752256), 0.4264648136643914], [0.2953569342936011, 0.3659500303594625, 0.6047549920909814, np.float64(0.4155260274056104), 1.0412127869525951], [0.6167907921731448, 0.18463588635622183, 1.5294785365014805, np.float64(0.460549990730005), 1.7444439879385827], [0.517376967819739, 0.07541279585211762, 1.9409091944775263, np.float64(0.9059306098166788), 1.6251203183255543], [0.5638173904655239, 0.37395397940961217, 0.89826244736222, np.float64(0.7622571582395441), 1.8252587949813641], [0.13958116275434468, 0.42719626898798646, 1.6888735150082952, np.float64(-0.18437593873550218), 0.17712428563628158], [0.18671133826199568, 0.4875616991428157, 1.7108055505649407, np.float64(-0.8453218273531848), 0.5539123303798772], [0.782938069412775, 0.046855476157732195, 0.589531422394793, np.float64(0.02187584573339664), 0.023696424513304004], [0.3925356832399642, 0.17739139049486743, 1.8716786147932314, np.float64(0.09964869772989343), 0.06183607949490222], [0.47031330206222344, 0.4332643553915542, 0.8664413120458836, np.float64(0.5524246961495529), 1.4989353929430698], [0.8011651267471145, 0.08214465617291684, 1.1462600779351972, np.float64(-0.007132188757407287), 0.8247493447365684], [0.3060994102202625, 0.48798773495581793, 1.2102165492924912, np.float64(0.7891452591607631), 1.3653621000088398], [0.00399882969264699, 0.9276180424799187, 1.991649132992543, np.float64(-0.5876832475550082), 1.748453430610656], [0.5222969889557657, 0.46793765336655935, 1.9282751744859152, np.float64(0.38801546874723725), 1.6942489525958737], [0.10830361514460254, 0.46948012081118123, 1.5336567625132957, np.float64(-1.0152104995956592), 1.7626205955701821], [0.5023802478647574, 0.13873536490036467, 1.7614800730224924, np.float64(0.737265270268444), 0.36625544985602065], [0.3042369767885765, 0.4236691547151846, 0.16665827603304684, np.float64(0.03390790782393477), 0.21888729954507213], [0.3466747177091465, 0.06622778204991378, 1.9874975043479999, np.float64(0.06893043805817389), 0.7312330492739738], [0.35820610982018725, 0.3426569173980422, 0.4639688996857916, np.float64(-0.19498521573118963), 0.13415107971212503], [0.4936346904206096, 0.283108888342331, 0.21090616711623134, np.float64(-0.010588859924052202), 0.17347336508952105], [0.3809677511562454, 0.3939151715604325, 0.4029628949014801, np.float64(0.3688654218725612), 1.2443335777766495], [0.33348140630186907, 0.09586823519507531, 1.939942672444102, np.float64(1.3679766934908717), 1.6994405241074764], [0.332920218958037, 0.39693912354456784, 1.6330542930873646, np.float64(-0.02742112031635921), 0.014359867080624422], [0.7100720807137153, 0.12999319943765003, 0.32546523127861837, np.float64(0.5226459226822942), 1.4711411695479462], [0.14749782530644379, 0.10066456031830967, 1.2635767067696342, np.float64(-0.0004606786115146644), 1.3561534152412895], [0.6786469761448588, 0.2814714237195239, 1.5662429024762834, np.float64(-0.009814961749722073), 1.0688090294016515], [0.41175940237229086, 0.5180677339749454, 0.05282729445019818, np.float64(-0.10747335304396946), 0.7389032278557499], [0.37900417894494526, 0.25726316947696887, 1.2975436691021063, np.float64(-0.019922285233516468), 0.10591065375928133], [0.21329443289138078, 0.3824543742217693, 0.943665859368576, np.float64(-0.3592477321150176), 0.802176261652614], [0.39248584947045156, 0.42488768478570216, 1.7367828948941442, np.float64(0.950191885475488), 0.5903916287151363], [0.39643645494672664, 0.20740996824395966, 0.878955812488238, np.float64(0.08117813545016149), 0.009192193776444535], [0.5613989593303148, 0.2339370305712483, 0.7918423109744743, np.float64(-0.4307292947884317), 1.647965254558636], [0.2125595542874605, 0.38607983379980637, 1.6322182974200652, np.float64(0.16681136863838764), 1.7556462747029364], [0.48534364232014165, 0.41240799076257556, 0.19944507611308904, np.float64(0.520543163143804), 1.7654989894979967], [0.07030791070963247, 0.4345804727585336, 1.988104889740077, np.float64(1.0971597652873715), 1.4307919175052428], [0.24845062668257015, 0.4543281382020148, 1.9224147089722698, np.float64(0.05252284659858297), 0.6228890876233206], [0.5388951988688745, 0.4445272766526614, 1.4315311807753552, np.float64(-1.2146986661836985), 1.6822588168068204], [0.36928628524833174, 0.2074102748829358, 0.42573785209207826, np.float64(-0.39760285257397565), 1.0712154765271242], [0.42107987829336113, 0.284625854669468, 0.7158613273175176, np.float64(0.06652641864863493), 0.020938521588204972], [0.13772788309397188, 0.3007142481803916, 0.21790630516364184, np.float64(-0.248299582484595), 0.5255475434993422], [0.5818216418492889, 0.37814076642593225, 0.821630496987328, np.float64(0.1794433990425347), 0.9875342742838447], [0.5031793325623714, 0.10710190737014635, 0.9291607418439691, np.float64(-0.8725307734560586), 1.6529814530496385], [0.16871974416159205, 0.5117412450144577, 1.4396935838593845, np.float64(-0.22926706174969946), 0.47535889442163093], [0.5877516691573587, 0.35927797873800077, 1.2930480756487517, np.float64(-0.4556978592929579), 1.6221954269794843], [0.7595470831297896, 0.10434478839910205, 1.1252820095911835, np.float64(0.6466101743347573), 0.8032661205650753], [0.3320904697213568, 0.3981091361962083, 0.045619568993587656, np.float64(-0.17109814321494915), 0.8929282339645246], [0.1352489077352564, 0.28573847311088313, 0.8370699005664852, np.float64(-0.48400981477323174), 0.9072547064615017], [0.3322297593342785, 0.4703936787525612, 1.5899711517601285, np.float64(-0.26857577879559064), 0.0904612789554422], [0.4076134570189294, 0.5422076537450176, 1.2687672851711982, np.float64(-0.023294509931201413), 1.8678838297615477], [0.23770575038778705, 0.34866512803915783, 1.0836998083729703, np.float64(-0.25264020026066714), 0.8344595484872188], [0.31066432093050494, 0.2847647685621477, 0.2689098074149967, np.float64(7.978640186638575e-05), 0.012376579391412124], [0.5380878837333941, 0.21141064098174234, 0.14144037324212322, np.float64(0.3703198907514363), 1.5456250659930137], [0.3459591689840704, 0.6385804230981139, 1.7026429054443015, np.float64(-0.3222928792748565), 1.2620283818320337], [0.2796180286556356, 0.2794595354383969, 1.1162043416209115, np.float64(0.36577224625176097), 0.3445464048215645], [0.4083006099648495, 0.18056745029495191, 0.08768801746026589, np.float64(0.21076664846239274), 1.6243113582108568], [0.37744115692403224, 0.3035192193450436, 1.916932459697617, np.float64(0.2963424642545056), 1.5577998625247516], [0.2529412733776454, 0.2240367308446115, 1.3293598094721137, np.float64(0.6477186306682372), 0.5931734472070427], [0.1921897322956072, 0.4511447049753618, 0.32691044324662166, np.float64(-0.28049935564988376), 0.950257383851105], [0.5034849512181531, 0.3882567126350668, 1.2176379983119752, np.float64(0.8604013472607417), 0.9336767732120987], [0.10785131678733441, 0.23555933663973871, 1.5441898518498265, np.float64(1.1128460826811306), 0.8270085962827145], [0.3830587212862429, 0.49771856107555434, 0.37831050354988904, np.float64(-0.2717489494596904), 1.571592874474108], [0.46171878326383614, 0.46056076397270657, 0.2685390249994337, np.float64(-0.22010359530426848), 1.6546377880046177]]   
    print("NOW\n", list_inputs)

    list_inputs.append([0.2, 0.4, 1.3, -0.8, 1.4])
    list_inputs.append([0.2, 0.4, 2.6, -1.6, 2.8])

    OMIX1_list = []; OMIX2_list = []; OMIX3_list = []; FMIX1_list = []; FMIX2_list = []; FMIX3_list = []
    for i, test_set in enumerate(list_inputs):
        Z1, Z2, chi11, chi12, chi22 = test_set
        OMIX, FMIX, chi, chi_eta, x_prime, eta, Z_opt, Z_stoic = Q2DFGeneralizer(Z1, Z2, chi11, chi12, chi22, "toluene", "air", "n-heptane", i+1).get_optimal_mapping()

        Z1, Z2, chi11, chi12, chi22, chi, chi_eta, x_prime, eta, Z_opt, Z_stoic, OMIX1, OMIX2, OMIX3, FMIX1, FMIX2, FMIX3 = [str(float(i)) + "_WP" for i in [Z1, Z2, chi11, chi12, chi22, chi, chi_eta, x_prime, eta, Z_opt, Z_stoic, OMIX[0], OMIX[1], OMIX[2], FMIX[0], FMIX[1], FMIX[2]]]
        Z1_CFD_list.append(Z1); Z2_CFD_list.append(Z2); chi11_CFD_list.append(chi11); chi12_CFD_list.append(chi12); chi22_CFD_list.append(chi22)
        chi_xi_truth.append(chi); chi_eta_truth.append(chi_eta); x_prime_truth.append(x_prime); eta_opt_truth.append(eta); Z_opt_truth.append(Z_opt); Z_stoic_truth.append(Z_stoic)
        OMIX1_list.append(OMIX1); OMIX2_list.append(OMIX2); OMIX3_list.append(OMIX3); FMIX1_list.append(FMIX1); FMIX2_list.append(FMIX2); FMIX3_list.append(FMIX3)

    OMIX_truth = OMIX1_list + OMIX2_list + OMIX3_list
    FMIX_truth = FMIX1_list + FMIX2_list + FMIX3_list

    fortran_code = {"Z1_CFD_list":Z1_CFD_list,
                    "Z2_CFD_list":Z2_CFD_list,
                    "chi11_CFD_list":chi11_CFD_list,
                    "chi12_CFD_list":chi12_CFD_list,
                    "chi22_CFD_list":chi22_CFD_list,
                    "OMIX_truth":OMIX_truth,
                    "FMIX_truth":FMIX_truth,
                    "chi_xi_truth":chi_xi_truth,
                    "chi_eta_truth":chi_eta_truth,
                    "x_prime_truth":x_prime_truth,
                    "eta_opt_truth":eta_opt_truth,
                    "Z_opt_truth":Z_opt_truth,
                    "Z_stoic_truth":Z_stoic_truth}

    print("\nRESULTS BELOW:")
    strToPrint = ""
    for var in fortran_code.keys():
        if var in ["OMIX_truth", "FMIX_truth"]:
            strToPrint += "\t" + var + " = reshape([" + ", ".join(fortran_code[var]) + "], shape(" + var + "))" + "\n"
        else:
            strToPrint += "\t" + var + " = [" + ", ".join(fortran_code[var]) + "]" + "\n"
    
    strToPrint = strToPrint.split("\n")
    LENGTH = 80
    for i in strToPrint:
        chunks = textwrap.wrap(i, width=LENGTH, break_long_words=False)  # [i[j:j + LENGTH] for j in range(0, len(i), LENGTH)]
        print(" &\n\t\t".join(chunks))

def get_rand():
    import random
    
    Z1 = random.random()
    Z2 = random.random()
    Z3 = random.random()
    sumZ = Z1 + Z2 + Z3
    Z1 /= sumZ
    Z2 /= sumZ
    Z3 /= sumZ

    chi11 = 2 * random.random()
    chi22 = 2 * random.random()
    maxCrossChi = np.sqrt(chi11 * chi22)
    chi12 = -maxCrossChi + 2 * maxCrossChi * random.random()

    return Z1, Z2, chi11, chi12, chi22

def main():
    import time

    myQ2DFGen = Q2DFGeneralizer(0.1, 0.3, 3, -2.4, 4, "toluene", "air", "n-heptane")
    currentTime = time.time()

    NUM_TIMES = 40
    extrema_count = 0; extrema_count2 = 0
    for j in range(NUM_TIMES):
        Z1, Z2, chi11, chi12, chi22 = get_rand()

        extrema = myQ2DFGen.get_extrema_chi_ratio(Z1, Z2, chi11, chi12, chi22)
        extrema = [float(extremum) for extremum in extrema]
        extrema2 = myQ2DFGen.get_extrema_chi_ratio2(Z1, Z2, chi11, chi12, chi22)
        extrema2 = [float(extremum) for extremum in extrema2 if 0<=extremum and extremum<=1]
        print(extrema, extrema2)
        for extremum in extrema + extrema2:
            myQ2DFGen.assert_local_extremum(Z1, Z2, chi11, chi12, chi22, extremum)
            extrema_count += 1
            extrema_count2 += 1

    print("Avg num. extrema", extrema_count/NUM_TIMES, extrema_count2/NUM_TIMES)
    print("Avg. Time spent", (time.time() - currentTime) / NUM_TIMES)


def test():
    Z1_CFD = 0.334
    Z2_CFD = 0.129
    chi11_CFD = 3.743
    chi12_CFD = -1.348
    chi22_CFD = 2.1893

    myQ2DFGen = Q2DFGeneralizer(Z1_CFD, Z2_CFD, chi11_CFD, chi12_CFD, chi22_CFD, "toluene", "air", "n-heptane")
    chi, OMIX, FMIX, x_prime = myQ2DFGen.get_optimal_mapping()

    print("OMIX", OMIX)
    print("FMIX", FMIX)
    print("chi", chi)
    print("xprime", x_prime)

if __name__ == "__main__":
    unit_test_fortran()
