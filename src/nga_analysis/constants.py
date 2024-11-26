# Below is the PDRs input file
# Zeta, added end time, removed use ajac

CLIPPING_VALUE_COLOR_MAP = 10e10
LIST_NONFUELS = ["O2", "N2"]
# If PDRs freezes on a filename, you can add to this list to skip it:
BAD_INPUT_FILES_RUNNING_LIST = ["TolAirHepZetaQ2DF3ZMIX0~08188FMIX0~24048"] #"TolAirHepManiQ2DFZMIX0~37356FMIX0~14037", "TolAirHepManiQ2DFZMIX0~7122FMIX0~23566.multicomponent", "TolAirHepZetaQ2DF3ZMIX0~08188FMIX0~24048"]

ZETA_INPUT_ZETA_FOLDER = '''
! ZMIX Value = ZMIX_INPUT
! ZSTAR Value = ZSTAR_INPUT
! FMIX Value = FMIX_INPUT
! CHI Value = CHI_INPUT
! C_ZZST Value = C_ZZST_INPUT
! C_ZST Value = C_ZST_INPUT
! BestQ2DFModel (possibly neglecting some models) = BestQ2DFModel_INPUT

Manifold type : Zeta
Mechanism file : /home/efeeroz/Documents/CombustionModelAnalysis/mechanisms/nc7h16_stanford_reduced/RedMechU.chmech
Thermodata file : /home/efeeroz/Documents/CombustionModelAnalysis/mechanisms/nc7h16_stanford_reduced/RedMechU.chthermo
Transdata file : /home/efeeroz/Documents/CombustionModelAnalysis/mechanisms/nc7h16_stanford_reduced/RedMechU.chtrans
Pressure : 1.0133e5
Oxidizer temperature : OxidTemp_INPUT
Oxidizer composition : OxidBC_INPUT
Fuel temperature : FuelTemp_INPUT
Fuel composition : FuelBC_INPUT
Fuel species : FuelSpecies_INPUT
Initial guess type : Equilibrium
End time : 1.0e+03
Relative tolerance : 1.0e-6 !1.0e-8
Absolute tolerance : 1.0e-12
Solver : Banded
PDRs output type : File
Use steady : .false.
Grid type : stretched
Beta Z : 4
Number of grid points : 31
Write to screen : .true.
Write mass fractions : .true.
Write flamelet file : .true.
Verbose mode : .true.
Output directory : outputs_pdrs/24_11_03Zeta
Output phi : 1
Return type : .Y
Scalar dissipation rate : ScalarDissRate_INPUT
Mechanism format : CK'''.strip()

ZETA_INPUT_Q2DF_FOLDER = ZETA_INPUT_ZETA_FOLDER.replace("outputs_pdrs/24_11_03Zeta", "outputs_pdrs/24_11_15TESTQ2DF")

MANI_Q2DF_INPUT = '''
! BestQ2DFModel (possibly neglecting some models) = BestQ2DFModel_INPUT

Manifold type : Q2DF
Mechanism file : /home/efeeroz/Documents/CombustionModelAnalysis/mechanisms/nc7h16_stanford_reduced/RedMechU.chmech
Thermodata file : /home/efeeroz/Documents/CombustionModelAnalysis/mechanisms/nc7h16_stanford_reduced/RedMechU.chthermo
Transdata file : /home/efeeroz/Documents/CombustionModelAnalysis/mechanisms/nc7h16_stanford_reduced/RedMechU.chtrans
Pressure : 1.0133e5

ZMIX : ZMIX_INPUT
ZSTAR : ZSTAR_INPUT
Scalar dissipation rate ZMIX : CHI_INPUT
Cross scalar dissipation rate ZMIX ZSTAR : C_ZZST_INPUT
Scalar dissipation rate ZSTAR : C_ZST_INPUT
Z1 temperature : 300.0
Z1 composition : A1CH3 1.0
Z2 temperature : 300.0
Z2 composition : O2 0.23292
    N2 0.76708
Z3 temperature : 300.0
Z3 composition : NXC7H16 1.0
Fuel species : A1CH3 NXC7H16
Initial guess type : Equilibrium
End time : 1.0e+03
Relative tolerance : 1.0e-6 !1.0e-8
Absolute tolerance : 1.0e-12
Solver : Banded
PDRs output type : File
Use steady : .false.
Grid type : stretched
Beta Z : 4
Number of grid points : 31
Write to screen : .true.
Write mass fractions : .true.
Write flamelet file : .true.
Verbose mode : .true.
Output directory : outputs_pdrs/24_11_15TESTQ2DF
Output interval : 100
Output phi : 1
Return type : .Y
Mechanism format : CK
'''.strip()