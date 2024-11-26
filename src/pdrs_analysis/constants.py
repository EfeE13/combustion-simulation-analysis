VARY_SDR_FMX = '''
Manifold type : Zeta FBC BD
Mechanism file : /home/efeeroz/Documents/CombustionModelAnalysis/mechanisms/nc7h16_stanford_reduced/RedMechU.chmech
Thermodata file : /home/efeeroz/Documents/CombustionModelAnalysis/mechanisms/nc7h16_stanford_reduced/RedMechU.chthermo
Transdata file : /home/efeeroz/Documents/CombustionModelAnalysis/mechanisms/nc7h16_stanford_reduced/RedMechU.chtrans
Pressure : 1.0133e5
Oxidizer temperature : 300.0
Oxidizer composition : O2 0.23292
    N2 0.76708
Fuel temperature : 300.0
Fuel composition : NXC7H16 1.0
Fuel 2 temperature : 300.0
Fuel 2 composition : A1CH3 1.0
Fuel species : NXC7H16 A1CH3
Initial guess type : Equilibrium
Relative tolerance : 1.0e-6 !1.0e-8
Absolute tolerance : 1.0e-12
Verbose mode : .true.
Solver : Banded
Use ajac : .false.
PDRs output type : File
Use steady : .false.
Grid type : stretched
Beta Z : 4
Number of grid points : 31
Write to screen : .true.
Write mass fractions : .true.
Write flamelet file : .true.
Output directory : ~/Documents/CombustionModelAnalysis/outputs_pdrs/pdrs_cartesian_product/24-10-15
Output interval : 100
Output phi : 1
Return type : .Y
Scalar dissipation rate {SDR}: 1.0|2.0|3.0|4.0|5.0|10.0|20.0|30.0|40.0|50.0
FMIX {FMX}: 0.0|0.05|0.1|0.2|0.3|0.4|0.5|0.6|0.7|0.8|0.9|0.95|1.0
Mechanism format : CK
'''