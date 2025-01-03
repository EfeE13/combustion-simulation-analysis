# Functionality

Various things for NGA and offline PDRs are implemented here, including:
- Running and analyzing PDRs with a Cartesian product of different inputs
    - Explicit implementation for wide range of 2D cartesian product color maps
    - Has opportunity to easily extend to 1D graphing by creating a child class (PDRsRuns1D) of PDRsRunsND (n-dimensional Cartesian product) module's class in the pdrs_analysis package
- Plotting the relation between any three NGA output quantities within a time step via colored scatterplot
    - Functionality for automatically running the three quasi-dimensional flamelet (Q2DF) models for a chosen subset of the cells in a time step
        - Functionality for determining which cells yield stoichiometric mixture fraction in [0, 1] for a chosen Q2DF model

To run a code, one can go to the folder `CombustionModelAnalysis` and run something like `python -m src.pdrs_analysis.MainExample`

# Requirements

In addition to the source code, one should have:
- mechanisms directory for provided PDRs mechanisms
- inputs_pdrs directory to store generated PDRs inputs
- outputs_pdrs directory to store generated PDRs outputs
- graphs directory to store generated graphs
- some well-known Python libraries (e.g., `numpy`, `matplotlib`, `scipy`, `math`, `copy`, `struct`, `random`, `subprocess`, `os`, `itertools`, `time`)

**Note:** `main` file examples will be included.