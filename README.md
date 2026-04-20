
# MixPI: Software for atomistic PIMD simulations with mixed time slicing

<img src="images/toc.png" alt="Logo" width="300" />


MixPI is an atomistic molecular dynamics software package which conducts path-integral molecular dynamics and ring polymer molecular dynamics simulations using the mixed-time slicing approach. Mixed-time slicing PIMD/RPMD allows for the quantization of specific atoms in the simulation to different degrees; 
this approach can significantly reduce the computational cost of PIMD/RPMD simulations, especially for atomistic simulation with large system sizes. MixPI is written in Fortran and interfaces with the molecular dynamics software
CP2K for force and energy evaluations. See the citations below for more information on the mixed-time slicing PIMD method.

# Installation
Detailed installation instructions are provided in the PDF MixPI Manual in the repository. MixPI requires the installation of CP2K, which can be found at https://github.com/cp2k/cp2k/blob/master/INSTALL.md. 
The main repository contains the most recent developer version of MixPI that may contain additional features or parameters that are not outlined in the MixPI manual. MixPI v2.1 (the most recent version published on the GitHub page) is the most recent version 
of MixPI to have undergone limited compile and runtime checks and is up to date with the PDF Manual in the repository. MixPI is a code under active development and may contain bugs or errors. 


# Citation
An introduction to MixPI and some example simulations have been published at 
Johnson, Britta A., Siyu Bu, Christopher J. Mundy, and Nandini Ananth. "MixPI: Mixed-Time Slicing Path Integral Software for Quantized Molecular Dynamics Simulations." 
arXiv preprint arXiv:2411.11988 (2024).

The MixPI software is citable at https://doi.org/10.5281/zenodo.11130634. 

# Acknowledgements
The authors are grateful to Greg Schenter for useful discussions. 
B.A.J. and C.J.M acknowledge support by the DOE Office of Science, Office of Basic Energy Sciences, Division of Chemical Sciences, Geosciences, and Biosciences, Condensed Phase and Interfacial Molecular Science program, FWP 16249. PNNL is operated by Battelle Memorial Institute for the United States Department of Energy under DOE Contract No. DE-AC05-76RL1830. Research reported in this publication was also partially supported by the National Institute of General Medical Sciences of the National Institutes of Health under award number 1R35GM158400-01 to N.A., S. B. thanks the Cornell Arts and Sciences New Frontiers Grant for seed funding. In addition, N.A. and B.A.J. acknowledge prior support from the U.S. Department of Energy, Office of Basic Energy Sciences, Division of Chemical Sciences, Geosciences and Biosciences under Award DE-FG02-12ER16362 (Nanoporous Materials Genome: Methods and Software to Optimize Gas Storage, Separations, and Catalysis).  This research used resources of the National Energy Research Scientific Computing Center (NERSC), a Department of Energy User Facility using NERSC award BES-ERCAP0036077.
