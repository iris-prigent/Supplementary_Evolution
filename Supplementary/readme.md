This file describes the files and procedures used to produce results from the article "The moulding of intra-specific trait variation by selection under ecological inheritance" by Iris Prigent (iris.prigent@unil.ch) and Charles Mullon (charles.mullon@unil.ch) for the journal Evolution.



The folder "mathematica_files" contains the Mathematica notebooks and the data used to perform the mathematical analyses of the model investigating the co-evolution the attack rate with dispersal. 

Files description:

- The notebook "Supplementary_material_Evolution.nb" contains the code used to make Figure 2 and Figure S4.A-B of the manuscript. 
- The notebook "Computation_hessian.nb" contains the code used to generate the .txt files used to make Figure 2.B and S4.B.
- The .txt files contain the elements of the Hessian matrix used to make Figure 2.B and S4.B, computed using the notebook "Computation_hessian.nb". 



The folders "code_cpp_fixed_demography" and "code_cpp_explicit_demography" contain the simulation programs used in the study. "code_cpp_fixed_demography" contains the programs used to simulate a population with fixed demography, which we used to make Figure 3-4, Figure S4.C-D, Figure S5 and Figure S6 of the manuscript. "code_cpp_explicit_demography" contains the programs used to simulate a population with fluctuating demography, which we used to make Figure S-6 of the manuscript. All the files required for compiling and running the program are given in those folders. The program has to be compiled using -std=c++11 and was written on Windows architecture. "code_cpp_fixed_demography" and "code_cpp_explicit_demography" contain the following files (relevant files are annotated). 


Files description:

- "main.cpp" is the main file of the program. It contains the sets of parameters and calls the functions used. You can use this file to change the parameters to run the simulations.
- "ranbin.cpp" contains the functions we use to simulate various distributions using a random numbers generator.
- "functions.cpp" contains all the functions used to simulate and sample the population. The function "iteration" iterates a generation, following the procedure described in Appendix F4.
- "mt.h" is the MersenneTwister random numbers generator.
- "header.h" contains all the prototypes of functions and structures used in the program.


Output files description:
 
- "variable.txt" contains for each generation the mean and variance in attack rates, the mean and variance in dispersal, and the mean and variance in resource abundance in the entire population. It also contains the total size of the population when demography fluctuates.
- "sampling_z_generation_x.txt" contains the values of the attack rate of all individuals in the population at generation x
- "sampling_d_generation_x.txt" contains the values of the dispersal of all individuals in the population at generation x
- "sampling_z_generation_x.txt" contains the values of the resource abundance of all patches in the population at generation x
- "sampling_i_generation_x.txt" contains the identity of the patch each individuals is in, for the entire population at generation x. This output is only generated when population size flutuates, in order to keep track of the size of each patch.

The folder "simulated_results" contains the output files of the simulations used to make Figure 3-4, S4-S5-S6. For each figure, the parameters used to generate the output are specified in the legend. 


Running the simulations:

First, you must make sure that the folders where the output files will be written exist (you need to give a correct path at the placeholder 'PATH' in the file "function.cpp").
Then the program was compiled using the following command: g++ -std=c++11 -o progam_name main.cpp functions.cpp ranbin.cpp
The program can then be launch from the working directory using the command line : program_name

Note: this might be a different command line if you are using a Ubuntu or Mac architecture. 


