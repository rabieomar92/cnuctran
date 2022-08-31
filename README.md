# CNUCTRAN: A C++ nuclides depletion solver based on the probabilistic method.

<div align="justify">
CNUCTRAN is a C++ library created to simulate various nuclear transmutations such as decays, fissions, and neutron absorptions. Download CNUCTRAN user's manual here: https://github.com/rabieomar92/cnuctran/edit/main/cnuctran_manual.pdf

## Features

- Friendly to all physicists! In order to understand the method implemented in CNUCTRAN, you only need to know <a href="https://en.wikipedia.org/wiki/Poisson_distribution" target=_blank>Poisson distribution</a> and <a href="https://en.wikipedia.org/wiki/Matrix_multiplication" target=_blank>matrix multiplications</a>! Free from math jargons, hard-to-understand algorithms and approximations.
- Capable of simulating complex transmutation chains.
- Helps nuclear physics students to understand transmutation processes through simulations. They can create, study and design any depletion chains and simulate the transmutations.

## Program Flow
The program package consists of four (4) important files in order to run smoothly. They are the input file (input.XML) which consists of the problem definition set by the user; cnuctran.exe which is the executable of the compiled CNUCTRAN C++ codes; and finally, mpfr.dll and mpir.dll which are the dependencies of CNUCTRAN to enable high-precision arithmetics. Firstly, cnuctran.exe will read the simulation parameters defined in input.XML, including the arithmetic precision level. Here, the precision is defined in terms of the number of significant digits, n_sig, which can be set via the ```<precision_digits>n_sig</precision_digits>``` tag, inside the  ```<simulation_params></simulation_params>``` node. 

In the second step, the program will iterate over all zones defined by the user via the ```<zone></zone>``` tag. For each zone, the program reads the species names involved in the calculation as well as the reaction rates and the initial species concentrations. Once all parameters have been acquired, the transfer matrix is constructed for the zone. Finally, the program will run the calculation, and these processes repeat for all zones. The calculated species concentrations will be printed in the output file specified by the user via the ```<output>output_location</output>``` tag inside the  ```<simulation_params></simulation_params>``` of the input file. 

