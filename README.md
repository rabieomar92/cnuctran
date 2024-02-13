# CNUCTRAN: A C++ nuclides depletion solver based on the probabilistic method.


CNUCTRAN is a C++ library created to simulate various nuclear transmutations such as decays, fissions, and neutron absorptions. Download CNUCTRAN user's manual here: https://github.com/rabieomar92/cnuctran/blob/main/cnuctran_manual.pdf

## Features
<div align="justify">

- Friendly to all physicists! In order to understand the method implemented in CNUCTRAN, you only need to know <a href="https://en.wikipedia.org/wiki/Poisson_distribution" target=_blank>Poisson distribution</a> and <a href="https://en.wikipedia.org/wiki/Matrix_multiplication" target=_blank>matrix multiplications</a>! Free from hard-to-understand algorithms and math jargons (of course, matrix multiplication is the minimum requirement).
- Capable of simulating complex transmutation chains.
- Based on the realistic transmutation processes.
- Helps nuclear physics students to understand transmutation processes through simulations. They can create, study and design any depletion chains and simulate the transmutations.

## How to Install
The current development of CNUCTRAN only provides a Windows executable. For other operating systems, you may have to compile it on your own. Unfortunately, you may have to replace the Microsoft Parallel Pattern Library (PPL) with any other compatible infrastructure, such as OpenMP. Note that CNUCTRAN requires PPL for its optimization.
The installation of CNUCTRAN involves the following simple steps:

**Step 1**: Visit https://github.com/rabieomar92/cnuctran/releases. Select the latest release and download the compressed file named cnuctran.zip.

**Step 2**: Extract the contents of <code>cnuctran.zip</code> into the same working directory. The content of <code>cnuctran.zip</code> is summarized below:

<small>
  
| File | Type | Purpose |
|----------|:----------|:----------|
cnuctran.exe | Executable File | The primary CNUCTRAN executable. |
<code>cmd.lnk</code> | Windows Link File | A shortcut file that opens a Command Prompt with the main CNUCTRAN directory set as the working directory. |
<code>mpfr.dll</code> | Dynamic Link Library (DLL) | A library file enabling multiple precision floating point arithmetic. |
<code>mpir.dll</code> | Dynamic Link Library (DLL) | A library file enabling multiple precision integer arithmetic. |
<code>chain_endfb71.xml</code> | Extensible Markup Language (XML) | The default ENDFBVII.1 data contains transmutation chain information such as daughter products, fission yields and half-lives. |
<code>input.xml</code> | Extensible Markup Language (XML) | A default input file for a test run, serving as the main model file requiring user modification. |
<code>CNUCTRAN Manual.pdf</code>	| PDF | CNUTRAN's user manual hanbook. |
<code>sample_inputs</code> | Directory | A directory containing four examples of test input XML files. |

</small>

**Step 3**: Double-click the <code>cmd</code> shortcut file. The command prompt console will pop up from the current working directory.

**Step 4**: Prepare the XML input file. You may rename the input XML file according to your needs. If you are a first-timer, you may skip this step because a predefined input file is already included in <code>cnuctran.zip</code>.

**Step 5**: If the input file name is <code>input.xml</code>, then type in <code>cnuctran</code> in the command prompt console, then press "enter". Or else, type in <code>cnuctran [anyname].xml</code>, then press enter.



## Program Flow
The program package consists of four (4) important files in order to run smoothly. They are the input file (<code>input.xml</code>) which consists of the problem definition set by the user; cnuctran.exe which is the executable of the compiled CNUCTRAN C++ codes; and finally, <code>mpfr.dll</code> and <code>mpir.dll</code> which are the dependencies of CNUCTRAN to enable high-precision arithmetics. Firstly, <code>cnuctran.exe</code> will read the simulation parameters defined in <code>input.xml</code>, including the arithmetic precision level. Here, the precision is defined in terms of the number of significant digits, n_sig, which can be set via the ```<precision_digits>n_sig</precision_digits>``` tag, inside the  ```<simulation_params></simulation_params>``` node. 

In the second step, the program will iterate over all zones defined by the user via the ```<zone></zone>``` tag. For each zone, the program reads the nuclide names involved in the calculation as well as the reaction rates and the initial nuclide concentrations. Once all parameters have been acquired, the transfer matrix is constructed for the zone. Finally, the program will run the calculation, and these processes repeat for all zones. The calculated species concentrations will be printed in the output file specified by the user via the ```<output>output_location</output>``` tag inside the  ```<simulation_params></simulation_params>``` of the input file. 

