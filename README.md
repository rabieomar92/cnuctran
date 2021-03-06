# CNUCTRAN: A C++ nuclides depletion solver based on the probabilistic method.

<div align="justify">
CNUCTRAN is a C++ library created to simulate various nuclear transmutations such as decays, fissions, and neutron absorptions. The code helps physicists avoid cumbersome numerical issues of solving the nuclide depletion equations (also known as the Bateman's equations). These issues include the stiffness of the Bateman's equations due to the complex decay chain problems. The Bateman's equation of a depletion problem is given as follows:
</div>

\
<img src="https://latex.codecogs.com/svg.image?\frac{dn_{i}}{dt}=&space;\underset{\mathrm{production}}{\underbrace{\sum_{j=1}^{m}b_{ij}\lambda_{j}n_{j}&plus;\sum_{k=1}^{m}y_{ik}\Lambda_{ik}n_{k}}}&space;&space;-&space;&space;&space;\underset{\mathrm{removal}}{\underbrace{\left&space;(&space;\lambda_{i}&plus;\Lambda_{i}&space;\right)n_{i}}&space;" title="\frac{dn_{i}}{dt}= \underset{\mathrm{production}}{\underbrace{\sum_{j=1}^{m}b_{ij}\lambda_{j}n_{j}+\sum_{k=1}^{m}y_{ik}\Lambda_{ik}n_{k}}} - \underset{\mathrm{removal}}{\underbrace{\left ( \lambda_{i}+\Lambda_{i} \right)n_{i}} " />

<img src="https://latex.codecogs.com/svg.image?\inline&space;n_i" title="\inline n_i" /> is the atom density of isotope-*i*;\
<img src="https://latex.codecogs.com/svg.image?\inline&space;\lambda_i" title="\inline \lambda_i" /> is the radioactive decay constant of isotope-*i* causing its removal from the system (the decay rate per second);\
<img src="https://latex.codecogs.com/svg.image?\inline&space;\lambda_{ij}" title="\inline \lambda_{ij}" /> is the decay constant of isotope-*j* causing the production of isotope-*i* in the system; \
<img src="https://latex.codecogs.com/svg.image?\inline&space;b_{ij}" title="\inline b_{ij}" /> is the decay branching ratio from isotope-*j* into isotope-*i*;\
<img src="https://latex.codecogs.com/svg.image?\inline&space;y_{ij}" title="\inline y_{ij}" /> is the production yield from isotope-*k* into isotope-*i*; and,\
<img src="https://latex.codecogs.com/svg.image?\inline&space;\Lambda_{ij}" title="\inline \Lambda_{ij}" /> is the production rate of isotope-*i* due to the removal of isotope-*k* from the system (per second).

<div align="justify">
It is also possible to simulate the actual transmutation processes using Monte Carlo method via iterations over a massive amount of nuclides. Alas, the simulation speed increases with increasing accuracy, which makes it not practical, and the simulation of billions of nuclides can take forever!  Thanks to the variance reduction technique in Monte Carlo method, the concept of isotope weight can be adopted. Here, the isotope weights can be adjusted and scaled using event probabilities, ??, which will be described later. The method is easy to understand, and it can be used by physicists from various mathematical background.
</div>

## Features

- Friendly to all physicists! In order to understand the method implemented in CNUCTRAN, you only need to know <a href="https://en.wikipedia.org/wiki/Poisson_distribution" target=_blank>Poisson distribution</a> and <a href="https://en.wikipedia.org/wiki/Matrix_multiplication" target=_blank>matrix multiplications</a>! Free from math jargons, hard-to-understand algorithms and approximations.
- Capable of simulating complex transmutation chains.
- Helps nuclear physics students to understands transmutation processes through simulations. They can create, study and design any depletion chains and simulate the transmutations.

If you don't prefer dealing with complicated mathematical methods to solve Bateman's equations,  give CNUCTRAN a try!

## Program Flow
The program package consists of four (4) important files in order to run smoothly. They are the input file (input.XML) which consists of the problem definition set by the user; cnuctran.exe which is the executable of the compiled CNUCTRAN C++ codes; and finally, mpfr.dll and mpir.dll which are the dependencies of CNUCTRAN to enable high-precision arithmetics. Firstly, cnuctran.exe will read the simulation parameters defined in input.XML, including the arithmetic precision level. Here, the precision is defined in terms of the number of significant digits, n_sig, which can be set via the ```<precision_digits>n_sig</precision_digits>``` tag, inside the  ```<simulation_params></simulation_params>``` node. 

In the second step, the program will iterate over all zones defined by the user via the ```<zone></zone>``` tag. For each zone, the program reads the species names involved in the calculation as well as the reaction rates and the initial species concentrations. Once all parameters have been acquired, the transfer matrix is constructed for the zone. Finally, the program will run the calculation, and these processes repeat for all zones. The calculated species concentrations will be printed in the output file specified by the user via the ```<output>output_location</output>``` tag inside the  ```<simulation_params></simulation_params>``` of the input file. 

## Summary of the Method

<div align="justify">
A transmutation process involves the removal of a nuclide from a system. Then it leads to the creation of another daughter nuclide. For instance, the decay of U-238 into Th-234 involves removing U-238 from the system via alpha decay, which in fact mutating the U-238 nucleus into Th-234. In reality, such a transmutation process occur at a certain rate, ??. CNUCTRAN works by first accumulating the removal parameters from the user. The removal parameters include the rate, parent isotope and the daughter isotope(s). Practically, any number of removals can be added to CNUCTRAN. Next, CNUCTRAN prepares the transfer matrix, <b>A</b>, a scaling factor that scales the initial isotope concentrations into the final concentrations. The construction of <b>A</b> requires the probabilities of various transmutations from one species to another.
</br></br>
</div>

<div align="justify">
  The simulation requires the division of timestep, <i>t</i> into <i>N</i> regular substeps, and the substep interval is given by ??<i>t</i>. Consider an isotope-<i>i</i> which is expecting to experience <i>J<sub>i</sub></i> removal events. For example, U-235 may experience absorption, fission and decay events, so, there are <i>J<sub>i</sub></i>=3 removal events. Let ??<sub><i>ij</i></sub> be the rate (per secs.) of removal event-<i>j</i> experienced by isotope <i>i</i>. For decay removal events, ??<sub><i>ij</i></sub> is the decay constant multiplied with the branching ratio of the decay branch, i.e. ??<sub><i>ij</i></sub> = b<sub><i>ij</i></sub></sub>??<sub><i>ij</i></sub>.  The removal probability (or probability of a removal event to occur) of isotope-<i>i</i> from a system due to <i>j</i>-th removal process can be derived from Poisson statistics, leading to an un-normalized joint Poisson distribution (later known as ??-distribution),
</div>

\
<img src="https://latex.codecogs.com/svg.image?\widetilde{\pi}_{il}&space;=&space;\prod_{j=1}^{J_i}\left\{&space;\delta_{lj}&plus;(-1)^{\delta_{lj}}&space;e^{-\Lambda_{j}\Delta&space;t}\right\}" title="\pi_{ij} = \prod_{l=1}^{J_i}\left\{ \delta_{lj}+(-1)^{\delta_{lj}} e^{-\Lambda_{l}\Delta t}\right\}" />

where ??<sub><i>il</i></sub> is the <a href="https://en.wikipedia.org/wiki/Kronecker_delta">Kronecker delta</a>. The probability of isotope-*i* for not being removed from the system is given by *l*=0

<img src="https://latex.codecogs.com/svg.image?\widetilde{\pi}_{i0}&space;=&space;\prod_{j=1}^{J_i}&space;e^{-\Lambda_{j}\Delta&space;t}" title="\pi_{i0} = \prod_{j=1}^{J_i} e^{-\lambda_{j}\Delta t}" />

The normalized probability of removal-*l* to occur is given by (*l*=0 is for no-removal, the case when no removal events happen):

<img src="https://latex.codecogs.com/svg.latex?\pi_{il}&space;=&space;\frac{\widetilde{\pi}_{il}&space;}{\sum_{j=0}^{J_{i}}\widetilde{\pi}_{ij}}" title="P_{il} = \frac{f_{il} }{\sum_{j=0}^{J_{i}}f_{ij}}" />

If event-<i>l</i> of species <i>i</i> is a fission reaction, then, the product's transfer probability must be scaled to the fission yield, <i>y<sub>ij</sub></i>,

<img src="https://latex.codecogs.com/svg.latex?\pi_{il}&space;=&space;\frac{\widetilde{\pi}_{il}&space;}{\sum_{j=0}^{J_{i}}\widetilde{\pi}_{ij}}y_{il}" title="P_{il} = \frac{f_{il} }{\sum_{j=0}^{J_{i}}f_{ij}}" />

Conveniently, the derived joint Poisson distribution is coined as the ??-distribution. At this point, we let <i>I</i> as the total number of species involved in the depletion problem, and we define the transfer matrix as 

<img src="https://latex.codecogs.com/svg.latex?\mathbf{A}&space;=&space;\begin{pmatrix}&space;\pi_{1\rightarrow&space;1}&space;&&space;\cdots&space;&&space;\pi_{I\rightarrow&space;1}&space;\\&space;\vdots&space;&&space;\ddots&space;&&space;\vdots&space;\\&space;\pi_{1\rightarrow&space;I}&space;&&space;\cdots&space;&&space;\pi_{I\rightarrow&space;I}&space;\end{pmatrix}" title="\mathbf{A} = \begin{pmatrix} \pi_{1\rightarrow 1} & \cdots & \pi_{I\rightarrow 1} \\ \vdots & \ddots & \vdots \\ \pi_{1\rightarrow I} & \cdots & \pi_{I\rightarrow I} \end{pmatrix}" />

where ??(k???i) is a scaling factor that scales the initial concentration of isotope-*i*. Also, ??(k???i) is defined as

<img src="https://latex.codecogs.com/svg.latex?\pi_{k&space;\rightarrow&space;i}&space;=&space;\sum_{j\in&space;R}^{}&space;\pi_{kj};\;&space;\;&space;\;\pi_{i&space;\rightarrow&space;i}&space;=&space;\pi_{i0}" title="\pi_{k \rightarrow i} = \sum_{j\in R}^{} \pi_{kj}" />

<div align="justify">
Here, R is a set of transmutation events that mutate species k into species i. Note that matrix A is a square matrix (<i>IxI</i>) with its columns as the parent species and the its rows as the daughter species. Now, let <b>w</b>(t) and <b>w</b>(0) be the column matrices representing the final and initial concentration of all species involved, respectively. Then, the final concentrations can be easily evaluated via,
</div>

\
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{w}(t)=\mathbf{A}^{t/\Delta&space;t}\mathbf{w}(0)&space;\quad&space;\textrm{where}&space;\quad&space;t/\Delta&space;t&space;\in&space;\mathbb{N}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{w}(t)=\mathbf{A}^{t/\Delta&space;t}\mathbf{w}(0)&space;\quad&space;\textrm{where}&space;\quad&space;t/\Delta&space;t&space;\in&space;\mathbb{N}" title="\mathbf{w}(t)=\mathbf{A}^{t/\Delta t}\mathbf{w}(0) \quad \textrm{where} \quad t/\Delta t \in \mathbb{N}" /></a>

It is important to remark that t/??t forms a very large natural number. So, the matrix exponentiation in the above equation can be evaluated at an incredible speed using the binary decomposition method, see https://cp-algorithms.com/algebra/binary-exp.html. Also, the transfer matrix A is a sparse matrix and most of its elements are zero. Alas, the calculation is very stiff because the transfer matrix A is constructed by using a wide range of nuclide half-lives. Thus, the computation requires high precision arithmetic. A special sparse matrix data structure is programmed in ```smatrix``` class for CNUCTRAN's specific use where the high-precision matrix elements are represented using the C++'s third party library, ```MPFR``` (see [MPFR](https://www.mpfr.org/)). The matrix power is also handled using ```smatrix::binpow(...)``` method.

## Derivation of ??-distribution

<div align="justify">
The Poisson probability distribution function appeared in many applications that model the number of events occurring within a given time interval, ??t. The abscissa of the distribution function is the number of chance(s) for the event to occur, <i>k</i>???{0,1,???,???}, occurring within ??t. The distribution is characterized by the parameter ??, which is the average number of events occurring within (t, t+??t). The formula for Poisson distribution is given by
 </div>

\
<img src="https://latex.codecogs.com/svg.latex?p_i\left&space;(&space;k;&space;\mu&space;\right&space;)&space;=&space;\frac{e^{-\mu&space;}&space;\mu^k}{k!}" title="p_i\left ( k; \Lambda \right ) = \frac{e^{-\Lambda } \Lambda^k}{k!}" />

<div align="justify">
Transmutations follow the Poisson process, and they consist of a series of nuclear events. The average time between events is known, but the exact timing of events is random. Most importantly, the events are independent. Here, one event does not affect the arrival of the other events. For instance, the decay event of a nuclide does not affect the decay event of other nuclides, since the decay rate is the intrinsic property of the species. 

 \
Suppose we select a nuclide of species i from a system during a particular time t. Depending on the environment, the nuclide can be associated to various possible transmutation events. It is convenient to index these possible events as j???{1,2,3,???,J_i}. Conveniently, one may assign event-1 as neutron absorption, event-2 as fission, event-3 as radioactive decay, and so on. Also, J_i is the total number of events defined for species i. For simplicity, j=0 is designated as the case when the nuclide does not undergo transmutation within the interval (t,t+??t).
 </div>
 
Now, consider that the rate, ??<sub><i>ij</i></sub>, for each possible event is known. Intuitively,

<img src="https://latex.codecogs.com/svg.latex?\mu_j&space;=&space;\Lambda_j&space;\Delta&space;t" title="\lambda_j = \Lambda_j \Delta t" />


<div align="justify">
The probability for a nuclide of species <i>i</i> to undergo reaction event-<i>j</i> depends on the number of chances for the reaction to occur, <i>k</i>. Of course, a single nuclide will only need at least one chance to occur, i.e., <i>k</i>???1. The probability of one or more chances for event-<i>j</i> to occur within the interval (t, t+??t) is given by the following summation,
</div>

\
<img src="https://latex.codecogs.com/svg.latex?p_{ij}&space;\equiv&space;\sum_{k=1}^{\infty&space;}&space;p_i(k;\Lambda_{ij})&space;=&space;\frac{\left&space;(&space;\Lambda_{ij}&space;\Delta&space;t&space;\right&space;)&space;^{1}}{1!}e^{-\Lambda_{ij}&space;\Delta&space;t}&space;&plus;&space;\frac{\left&space;(&space;\Lambda_{ij}&space;\Delta&space;t&space;\right&space;)&space;^{2}}{2!}&space;e^{-\Lambda_{ij}&space;\Delta&space;t}&space;&plus;&space;..." title="p_{ij} \equiv \sum_{k=1}^{\infty } p_i(k;\Lambda_{ij}) = \frac{\left ( \Lambda_{ij} \Delta t \right ) ^{1}}{1!}e^{-\Lambda_{ij} \Delta t} + \frac{\left ( \Lambda_{ij} \Delta t \right ) ^{2}}{2!} e^{-\Lambda_{ij} \Delta t} + ..." />

Since Poisson distribution normalizes to unity, then the simple representation of the above equation can be written as,

<img src="https://latex.codecogs.com/svg.latex?p_{ij}&space;=&space;1&space;-&space;\frac{\left&space;(&space;\Lambda_{ij}&space;\Delta&space;t&space;\right&space;)&space;^{0}}{0!}&space;e^{-\Lambda_{ij&space;\Delta&space;t}}&space;=&space;1&space;-&space;e^{-\Lambda_{ij}&space;\Delta&space;t}" title="p_l = 1 - \frac{\left ( \Lambda_l \Delta t \right ) ^{0}}{0!} e^{-\Lambda_{il \Delta t}} = 1 - e^{-\Lambda_{il} \Delta t}" />

Also, the probability of zero-chance of event-<i>j</i> to occur is given by,

<img src="https://latex.codecogs.com/svg.latex?q_{ij}&space;\equiv&space;p(k=0;\Lambda_{ij})&space;=&space;\frac{\left&space;(&space;\Lambda_{ij}&space;\Delta&space;t&space;\right&space;)&space;^{0}}{0!}&space;e^{-\Lambda_{ij&space;\Delta&space;t}}&space;=&space;e^{-\Lambda_{ij}&space;\Delta&space;t}" title="q_l \equiv p(k=0;\Lambda_l) = \frac{\left ( \Lambda_l \Delta t \right ) ^{0}}{0!} e^{-\Lambda_{il \Delta t}} = e^{-\Lambda_{il} \Delta t}" />


<div align="justify">
So far, we have established the event probability of a single event. Recall that an isotope species can be associated with more than one transmutation events, <i>j</i>???{0,1,2,???,<i>J<sub>i</sub></i>}. For simplicity, let E<sub>1</sub> be the condition where event-<i>l</i>, <i>l</i>???{0,1,2,???,<i>J<sub>i</sub></i>}, is occurring, and E<sub>2</sub> be a condition where all events except removal event-<i>l</i> are not occurring. The joint probability of both E<sub>1</sub> and E<sub>2</sub> are true is given by,
</div>

\
<img src="https://latex.codecogs.com/svg.latex?P(E_1&space;\cap&space;E_2)&space;=&space;p_{il}&space;\times&space;\prod_{\forall&space;j&space;\ne&space;l}^{}&space;q_{ij}" title="P(E_1 \cap E_2) = p_l \times \prod_{\forall j \ne l}^{} q_j" />

By using the given definition of <i>p<sub>ij</sub></i> and <i>q<sub>ij</sub></i>, the un-normalized probability distribution can be expressed as,

<img src="https://latex.codecogs.com/svg.latex?\widetilde{\pi}_{il}&space;=&space;\left&space;(&space;1&space;-&space;e^{-\Lambda_{il}&space;\Delta&space;t}&space;\right&space;)\times&space;\prod_{\forall&space;l&space;\ne&space;j}^{}&space;e^{-\Lambda_{il}&space;\Delta&space;t}" title="\widetilde{\pi}_{il} = \left ( 1 - e^{-\Lambda_{il} \Delta t} \right )\times \prod_{\forall l \ne j}^{} e^{-\Lambda_{il} \Delta t}" />

The probability of no-removal, <i>l</i>=0 is given by the product of zero-chance probability of all possible events,

<img src="https://latex.codecogs.com/svg.latex?\widetilde{\pi}_{i0}&space;=&space;\prod_{\forall&space;l}^{}&space;e^{-\Lambda_{il}&space;\Delta&space;t}" title="\widetilde{\pi}_{il} = \prod_{\forall l}^{} e^{-\Lambda_l \Delta t}" />

Therefore, the normalized joint Poisson distribution can be formally written as,


<img src="https://latex.codecogs.com/svg.latex?\pi_{il}&space;\equiv&space;\pi_i(l;&space;\Lambda_1,&space;...,&space;\Lambda_{J_i})&space;=&space;c_i&space;\prod_{j=1}^{J_i}&space;\left&space;(&space;\delta_{jl}&space;&plus;&space;\left&space;(&space;-1&space;\right&space;)^{\delta_{jl}}&space;e^{-\Lambda_j&space;\Delta&space;t}&space;\right&space;)" title="\pi_{il} \equiv \pi_i(l; \Lambda_1, ..., \Lambda_{J_i}) = c_i \prod_{j=1}^{J_i} \left ( \delta_{jl} + \left ( -1 \right )^{\delta_{jl}} e^{-\Lambda_j \Delta t} \right )" />

where c<sub><i>i</i></sub> is the normalization constant and ??<sub><i>jl</i></sub> is the Kronecker delta. The value of c<sub><i>i</i></sub> is given by,


<img src="https://latex.codecogs.com/svg.latex?c_i&space;=&space;\left&space;(&space;\sum_{l=0}^{J_i}&space;\widetilde{\pi}_{il}&space;\right&space;)^{-1}" title="c_i = \left ( \sum_{l=0}^{J_i} \widetilde{\pi}_l \right )^{-1}" />

## Input Example.

```xml
<?xml version="1.0" encoding="utf-8"?>

<problem>

  <!-- Calculation zones. -->
	<zone name="myzone">
		<!--
			 ZONE DEFINITIONS.
			 
			 Each zone has its own problem definition. For instance, a reactor fuel pin can be considered as
			 a zone. Users are allowed to define as many zones as they wish.
	
		-->
		
		<species>
			<!--
				SPECIES DEFINITIONS.
			
				Each zone must have a list of species names to enable the calculation.

				DATA SOURCE.
				The user can provide the nuclides data via the following attribute:
				<species source=".\nuclides\data\source\location"></species>

				SPECIES NAME SPECIFICATION.
				The user may specify the names manually, OR, specify the range of species names:
				
				For example, if the user wishes to include all nuclides with atomic number within
				[100,300], then add the following attribute,
				
    				<species source="nuclides\data\source\" amin="100" amax="300">...</species>

				However, the nuclides data source file must be specified. Specifying amin and amax 
				attributes will ignore the species names listed in the species node. 
			-->
				
			Np237 Pa233 U233 Th229 Ra225 Ac225 Fr221 At217 Bi213 Po213 Tl209 Pb209 Bi209
		</species>
		
		<initial_concentrations>
			<!--
				INITIAL CONCENTRATIONS SPECIFICATION.
				
				The user may define the initial concentration of each individual nuclide manually via
				the <n0 species="species_name" value="#####" /> tag.
 				The user can also import the previously calculated concentrations stored in an XML file.
				Note that the output XML of the previous calculations can serve as the initial concentration
				source file. To import the existing concentration file, the following attribute must exist,

				<initial_concentrations source=".\path\to\initial\concentration\file.xml"></initial_concentrations>
			-->
			
			<concentration species="Np237" value="1.0" />
		</initial_concentrations>
		
		<reaction_rates>
			<reaction species="U238" type="(n,gamma)" rate="1e-4" />
			<reaction species="U238" type="fission" rate="1e-5" />
		</reaction_rates>
		
		<removals>
			<removal rate="1.027082E-14" parent="Np237" daughters="Pa233" yields="1.0" />
			<removal rate="2.971055E-07" parent="Pa233" daughters="U233" yields="1.0" />
			<removal rate="1.380625E-13" parent="U233" daughters="Th229" yields="1.0" />
			<removal rate="2.994544E-12" parent="Th229" daughters="Ra225" yields="1.0" />
			<removal rate="5.384253E-07" parent="Ra225" daughters="Ac225" yields="1.0" />
			<removal rate="8.022537E-07" parent="Ac225" daughters="Fr221" yields="1.0" />
			<removal rate="2.406761E-03" parent="Fr221" daughters="At217" yields="1.0" />
			<removal rate="2.166085E+02" parent="At217" daughters="Bi213" yields="1.0" />
			<removal rate="0.00000519239" parent="Bi213" daughters="Tl209" yields="1.0" />
			<removal rate="0.00024324741" parent="Bi213" daughters="Po213" yields="1.0" />
			<removal rate="1.863299E+05" parent="Po213" daughters="Pb209" yields="1.0" />
			<removal rate="5.251115E-03" parent="Tl209" daughters="Pb209" yields="1.0" />
			<removal rate="5.924335E-05" parent="Pb209" daughters="Bi209" yields="1.0" />
		</removals>
	</zone>

  <!-- Simulation parameters. -->
	<simulation_params>
		<!-- 
		
			 SIMULATION PARAMETERS DESCRIPTION.
			 
			 n is the order of the calculation. This simply means that the substep size is below 10^-n sec.
			 time_step is the time step of the calculation.
			 precision_digits is the minimum number of accurate digits maintained in the arithmetics.
			 output_digits is the number of decimal points to be displayed in the output.
			 verbosity is the verbosity level. A higher verbosity means more console output messages will be displayed.
			 output is the location of the output file on the disk.
			 max_rate is the maximum removal rate to be considered in the calculation.*
			 min_rate is the minimum removal rate to be considered in the calculation.**
			 
			 *Any removal rates that fall above max_rate will be ignored by the program.
			 **Any removal rates that fall below min_rate will be ignored by the program.
			 
		-->
		<n>12</n>
		<time_step>31556952000000</time_step>
		<precision_digits>40</precision_digits>
		<output_digits>15</output_digits>
		<verbosity>2</verbosity>
		<output>.\output.xml</output>
		<max_rate>1000000000</max_rate>
		<min_rate>1e-200</min_rate>
	</simulation_params>
</problem>
```

