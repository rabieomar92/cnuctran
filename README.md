# CNUCTRAN: A Probabilistic Solver for Nuclear Transmutation Problems.

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
It is also possible to simulate the actual transmutation processes using Monte Carlo method via iterations over a massive amount of nuclides. Alas, the simulation speed increases with increasing accuracy, which makes it not practical. Thanks to the variance reduction technique in Monte Carlo method, the concept of isotope weight is adopted. Here, the isotope weights can be adjusted and scaled using event probabilities, π, which will be described later. The method is easy to understand, and it can be used by physicists from various mathematical background.
</div>

## Features

- Friendly to all physicists! In order to understand the method implemented in CNUCTRAN, you only need to know <a href="https://en.wikipedia.org/wiki/Poisson_distribution" target=_blank>Poisson distribution</a> and <a href="https://en.wikipedia.org/wiki/Matrix_multiplication" target=_blank>matrix multiplications</a>! Free from math jargons, hard-to-understand algorithms and approximations.
- Capable of simulating complex transmutation chains.
- Helps nuclear physics students to understands transmutation processes through simulations. They can create, study and design any depletion chains and simulate the transmutations.

If you don't prefer dealing with complicated mathematical methods to solve Bateman's equations,  give CNUCTRAN a try!

## Summary of the Method

<div align="justify">
A transmutation process involves the removal of a nuclide from a system. Then it leads to the creation of another daughter nuclide. For instance, the decay of U-238 into Th-234 involves removing U-238 from the system via alpha decay, which in fact mutating the U-238 nucleus into Th-234. In reality, such a transmutation process occur at a certain rate, Λ. CNUCTRAN works by first accumulating the removal parameters from the user. The removal parameters include the rate, parent isotope and the daughter isotope(s). Practically, any number of removals can be added to CNUCTRAN. Next, CNUCTRAN prepares the transfer matrix, <b>A</b>, that stores the probabilities of various transmutations from one species to another.
</br></br>
</div>

<div align="justify">
  The simulation requires the division of timestep, <i>t</i> into <i>N</i> regular substeps, and the substep interval is Δ<i>t</i>. Consider an isotope-<i>i</i> which is expecting to experience <i>J<sub>i</sub></i> removal events. For example, U-235 may experience absorption, fission and decay events, so, there are <i>J<sub>i</sub></i>=3 removal events. Let Λ<sub><i>ij</i></sub> be the rate (per secs.) of removal event-<i>j</i> experienced by isotope <i>i</i>. For decay removal events, Λ<sub><i>ij</i></sub> is the decay constant multiplied with the branching ratio of the decay branch, i.e. Λ<sub><i>ij</i></sub> = b<sub><i>ij</i></sub></sub>λ<sub><i>ij</i></sub>.  The removal probability (or probability of a removal event to occur) of isotope-<i>i</i> from a system due to <i>j</i>-th removal process can be derived from Poisson statistics, leading to an un-normalized joint Poisson distribution (later known as π-distribution),
</div>

\
<img src="https://latex.codecogs.com/svg.image?\widetilde{\pi}_{il}&space;=&space;\prod_{j=1}^{J_i}\left\{&space;\delta_{lj}&plus;(-1)^{\delta_{lj}}&space;e^{-\Lambda_{j}\Delta&space;t}\right\}" title="\pi_{ij} = \prod_{l=1}^{J_i}\left\{ \delta_{lj}+(-1)^{\delta_{lj}} e^{-\Lambda_{l}\Delta t}\right\}" />

where δ<sub><i>il</i></sub> is the <a href="https://en.wikipedia.org/wiki/Kronecker_delta">Kronecker delta</a>. The probability of isotope-*i* for not being removed from the system is given by *l*=0

<img src="https://latex.codecogs.com/svg.image?\widetilde{\pi}_{i0}&space;=&space;\prod_{j=1}^{J_i}&space;e^{-\Lambda_{j}\Delta&space;t}" title="\pi_{i0} = \prod_{j=1}^{J_i} e^{-\lambda_{j}\Delta t}" />

The normalized probability of removal-*l* to occur is given by (*l*=0 is for no-removal, the case when no removal events happen):

<img src="https://latex.codecogs.com/svg.latex?\pi_{il}&space;=&space;\frac{\widetilde{\pi}_{il}&space;}{\sum_{j=0}^{J_{i}}\widetilde{\pi}_{ij}}" title="P_{il} = \frac{f_{il} }{\sum_{j=0}^{J_{i}}f_{ij}}" />

If event-<i>l</i> of species <i>i</i> is a fission reaction, then, the product's transfer probability must be scaled to the fission yield, <i>y<sub>ij</sub></i>,

<img src="https://latex.codecogs.com/svg.latex?\pi_{il}&space;=&space;\frac{\widetilde{\pi}_{il}&space;}{\sum_{j=0}^{J_{i}}\widetilde{\pi}_{ij}}y_{il}" title="P_{il} = \frac{f_{il} }{\sum_{j=0}^{J_{i}}f_{ij}}" />

Conveniently, the derived joint Poisson distribution is coined as the π-distribution. At this point, we let <i>I</i> as the total number of species involved in the depletion problem, and we define the transfer matrix as 

<img src="https://latex.codecogs.com/svg.latex?\mathbf{A}&space;=&space;\begin{pmatrix}&space;\pi_{1\rightarrow&space;1}&space;&&space;\cdots&space;&&space;\pi_{I\rightarrow&space;1}&space;\\&space;\vdots&space;&&space;\ddots&space;&&space;\vdots&space;\\&space;\pi_{1\rightarrow&space;I}&space;&&space;\cdots&space;&&space;\pi_{I\rightarrow&space;I}&space;\end{pmatrix}" title="\mathbf{A} = \begin{pmatrix} \pi_{1\rightarrow 1} & \cdots & \pi_{I\rightarrow 1} \\ \vdots & \ddots & \vdots \\ \pi_{1\rightarrow I} & \cdots & \pi_{I\rightarrow I} \end{pmatrix}" />

where π(k→i) is the transfer probability which is defined as

<img src="https://latex.codecogs.com/svg.latex?\pi_{k&space;\rightarrow&space;i}&space;=&space;\sum_{j\in&space;R}^{}&space;\pi_{kj};\;&space;\;&space;\;\pi_{i&space;\rightarrow&space;i}&space;=&space;\pi_{i0}" title="\pi_{k \rightarrow i} = \sum_{j\in R}^{} \pi_{kj}" />

<div align="justify">
Here, R is a set of transmutation events that mutate species k into species i. Note that matrix A is a square matrix (<i>IxI</i>) with its columns as the parent species and the its rows as the daughter species. Now, let <b>w</b>(t) and <b>w</b>(0) be the column matrices representing the final and initial concentration of all species involved, respectively. Then, the final concentrations can be easily evaluated via,
</div>

\
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{w}(t)=\mathbf{A}^{t/\Delta&space;t}\mathbf{w}(0)&space;\quad&space;\textrm{where}&space;\quad&space;t/\Delta&space;t&space;\in&space;\mathbb{N}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\mathbf{w}(t)=\mathbf{A}^{t/\Delta&space;t}\mathbf{w}(0)&space;\quad&space;\textrm{where}&space;\quad&space;t/\Delta&space;t&space;\in&space;\mathbb{N}" title="\mathbf{w}(t)=\mathbf{A}^{t/\Delta t}\mathbf{w}(0) \quad \textrm{where} \quad t/\Delta t \in \mathbb{N}" /></a>

It is important to remark that the matrix power in the above equation can be evaluated at an incredible speed using the binary decomposition method, see https://cp-algorithms.com/algebra/binary-exp.html. Also, the transfer matrix A is a sparse matrix and most of its elements are zero. The ```scipy.sparse``` class could be used to cater sparse matrix multiplication and power, however, the class lacks arithmetic precisions. Thus, a special sparse matrix data structure is programmed in ```pynuctran.sparse.smatrix``` class, for CNUCTRAN's specific use where the high-precision matrix elements are represented using the C++'s ```decimal``` class, see https://docs.python.org/3/library/decimal.html. The matrix power is also handled using ```pynuctran.sparse.smatrix.bpow(...)``` method.

## Derivation of π-distribution

<div align="justify">
The Poisson probability distribution function appeared in many applications that model the number of events occurring within a given time interval, Δt. The abscissa of the distribution function is the number of chance(s) for the event to occur, <i>k</i>∈{0,1,…,∞}, occurring within Δt. The distribution is characterized by the parameter μ, which is the average number of events occurring within (t, t+Δt). The formula for Poisson distribution is given by
 </div>

\
<img src="https://latex.codecogs.com/svg.latex?p_i\left&space;(&space;k;&space;\mu&space;\right&space;)&space;=&space;\frac{e^{-\mu&space;}&space;\mu^k}{k!}" title="p_i\left ( k; \Lambda \right ) = \frac{e^{-\Lambda } \Lambda^k}{k!}" />

<div align="justify">
Transmutations follow the Poisson process, and they consist of a series of nuclear events. The average time between events is known, but the exact timing of events is random. Most importantly, the events are independent. Here, one event does not affect the arrival of the other events. For instance, the decay event of a nuclide does not affect the decay event of other nuclides, since the decay rate is the intrinsic property of the species. 

 \
Suppose we select a nuclide of species i from a system during a particular time t. Depending on the environment, the nuclide can be associated to various possible transmutation events. It is convenient to index these possible events as j∈{1,2,3,…,J_i}. Conveniently, one may assign event-1 as neutron absorption, event-2 as fission, event-3 as radioactive decay, and so on. Also, J_i is the total number of events defined for species i. For simplicity, j=0 is designated as the case when the nuclide does not undergo transmutation within the interval (t,t+Δt).
 </div>
 
Now, consider that the rate, Λ<sub><i>ij</i></sub>, for each possible event is known. Intuitively,

<img src="https://latex.codecogs.com/svg.latex?\mu_j&space;=&space;\Lambda_j&space;\Delta&space;t" title="\lambda_j = \Lambda_j \Delta t" />


<div align="justify">
The probability for a nuclide of species <i>i</i> to undergo reaction event-<i>j</i> depends on the number of chances for the reaction to occur, <i>k</i>. Of course, a single nuclide will only need at least one chance to occur, i.e., <i>k</i>≥1. The probability of one or more chances for event-<i>j</i> to occur within the interval (t, t+Δt) is given by the following summation,
</div>

\
<img src="https://latex.codecogs.com/svg.latex?p_{ij}&space;\equiv&space;\sum_{k=1}^{\infty&space;}&space;p_i(k;\Lambda_{ij})&space;=&space;\frac{\left&space;(&space;\Lambda_{ij}&space;\Delta&space;t&space;\right&space;)&space;^{1}}{1!}e^{-\Lambda_{ij}&space;\Delta&space;t}&space;&plus;&space;\frac{\left&space;(&space;\Lambda_{ij}&space;\Delta&space;t&space;\right&space;)&space;^{2}}{2!}&space;e^{-\Lambda_{ij}&space;\Delta&space;t}&space;&plus;&space;..." title="p_{ij} \equiv \sum_{k=1}^{\infty } p_i(k;\Lambda_{ij}) = \frac{\left ( \Lambda_{ij} \Delta t \right ) ^{1}}{1!}e^{-\Lambda_{ij} \Delta t} + \frac{\left ( \Lambda_{ij} \Delta t \right ) ^{2}}{2!} e^{-\Lambda_{ij} \Delta t} + ..." />

Since Poisson distribution normalizes to unity, then the simple representation of the above equation can be written as,

<img src="https://latex.codecogs.com/svg.latex?p_{ij}&space;=&space;1&space;-&space;\frac{\left&space;(&space;\Lambda_{ij}&space;\Delta&space;t&space;\right&space;)&space;^{0}}{0!}&space;e^{-\Lambda_{ij&space;\Delta&space;t}}&space;=&space;1&space;-&space;e^{-\Lambda_{ij}&space;\Delta&space;t}" title="p_l = 1 - \frac{\left ( \Lambda_l \Delta t \right ) ^{0}}{0!} e^{-\Lambda_{il \Delta t}} = 1 - e^{-\Lambda_{il} \Delta t}" />

Also, the probability of zero-chance of event-<i>j</i> to occur is given by,

<img src="https://latex.codecogs.com/svg.latex?q_{ij}&space;\equiv&space;p(k=0;\Lambda_{ij})&space;=&space;\frac{\left&space;(&space;\Lambda_{ij}&space;\Delta&space;t&space;\right&space;)&space;^{0}}{0!}&space;e^{-\Lambda_{ij&space;\Delta&space;t}}&space;=&space;e^{-\Lambda_{ij}&space;\Delta&space;t}" title="q_l \equiv p(k=0;\Lambda_l) = \frac{\left ( \Lambda_l \Delta t \right ) ^{0}}{0!} e^{-\Lambda_{il \Delta t}} = e^{-\Lambda_{il} \Delta t}" />


<div align="justify">
So far, we have established the event probability of a single event. Recall that an isotope species can be associated with more than one transmutation events, <i>j</i>∈{0,1,2,…,<i>J<sub>i</sub></i>}. For simplicity, let E<sub>1</sub> be the condition where event-<i>l</i>, <i>l</i>∈{0,1,2,…,<i>J<sub>i</sub></i>}, is occurring, and E<sub>2</sub> be a condition where all events except removal event-<i>l</i> are not occurring. The joint probability of both E<sub>1</sub> and E<sub>2</sub> are true is given by,
</div>

\
<img src="https://latex.codecogs.com/svg.latex?P(E_1&space;\cap&space;E_2)&space;=&space;p_{il}&space;\times&space;\prod_{\forall&space;j&space;\ne&space;l}^{}&space;q_{ij}" title="P(E_1 \cap E_2) = p_l \times \prod_{\forall j \ne l}^{} q_j" />

By using the given definition of <i>p<sub>ij</sub></i> and <i>q<sub>ij</sub></i>, the un-normalized probability distribution can be expressed as,

<img src="https://latex.codecogs.com/svg.latex?\widetilde{\pi}_{il}&space;=&space;\left&space;(&space;1&space;-&space;e^{-\Lambda_{il}&space;\Delta&space;t}&space;\right&space;)\times&space;\prod_{\forall&space;l&space;\ne&space;j}^{}&space;e^{-\Lambda_{il}&space;\Delta&space;t}" title="\widetilde{\pi}_{il} = \left ( 1 - e^{-\Lambda_{il} \Delta t} \right )\times \prod_{\forall l \ne j}^{} e^{-\Lambda_{il} \Delta t}" />

The probability of no-removal, <i>l</i>=0 is given by the product of zero-chance probability of all possible events,

<img src="https://latex.codecogs.com/svg.latex?\widetilde{\pi}_{i0}&space;=&space;\prod_{\forall&space;l}^{}&space;e^{-\Lambda_{il}&space;\Delta&space;t}" title="\widetilde{\pi}_{il} = \prod_{\forall l}^{} e^{-\Lambda_l \Delta t}" />

Therefore, the normalized joint Poisson distribution can be formally written as,


<img src="https://latex.codecogs.com/svg.latex?\pi_{il}&space;\equiv&space;\pi_i(l;&space;\Lambda_1,&space;...,&space;\Lambda_{J_i})&space;=&space;c_i&space;\prod_{j=1}^{J_i}&space;\left&space;(&space;\delta_{jl}&space;&plus;&space;\left&space;(&space;-1&space;\right&space;)^{\delta_{jl}}&space;e^{-\Lambda_j&space;\Delta&space;t}&space;\right&space;)" title="\pi_{il} \equiv \pi_i(l; \Lambda_1, ..., \Lambda_{J_i}) = c_i \prod_{j=1}^{J_i} \left ( \delta_{jl} + \left ( -1 \right )^{\delta_{jl}} e^{-\Lambda_j \Delta t} \right )" />

where c<sub><i>i</i></sub> is the normalization constant and δ<sub><i>jl</i></sub> is the Kronecker delta. The value of c<sub><i>i</i></sub> is given by,


<img src="https://latex.codecogs.com/svg.latex?c_i&space;=&space;\left&space;(&space;\sum_{l=0}^{J_i}&space;\widetilde{\pi}_{il}&space;\right&space;)^{-1}" title="c_i = \left ( \sum_{l=0}^{J_i} \widetilde{\pi}_l \right )^{-1}" />

