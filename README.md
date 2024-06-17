# PALT
**PALT** (<ins>P</ins>eriodic-inversion <ins>_a_</ins>_b_ _initio_ <ins>L</ins>aser <ins>T</ins>heory) solves Maxwell-Bloch equations for [limit cycles](https://en.wikipedia.org/wiki/Limit_cycle) with dynamic inversion, as well as linear cavity modes and fixed points. This program implement PALT for one-dimentional laser systems under outgoing boundary conditon. It is recommended to run on MATLAB R2020a or a later version.

The [PALT program](https://github.com/complexphoton/PALT/tree/main/PALT) includes the following components,
 - Four executable scripts,
   - [Dth0_Dth1.m](https://github.com/complexphoton/PALT/blob/main/PALT/Dth0_Dth1.m)
   - [Dth1_Dth2.m](https://github.com/complexphoton/PALT/blob/main/PALT/Dth1_Dth2.m)
   - [Dth2_Dth3.m](https://github.com/complexphoton/PALT/blob/main/PALT/Dth2_Dth3.m)
   - [FieldProfile.m](https://github.com/complexphoton/PALT/blob/main/PALT/Dth2_Dth3.m)
-	functions folder containing the core equations of PALT solver,
-	data folder that stores the inputs and outputs of the program.
  
Users can edit the four executable scripts to build the laser system and to generate the PALT data as presented in the manuscript. It is not recommended to edit any subfunction in the “functions” folder unless the computation method or the fundamental equations need to be revised.

## Dth0_Dth1.m
This script has two blocks.
The first block is for constructing the laser system. Users can build a multi-layer structure with one layer being active. The active layer is defined with gain parameters ($\omega_{ab}$ , $\gamma_{\parallel}$ , $\gamma_{\perp}$) and a linear refractive index $n_a$. Other passive layers are defined by layer thickness, refractive index, and absorption rate. All the parameters are normalized.
 - Spatial coordinates are normalized by the thickness of the active layer $L$.
 - Frequency/rate quantities, such as $\omega_{ab}$, $\gamma_{\parallel}$  and $\gamma_{\perp}$, are normalized by $c/L$, where $c$ is the speed of light. 
 - The absorption rate in physical unit is related to the “s” variable in the script via σ⁄ε_0 =(cn_p^2/L)s, where n_p is the refractive index of the passive cavity.

The default set-up is the EP laser in Fig.3a of the manuscript. 

The second block solves for the complex eigen frequencies of the system with outgoing boundaries. Gain saturation is not considered. Therefore, the solution is valid until one eigen frequency on the complex-frequency plane reaches the real axis from below. This fact is used by the solver to determine the first lasing threshold, $D_1^{th}$.

The eigenvalue solver requires initial guess (line 85). Different initial guess leads to different solutions. It is recommended to scan the function value of M0.m on the complex-frequency plane for the estimation of initial guess. This can be done by enabling line 114-128, which by default is disabled to save computation time. 

<img src= "https://github.com/complexphoton/PALT/assets/172996975/71ca0f91-b0c2-4a72-97d8-deb831ce56e4" width="400">

The solver exports the first threshold Dth1 and the eigenvalues of the system at Dth1. Please double check on the above plot that (a) Dth1_w(1) has negligible imaginary part, (b) all the other eigenvalues are below Dth1_w(1). If either of the two conditions is unsatisfied, adjust Dmax or the initial guess and re-run the script.

## Dth1_Dth2.m
This script calculates the lasing field and lasing frequency above the first threshold. It also solves the PALT perturbation equations for complex ω_d. The second threshold is then found as the Dmax value at which Im(ω_d )=0. The script also exports the lasing field, the lasing frequency ω_0, the line spacing ω_d and the corresponding eigenvectors. These data will be taken in Dth2_Dth3.m as the initial guess to solve for limit cycles slightly above the second threshold.

<img src= "https://github.com/complexphoton/PALT/assets/172996975/be58ca19-0101-48c1-baf1-57db8f8f6d41" width="400">

## Dth2_Dth3.m
This script calculates limit cycles, or frequency combs with periodic inversion above the second threshold. An appropriate initial guess is the key for the solver to converge correctly. The script provides two methods to form the initial guess.
  1. For $D_{max} \gtrsim D_2^{th}$,  use the data from Dth1_Dth2.m to form the initial guess. Use this method on two different $D_{max}$ values near $D_2^{th}$, then switch to the second method.
  2. If there are existing limit cycle solutions saved in the data folder, use them as the initial guess for a nearby point in the phase space.
	
The data folder has a sample of limit cycles. It is the EP comb at $D_{max}=0.2$, $\sigma⁄\varepsilon_0 =0.0041(cn_c^2/L)$, shown in Fig.3g of the manuscript. Using this comb as an initial guess, then run the script for Dmax = 0.21. The solution will be found within about 15~20 minutes, depending on the performance of hardware. A plot of the comb spectrum will be generated,

<img src= "https://github.com/complexphoton/PALT/assets/172996975/0a70b414-88c3-49b7-9a2b-8af684ff79ab" width="400">

## FieldProfile.m
Taking an existing limit cycle solution, this script calculates the electrical field, polarization, and population inversion of each comb line at any location using Green’s function interpolation. It can be used to plot the spatial profile of the field or other variables both in and out of the gain cavity. By default, the script reads the EP-comb sample, then plots the intensities of the first three comb lines as below, 

<img src= "https://github.com/complexphoton/PALT/assets/172996975/4a9cfcc0-a170-4945-810a-8dbaf7010ac0" width="450">

# FDTD simulation
The [FDTD program](https://github.com/complexphoton/PALT/tree/main/FDTD) iterates Maxwell-Bloch equations with [FDTD](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) method. It is used to demonstrate PALT. The folder includes,
- [**FDTD.jl**](https://github.com/complexphoton/PALT/blob/main/FDTD/FDTD.jl), the FDTD simulator.  
- **an example of initial condition**: xB.mat, xD.mat, xE.mat, xP_i.mat, xP_r.mat,
- [**DataProcessing.m**](https://github.com/complexphoton/PALT/blob/main/FDTD/DataProcessing.m). This MATLAB script is used to read the transient simulation results and to compute the spectrum using FFT.

## FDTD.jl
The FDTD simulator is developed in Julia. There are two build-in methods of initialization in the script (line 99 – line 117). Users can customize the initial condition as well. 

Here is the list of the outputs generated at the end of the simulation,
-	Variables’ values at each spatial pixel at the last time step, namely electrical field xE, magnetic field xB, population inversion xD, the real part of polarization xP_r, and the imaginary part of polarization xP_i.
-	The time-domain electrical field E, and population inversion D at a customized position.

All the outputs are in (*.mat) format. 

## Running the simulation
In FDTD.jl, the default setup is the EP laser presented in Fig.3a of the manuscript. The space dimension is in the unit of 1 nm. The time is in the unit of c/(1 nm) with c being the light speed. The spatial resolution is 4/nm. Total simulation time spans about 5 periods of the limit cycle.

As a test, the simulation result (xE, xB, xD, xP_i, xP_r) at a smaller pumping strength is used as the initial condition. The PALT solution for this EP laser is provided for comparison.
1. Run FDTD.jl in Julia. This simulation takes 20~30 hours, depending on the hardware performance.
2. Run DataProcessing.m in MATLAB. It generates the following plots,

<img src= "https://github.com/complexphoton/PALT/assets/172996975/5d2819cc-1b9f-4bb6-8030-456deb1bdb64" width="450">

<img src= "https://github.com/complexphoton/PALT/assets/172996975/c82fc7fe-2a27-4e95-9d4d-d732da568c58" width="450">
 

