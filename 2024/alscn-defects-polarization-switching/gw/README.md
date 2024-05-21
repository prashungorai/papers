# Steps to run GW calculations using fixed GGA(+U) waefunctions

The calculations were done using VASP 5.4.4 and we generally followed the steps described in [VASP Wiki](https://www.vasp.at/wiki/index.php/Practical_guide_to_GW_calculations).
There are three steps and the associated VASP inputs are provided as INCAR<sub>1</sub>, INCAR<sub>2</sub>, and INCAR<sub>3</sub> 

1. INCAR<sub>1</sub>

* a standard ground-state calculation with few unoccupied orbitals only 

2. INCAR<sub>2</sub>  and  WAVECAR froms step 1 

* a calculation of a large number of unoccupied orbitals and have the flag LOPTICS=.TRUE. to write WAVEDER

3. INCAR<sub>3</sub> with WAVECAR and WAVEDER from step2

* GW calculations (fixed wavefunction and update eigenvalues only) 
* repeat till convergence or set larger NELM

