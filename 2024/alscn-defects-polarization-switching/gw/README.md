# Steps to run GW calculations with fixed GGA(+U) wavefunctions

The calculations were done using VASP 5.4.4 and we generally followed the steps described in [VASP Wiki](https://www.vasp.at/wiki/index.php/Practical_guide_to_GW_calculations).
There are three steps and the associated VASP inputs are provided as INCAR\_3, INCAR\_2, and INCAR\_3 

1. INCAR\_1
* a standard ground-state calculation with few unoccupied orbitals only 

2. INCAR\_2  and  WAVECAR froms step 1 
* a calculation of a large number of unoccupied orbitals and have the flag LOPTICS=.TRUE. to write WAVEDER

3. INCAR\_3 with WAVECAR and WAVEDER from step 2
* GW calculations (fixed wavefunction and update eigenvalues only) 
* repeat till convergence or set larger NELM

