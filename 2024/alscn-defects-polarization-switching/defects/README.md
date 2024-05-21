# VASP inputs and relaxed defected supercells

NOTE: 
1. Gamma-point was used. An example VASP input (INCAR\_relax)is provided using V<sub>Al</sub><sup>-3</sup> as an example.

2. Within each composition there are three folders:
* bulk/ contains the undefected structure
* defect\_nominal\_charge\_all\_sites/ contains the relaxed structures for defects in nominal charge (e.g., V<sub>Al</sub><sup>-3</sup>, V<sub>Al</sub><sup>-3</sup>, and V<sub>Al</sub><sup>-3</sup>) on all the uniqe lattice sites.   
* defect\_upper\_lower\_energy\_sites/ contains the relaxed structures for defects in different charge state but on the lattice sites that have highest and lowest energies when in nominal charge state.
