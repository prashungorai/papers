# VASP inputs and relaxed defected supercells

NOTE: 
1. Gamma-point was used. An example VASP input (INCAR\_relax) is provided. It takes V<sub>Al</sub><sup>-3</sup> as an example.

2. Within each composition there are three folders:
* bulk/ contains the undefected SQS supercell
* defect\_nominal\_charge\_all\_sites/ contains the relaxed structures for defects only in nominal their charge (e.g., V<sub>Al</sub><sup>-3</sup>, V<sub>Sc</sub><sup>-3</sup>, and V<sub>N</sub><sup>+3</sup>) but on all the unique lattice sites.   
* defect\_upper\_lower\_energy\_sites/ contains the relaxed structures for defects in different charge states but only on the lattice sites that have highest and lowest energies when in nominal charge states.
