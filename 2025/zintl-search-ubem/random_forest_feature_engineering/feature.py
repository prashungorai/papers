import numpy as np
from pylada.crystal import read, neighbors
import pylada.periodic_table as PT



# Covalent radii (in angstrom) are taken from doi.org/10.1021/jp5065819 | J. Phys. Chem. A 2015, 119, 2326−2337
covalent_radii = {'P': 1.11, 'As': 1.21, 'Sb': 1.40, 'Bi': 1.60, 'Si': 1.16, 'Ge': 1.21, 'Sn': 1.40, 'Pb': 1.44,
        'Al': 1.26, 'Ga': 1.24, 'In': 1.42, 'Tl': 1.44, 'Zn': 1.20, 'Cd': 1.44, 'Hg': 1.44, 'Mn': 1.19, 'Be': 1.02, 'Mg': 1.39}

principal_quantum_no = {'P': 3, 'As': 4, 'Sb': 5, 'Bi': 6, 'Si': 3, 'Ge': 4, 'Sn': 5, 'Pb': 6,
        'Al': 3, 'Ga': 4, 'In': 5, 'Tl': 6, 'Mn': 4 ,'Zn': 4, 'Cd': 5, 'Hg': 6, 'Be': 2, 'Mg': 3}

charge = {'Li': 1, 'Na': 1,  'K': 1, 'Rb': 1, 'Cs': 1, 'Be': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2, 'Mn': 2,'Zn': 2,
          'Cd': 2, 'Hg': 2, 'Al': 3, 'Ga': 3, 'In': 3, 'Tl': 1, 'Si': 2, 'Ge': 2, 'Sn': 2, 'Pb': 2, 'P': 3, 'As': 3,
          'Sb': 3, 'Bi': 3}

#Ionic radii (in angstrom) are taken from doi.org/10.1107/S0567739476001551 | Acta Cryst. (1976). A32, 751-767 
ionic_radii = {'Li': 0.9, 'Na': 1.16,  'K': 1.52, 'Rb': 1.66, 'Cs': 1.81, 'Be': 0.59, 'Mg': 0.86, 'Ca': 1.14, 'Sr': 1.32, 'Ba': 1.49, 'Mn': 0.97,'Zn': 0.88,
          'Cd': 1.09, 'Hg': 1.16, 'Al': 0.675, 'Ga': 0.76, 'In': 0.94, 'Tl': 1.64, 'Si': 0.54, 'Ge': 0.87, 'Sn': 0.83, 'Pb': 1.33, 'P': 0.58, 'As': 0.72,
          'Sb': 0.9, 'Bi': 1.17}

#Valence electrons used by the PAW  (Projected Augmented Wave) potentials
valence = {'Li': 3, 'Na': 7,  'K': 9, 'Rb': 9, 'Cs': 9, 'Be': 2, 'Mg': 2, 'Ca': 8, 'Sr': 10, 'Ba': 10, 'Mn': 7,'Zn': 12,
          'Cd': 12, 'Hg': 12, 'Al': 3, 'Ga': 13, 'In': 13, 'Tl': 13, 'Si': 4, 'Ge': 14, 'Sn': 14, 'Pb': 14, 'P': 5, 'As': 5,
          'Sb': 5, 'Bi': 15}

def get_atom_types_counts(structure):

    """
    Function to get counts of each atom type in a structure.
    
    Args:
        - structure: Pylada structure object.

    Returns:
        - atom_types_counts: Dictionary containing counts of each atom type.
    """
    atom_types_counts = {}
    for site in structure:
        atom_type = site.type
        atom_types_counts[atom_type] = atom_types_counts.get(atom_type, 0) + 1
    return atom_types_counts


def calculate_difference(structure,en_tol=0.90):
    
    """
    Function to group cations and anions. 
    It also calculates the difference in electronegativity between cations and anions in a structure.

    Args:
        - structure: Pylada structure object.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.

    Returns:
        - tuple: A tuple containing the following:
               - The difference in average electronegativity between anion and cation groups.
               - The difference in electronegativity between highest and lowest electronegative elements in the structure.
               - A dictionary containing counts of cations, where the keys are atom types.
               - A dictionary containing counts of anions, where the keys are atom types.
               - The total number of cations.
               - The total number of anions.
    """
    atom_types_counts = get_atom_types_counts(structure)

    # default is 0.90, user can specify a different tolerance 
    electronegativity_tolerance = en_tol

    cations = {}
    anions = {}

    # Identify cations and anions based on the original electronegativity
    max_en_atom_type = max(atom_types_counts, key=lambda x: getattr(PT, x).pauling)
    min_en_atom_type = min(atom_types_counts, key=lambda x: getattr(PT, x).pauling)

    EN_diff_max_min = getattr(PT, max_en_atom_type).pauling - getattr(PT, min_en_atom_type).pauling
    for atom_type, count in atom_types_counts.items():
        if atom_type == min_en_atom_type:
            cations[atom_type] = count
        elif atom_type == max_en_atom_type:
            anions[atom_type] = count

    remaining_atom_types = set(atom_types_counts.keys()) - {min_en_atom_type, max_en_atom_type}

    # Categorize remaining atom types based on electronegativity difference
    for atom_type in remaining_atom_types:
        en_difference = abs(getattr(PT, atom_type).pauling - getattr(PT, max_en_atom_type).pauling)
        if en_difference < electronegativity_tolerance:
            anions[atom_type] = atom_types_counts[atom_type]
        else:
            cations[atom_type] = atom_types_counts[atom_type]

    # Calculate the total number of cations and anions
    total_cations = sum(cations.values())
    total_anions = sum(anions.values())

    # Calculate the average electronegativity for cations and anions
    avg_cation_en = sum(getattr(PT, atom_type).pauling * count for atom_type, count in cations.items()) / total_cations
    avg_anion_en = sum(getattr(PT, atom_type).pauling * count for atom_type, count in anions.items()) / total_anions

    return avg_anion_en - avg_cation_en, EN_diff_max_min, cations, anions, total_cations, total_anions

def calculate_difference_avg(structure,en_tol=0.90):
    
    """
    Function to calculate the difference in average electronegativity between anion and cation groups.

    Args:
        - structure: Pylada structure object.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.
    Returns:
        - float: The difference in average electronegativity between anions and cations.
    """
    atom_types_counts = get_atom_types_counts(structure)
    difference, _, _, _, _, _ = calculate_difference(structure,en_tol)
    return difference


def calculate_EN_diff_max_min(structure,en_tol=0.90):
    
    """
    Function to calculate the difference in electronegativity between highest and lowest electronegative elements in the structure.

    Args:
        - structure: Pylada structure object.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.
    Returns:
        - float: The difference in electronegativity between highest and lowest electronegative elements in the structure.
    """
    atom_types_counts = get_atom_types_counts(structure)
    _, EN_diff_max_min, _, _, _, _ = calculate_difference(structure,en_tol)
    return EN_diff_max_min

def calculate_en_std_dev(structure):
    
    """
    Function to calculate the standard deviation of electronegativities for atoms in a given structure.

    Args:
        - structure: Pylada structure object.

    Returns:
        - float: The standard deviation of electronegativities.
    """
    atom_types_counts = get_atom_types_counts(structure)
    electronegativities = [getattr(PT, atom_type).pauling for atom_type in atom_types_counts]
    return np.std(electronegativities)


def calculate_average_electronegativity(structure):
    
    """
    Function to calculate the average of electronegativities for atoms in a given structure.

    Args:
        - structure: Pylada structure object.

    Returns:
        - float: The average of electronegativities.
    """
    atom_types_counts = get_atom_types_counts(structure)
    electronegativities = [getattr(PT, atom_type).pauling for atom_type in atom_types_counts]
    return np.mean(electronegativities)

def calculate_sum_ionic_radii(structure):
    
    """
    Function to calculate the sum of ionic radii for the most electronegative and least electronegative elements in a given structure.

    Args:
        - structure: Pylada structure object.
    Returns:
        - float: The sum of ionic radii for the most electronegative and least electronegative elements.
    """
    atom_types_counts = get_atom_types_counts(structure)
    max_en_element = max(atom_types_counts, key=lambda x: getattr(PT, x).pauling)
    min_en_element = min(atom_types_counts, key=lambda x: getattr(PT, x).pauling)
    max_ionic_radius = ionic_radii.get(max_en_element)
    min_ionic_radius = ionic_radii.get(min_en_element)

    assert max_ionic_radius is not None, f"Ionic radius not found for element: {max_en_element} in ionic_radii dictionary"
    assert min_ionic_radius is not None, f"Ionic radius not found for element: {min_en_element} in ionic_radii dictionary"
    
    return max_ionic_radius + min_ionic_radius



def calculate_weighted_avg_ionic_radii(structure):
    
    """
    Function to calculate the weighted average of ionic radii for all atom types in a given structure.

    Args:
        - structure: Pylada structure object.

    Returns:
        - float: The weighted average of ionic radii for all atom types in the structure.
    """
    atom_types_counts = get_atom_types_counts(structure)
    weighted_sum = 0
    total_count = 0
    
    for atom_type, count in atom_types_counts.items():
        ionic_radius = ionic_radii.get(atom_type)
        
        assert ionic_radius is not None, f"Ionic radius not found for element: {atom_type} in ionic_radii dictionary"
        weighted_sum += count*ionic_radius
        total_count += count
    
    return weighted_sum/total_count


def calculate_cation_charge_density(structure, en_tol=0.90):
    
    """
    Function to calculate the average charge density contributed by cations in a given structure.

    Args:
        - structure: Pylada structure object.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.
    Returns:
        - float: The average cation charge density per atom of the structure.

    """
    atom_types_counts = get_atom_types_counts(structure)
    _, _, cations, _, _, _ = calculate_difference(structure, en_tol)
    
    total_charge = 0
    total_atoms = 0
    
    for atom_type in cations:
        atom_count = atom_types_counts.get(atom_type, 0)
        atom_charge = charge.get(atom_type)
        
        assert atom_charge is not None, f"Charge not found for cation: {atom_type} in charge dictionary"
            
        total_charge += atom_count*atom_charge
        total_atoms += atom_count
    
    return total_charge/len(structure)

def calculate_charge_per_size_cations(structure, en_tol=0.90):
    
    """
    Function to calculate the charge per size for cations in a given structure.

    Args:
        - structure: Pylada structure object.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.
    Returns:
        - float: The average charge per size for cations.
    """
    atom_types_counts = get_atom_types_counts(structure)
    _, _, cations, _, _, _ = calculate_difference(structure, en_tol)
    
    num_cation_types = len(cations)
    sum_charge_size_product = 0
    
    for atom_type in cations:
        atom_charge = charge.get(atom_type)
        ionic_radius = ionic_radii.get(atom_type)
        
        assert atom_charge is not None, f"Charge not found for cation: {atom_type} in charge dictionary"
        
        assert ionic_radius is not None, f"Ionic radius not found for element: {atom_type} in ionic_radii dictionary"
        
        sum_charge_size_product += atom_charge/ionic_radius
    
    return sum_charge_size_product/num_cation_types


def calculate_total_count_per_type(structure, cov_radii_tol=0.1, en_tol=0.90):
    
    """
    Function to get the total count of each bond type in a given structure.

    Args:
        - structure: Pylada structure object.
        - cov_radii_tol (float, optional): The tolerance for covalent radii.
          Defaults to 0.1.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.

    Returns:
        - dictionary: A dictionary containing the total count of each bond type.
    """
    atom_types_counts = get_atom_types_counts(structure)
    _, _, _, anions, _, _ = calculate_difference(structure, en_tol)
    covalent_radii_tolerance = cov_radii_tol
    bond_count = {}
    processed_bonds = set()

    for i in range(len(structure)):
        if structure[i].type in anions:
            ngh = neighbors(structure, 2, structure[i], 0.3)
            for j in range(len(ngh)):
                bond_length = ngh[j][-1]
                for k in range(len(structure)):
                    if structure[k].type in anions:  
                        radii_1 = covalent_radii.get(structure[i].type)
                        radii_2 = covalent_radii.get(structure[k].type)
                    
                        assert radii_1 is not None, f"Covalent radius not found for atom type: {structure[i].type} in Covalent_radii dictionary"
                        
                        assert radii_2 is not None, f"Covalent radius not found for atom type: {structure[k].type} in Covalent_radii dictionary"
                    
                        radii_1 += covalent_radii_tolerance
                        radii_2 += covalent_radii_tolerance

                        if i != k and ngh[j][0] == structure[k] and bond_length < (radii_1 + radii_2):
                            bond_pair = (structure[i].type, structure[k].type, i,k)
                            reverse_bond_pair = (structure[k].type, structure[i].type, k, i)
                            if bond_pair not in processed_bonds and reverse_bond_pair not in processed_bonds:
                                bond_count[bond_pair] = bond_count.get(bond_pair, 0) + 1
                                processed_bonds.add(bond_pair)


    total_count_per_type = {}
    for bond_type, count in bond_count.items():
        bond_type_without_indices = tuple(sorted([bond_type[0], bond_type[1]]))
        total_count_per_type[bond_type_without_indices] = total_count_per_type.get(bond_type_without_indices, 0) + 1

    return total_count_per_type


def calculate_total_bond_density(structure,cov_radii_tol=0.1,en_tol=0.90):
    
    """
    Function to calculate the total covelant bond density for a given structure.

    Args:
        - structure: Pylada structure object.
        - cov_radii_tol (float, optional): The tolerance for covalent radii.
          Defaults to 0.1.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.
    
    Returns:
        - tuple: A tuple containing the total bond density and a dictionary of bond densities for each bond type.
    """
    
    _, _, _, anions, _, _ = calculate_difference(structure,en_tol)
    total_count_per_type = calculate_total_count_per_type(structure,cov_radii_tol,en_tol)

    assert total_count_per_type is not None, "Unable to calculate total bond density"

    total_bond_density = 0
    bond_densities = {}

    for bond_type, count in total_count_per_type.items():
        anion1, anion2 = bond_type
        anion1_count = anions.get(anion1, 0)
        anion2_count = anions.get(anion2, 0)
        bond_type_sorted = tuple(sorted([anion1, anion2]))  # Introducing bond_type_sorted

        if anion1 == anion2:
            bond_densities[bond_type_sorted] = count/anion1_count if anion1_count > 0 else 0
            total_bond_density += bond_densities[bond_type_sorted]
        else:
            total_anion_count = anion1_count + anion2_count
            bond_densities[bond_type_sorted] = count/total_anion_count if total_anion_count > 0 else 0
            total_bond_density += bond_densities[bond_type_sorted]

    return total_bond_density, bond_densities

def calculate_orbital_energy(Z, n):
    
    """
    Function to calculate the orbital energy of an element based on the Bohr model.

    Args:
        - Z (int): Atomic number of the element.
        - n (int): Principal quantum number of the electron.

    Returns:
        - float: The calculated orbital energy of an element.
    """
    return (Z**2)*(1/n**2)

def calculate_orbital_energies(structure, cov_radii_tol=0.1, en_tol=0.90):
    
    """
    Function to calculate the orbital energy parameter for a bond in a given structure.

    Args:
        - structure: Pylada structure object.
        - cov_radii_tol (float, optional): The tolerance for covalent radii.
          Defaults to 0.1.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.
    Returns:
        - tuple: A tuple containing a dictionary of orbital energy parameters for each bond type 
                 and a float representing the total orbital energy parameter corresponding to all bond types.
    """
    total_count_per_type = calculate_total_count_per_type(structure, cov_radii_tol, en_tol)
    assert total_count_per_type is not None, "Unable to calculate orbital energy parameter"
    orbital_energies = {}
    Total_Orbital_E_Diff = 0
    for bond_type, count in total_count_per_type.items():
        anion1, anion2 = bond_type
        if anion1 in principal_quantum_no and anion2 in principal_quantum_no:
            n1 = principal_quantum_no[anion1]
            n2 = principal_quantum_no[anion2]
            Z1 = getattr(PT, anion1).atomic_number
            Z2 = getattr(PT, anion2).atomic_number
            orbital_energy1 = calculate_orbital_energy(Z1, n1)
            orbital_energy2 = calculate_orbital_energy(Z2, n2)
            if anion1 == anion2:
                energy_difference = 1/n1
            else:
                energy_difference = 1/abs(orbital_energy1 - orbital_energy2)
            bond_type_sorted = tuple(sorted([anion1, anion2]))
            orbital_energies[bond_type_sorted] = energy_difference
            Total_Orbital_E_Diff += energy_difference
        else:
            assert False, f"Principal quantum number not found for {anion1} or {anion2} in principal_quantum_no dictionary"
    return orbital_energies, Total_Orbital_E_Diff


def calculate_EN_Diff_dict(structure,cov_radii_tol=0.1,en_tol=0.90):
    
    """
    Function to calculate the electronegativity difference parameter for each bond type in a given structure.

    Args:
        - structure: Pylada structure object.
        - cov_radii_tol (float, optional): The tolerance for covalent radii.
          Defaults to 0.1.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.
    Returns:
        - tuple: A tuple containing a dictionary of electronegativity difference parameters for each bond type 
                 and the total electronegativity difference parameter corresponding to all bond types.
    """
    total_count_per_type = calculate_total_count_per_type(structure,cov_radii_tol,en_tol)
    assert total_count_per_type is not None, "Unable to calculate electronegativity parameter"
    EN_Diff_dict = {}
    Total_EN_Diff = 0
    for bond_type in total_count_per_type.keys():
        anion1, anion2 = bond_type
        bond_type_sorted = tuple(sorted([anion1, anion2]))
        if anion1 == anion2:
            EN_Diff = 1 - 0
        else:
            electroneg_diff = abs(getattr(PT, anion1).pauling - getattr(PT, anion2).pauling)
            EN_Diff = 1 - electroneg_diff
        EN_Diff_dict[bond_type_sorted] = EN_Diff
        Total_EN_Diff += EN_Diff
    return EN_Diff_dict, Total_EN_Diff

def calculate_weighted_bond_densities(structure,cov_radii_tol=0.1,en_tol=0.90):
    
    """
    Function to calculate the weighted bond densities for each bond type in the structure, 
    taking into account the electronegativity difference and orbital energy difference parameters.

    Args:
        - structure: Pylada structure object.
        - cov_radii_tol (float, optional): The tolerance for covalent radii.
          Defaults to 0.1.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.
    
    Returns:
        - float: The total covalent bond density weighted by electronegativity 
                 and orbital energy parameters.
    """
    _, bond_densities = calculate_total_bond_density(structure,cov_radii_tol,en_tol)
    assert bond_densities is not None, "Unable to calculate bond density"
    
    EN_Diff_dict, _ = calculate_EN_Diff_dict(structure,cov_radii_tol,en_tol)
    assert EN_Diff_dict is not None, "Unable to calculate electronegativity difference parameter"
    
    orbital_energies, _ = calculate_orbital_energies(structure,cov_radii_tol,en_tol)
    assert orbital_energies is not None, "Unable to calculate orbital energy parameter"

    weighted_bond_densities = {}
    weighted_total_bond_density = 0

    for bond_type_sorted, EN_Diff in EN_Diff_dict.items():
        energy_difference = orbital_energies.get(bond_type_sorted, 0)
        bond_density = bond_densities.get(bond_type_sorted, 0)

        weighted_bond_density = bond_density * EN_Diff * energy_difference
        weighted_bond_densities[bond_type_sorted] = weighted_bond_density
        weighted_total_bond_density += weighted_bond_density
    return weighted_total_bond_density

def calculate_weighted_bond_densities_2(structure,cov_radii_tol=0.1,en_tol=0.90):
    
    """
    Function to calculate the weighted bond densities for each bond type in the structure, 
    taking into account the electronegativity difference parameter alone.

    Args:
        - structure: Pylada structure object.
        - cov_radii_tol (float, optional): The tolerance for covalent radii.
          Defaults to 0.1.
        - en_tol (float, optional): The tolerance for considering atoms as cations or anions based on electronegativity.
          Defaults to 0.90.
    
    Returns:
        - float: Total covalent bond density weighted by electronegativity parameter.
    """
    _, bond_densities = calculate_total_bond_density(structure,cov_radii_tol,en_tol)
    assert bond_densities is not None, "Unable to calculate bond density"

    EN_Diff_dict, _ = calculate_EN_Diff_dict(structure,cov_radii_tol,en_tol)
    assert EN_Diff_dict is not None, "Unable to calculate electronegativity difference parameter"

    weighted_bond_densities_2 = {}
    weighted_total_bond_density_2 = 0

    for bond_type_sorted, EN_Diff in EN_Diff_dict.items():
        bond_density = bond_densities.get(bond_type_sorted, 0)

        weighted_bond_density_2 = bond_density * EN_Diff
        weighted_bond_densities_2[bond_type_sorted] = weighted_bond_density_2
        weighted_total_bond_density_2 += weighted_bond_density_2

    return weighted_total_bond_density_2

def calculate_weighted_average_valence(structure):
    
    """
    Function to calculate the weighted average valence of atoms in a given structure.

    Args:
        - structure: Pylada structure object.

    Returns:
        - float: The weighted average valence of atoms considered by PAW potentials in the structure.
    """
    atom_types_counts = get_atom_types_counts(structure)
    total_atoms = sum(atom_types_counts.values())
    for atom_type in atom_types_counts:
        atom_valence = valence.get(atom_type)
    assert atom_valence is not None, f"Valence not found for atom: {atom_type} in valence dictionary"

    weighted_sum = sum(atom_types_counts[atom_type]*valence[atom_type] for atom_type in atom_types_counts)
    weighted_average = weighted_sum/total_atoms
    return weighted_average


def calculate_atomic_mass(structure):
    
    """
    Function to calculate the atomic mass of a given structure in atomic units.

    Args:
        - structure: Pylada structure object.

    Returns:
        - float: The total atomic mass of the structure.
    """
    atom_types_counts = get_atom_types_counts(structure)
    atomic_mass = sum(getattr(PT, atom_type).atomic_weight*atom_types_counts[atom_type] for atom_type in atom_types_counts)/len(structure)
    return atomic_mass

def calculate_atomic_density(structure):
    
    """
    Function to calculate the atomic density of a given structure in atomic units per Angstrom cubed.

    Args:
        - structure: Pylada structure object.

    Returns:
        - float: The atomic density of the structure.
    """
    atom_types_counts = get_atom_types_counts(structure)
    atomic_mass = calculate_atomic_mass(structure)
    atomic_volume = structure.volume/len(structure)
    atomic_density = atomic_mass/atomic_volume
    return atomic_density


