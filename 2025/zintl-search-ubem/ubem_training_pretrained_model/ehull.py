
from typing import List
import pandas as pd

from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element


# ---FERE reference chemical potentials------
fere_chempot = {'Ag': -0.83, 'Al': -3.02, 'As': -5.06, 'Au': -2.23, 'Ba': -1.39, 'Be': -3.40, 'Bi': -4.39, 'Ca': -1.64,
                'Cd': -0.56, 'Cl': -1.63, 'Co': -4.75, 'Cr': -7.22, 'Cu': -1.97, 'F': -1.70, 'Fe': -6.15, 'Ga': -2.37,
                'Ge': -4.14, 'Hf': -7.40, 'Hg': -0.12, 'In': -2.31, 'Ir': -5.96, 'K': -0.80, 'La': -3.66, 'Li': -1.65,
                'Mg': -0.99, 'Mn': -7.00, 'N': -8.51, 'Na': -1.06, 'Nb': -6.69, 'Ni': -3.57, 'O': -4.76, 'P': -5.64,
                'Pd': -3.12, 'Pt': -3.95, 'Rb': -0.68, 'Rh': -4.76, 'S': -4.00, 'Sb': -4.29, 'Sc': -4.63, 'Se': -3.55,
                'Si': -4.99, 'Sn': -3.79, 'Sr': -1.17, 'Ta': -8.82, 'Te': -3.25, 'Ti': -5.52, 'V': -6.42, 'Y': -4.81,
                'Zn': -0.84, 'Zr': -5.87}

ferev2_chempot = {'Ag': -0.79, 'Al': -3.27, 'Al_anion': -3.55, 'As': -4.95, 'As_cation': -4.42, 'Au': -1.96, 'B': -6.73,
                  'B_anion': -6.44, 'Ba': -1.44, 'Be': -3.50, 'Bi_cation': -4.22, 'Bi': -4.19, 'Br': -1.54, 'C': -8.94,
                  'Ca': -1.78, 'Cd': -0.64, 'Cl': -1.74, 'Co': -4.67, 'Cr': -7.08, 'Cu': -1.87, 'F': -1.44, 'Fe': -6.00,
                  'Ga': -2.53, 'Ge': -4.34, 'Ge_anion': -4.84, 'Hf': -7.38, 'Hg': -0.10, 'I': -1.53, 'In': -2.39,
                  'Ir': -6.31, 'K': -0.79, 'La': -3.76, 'Li': -1.58, 'Mg': -1.23, 'Mo': -7.37, 'Mn': -6.86, 'N': -8.46,
                  'Na': -1.02, 'Nb': -6.92, 'Ni': -3.65, 'O': -4.80, 'P': -5.17, 'P_cation': -5.14, 'Pb': -3.85,
                  'Pb_anion': -4.29, 'Pd': -3.00, 'Pt': -3.88, 'Rb': -0.58, 'Rh': -4.66, 'Ru': -6.14, 'S': -4.01,
                  'Sb_cation': -4.13, 'Sb': -4.16, 'Sc': -4.42, 'Se': -3.54, 'Si': -5.30, 'Si_anion': -5.40,
                  'Sn': -3.87, 'Sn_anion': -3.71, 'Sr': -1.32, 'Ta': -8.82, 'Te': -3.18, 'Te_cation': -2.75,
                  'Ti': -5.39, 'V': -6.40, 'W': -9.61, 'Y': -4.96, 'Zn': -0.94, 'Zr': -6.39,'Cs': -0.96, 'Tl':-1.938}
fere_entries = [PDEntry(e, energy) for e, energy in ferev2_chempot.items() if "ion" not in e]


def convex_hull_stability(comp: Composition,
                          predicted_energy: float,
                          competing_phases: List[PDEntry],
                          ):
    """ compute the decomposition energy for a given composition 
    :param comp: composition such as Li1Sc1F4
    :param predicted_energy: predicted eV/atom for the structure corresponding to this composition
    :param competing_phases: list of competing phases used to 
        construct the convex hull for the elements of the given composition
    """
    # first make sure the composition is in the reduced form
    comp = comp.reduced_composition
    # convert the predicted energyperatom to total energy
    energy = predicted_energy * comp.num_atoms
    entry = PDEntry(comp, energy)
    elements = set(comp.elements)

    # get the entries that have the same or a subset of elements.
    # Since we're including DFT relaxed energies in the competing phases,
    # also don't include entries that match the composition and predicted energy exactly
    # since we are trying to compute the decomposition energy for that structure
    curr_entries = [e for e in competing_phases
                    if len(set(e.composition.elements) - elements) == 0
                    and e != entry]
    phase_diagram = PhaseDiagram(curr_entries, elements=elements)
    #print(phase_diagram.qhull_entries)
    #print(phase_diagram.all_entries)
    decomp, decomp_energy = phase_diagram.get_decomp_and_phase_separation_energy(
        entry,  
        # docs say: "if you have a huge proportion of unstable entries,
        # then this check can slow things down."
        stable_only=True)
    # if this structure is not stable, then skip the stability range calculation

    # since the decomposition energy is negative,
    # add the current entry to the hull to get the electrochem_stability_window
    #print(decomp)
    if decomp_energy > 0:
        return decomp_energy, None

    # since the decomposition energy is negative,
    # add the current entry to the hull to get the electrochem_stability_window
    phase_diagram_w_comp = PhaseDiagram(phase_diagram.entries + [entry],
                                        elements=comp.elements)
    # return the decomposition energy
    return decomp_energy


def setup_competing_phases(competing_phases_files):
    if not isinstance(competing_phases_files, list):
        competing_phases_files = [competing_phases_files]
    all_competing_phases = [load_competing_phases(f) for f in competing_phases_files]

    # also add the individual elements
    competing_phases = pd.concat([pd.Series(fere_entries)] + all_competing_phases).reset_index()[0]
    return competing_phases


def load_competing_phases(competing_phases_file):
    """ Convert a file with the energy of compositions 
    into a list of pymatgen phase diagram entries.
    
    Can load either the NREL MatDB file, or a file with predicted energies
    """
    print(f"Reading {competing_phases_file}")
    df = pd.read_csv(competing_phases_file)
    print(f"\t{len(df)} lines")
    print(df.head(2))

    assert ('sortedformula' in df.columns or 'comp' in df.columns) \
        and ('energyperatom' in df.columns or 'predicted_energy' in df.columns)
    if 'sortedformula' not in df.columns:
        df.rename(columns={'comp': 'sortedformula'}, inplace=True)
    if 'energyperatom' not in df.columns:
        df.rename(columns={'predicted_energy': 'energyperatom'}, inplace=True)
    # print("columns after renaming:", df.columns)

    df['energy'] = (
        df.energyperatom *
        df.sortedformula.apply(lambda x: Composition(x).num_atoms)
    )
    # convert the dataframe to a list of PDEntries used to create the convex hull
    pd_entries = df.apply(
        lambda row: PDEntry(Composition(row.sortedformula),
                            row.energy),
        axis=1
    )
    print(f"\t{len(pd_entries)} entries")
    return pd_entries

