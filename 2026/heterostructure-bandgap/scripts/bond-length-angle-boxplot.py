import pandas as pd
from glob import iglob
from pylada.crystal import read, neighbors
import numpy as np
import os

anions = ['N', 'S', 'Se', 'O', 'I', 'C', 'Te']

def myneighbors(strc=None, howmany=12, pos=None, tol=0.2):
    nghs = []
    for ng in neighbors(strc, howmany, pos, 0.1):
        if len(nghs) == 0:
            nghs.append(ng)
            continue
        elif nghs[0][2] <= ng[2] <= (1. + tol) * nghs[0][2]:
            nghs.append(ng)
    return nghs

def parse_folder_name(name):
    folder = name.split('/')[-1]
    parts  = folder.split('_')
    return {
        'base_strc': parts[-3],
        'alloy':     parts[2],
        'conc':      parts[0],
        'version':   parts[-1],
    }

def calc_angle(v1, v2):
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))

def analyse_structure(strc, meta, tol=0.25):
    bond_rows  = []
    angle_rows = []
    cation_no  = 0

    for atom in strc:
        if atom.type in anions:
            continue

        cation_no += 1
        nghs = myneighbors(strc, 4, atom.pos, tol)

        # ng[1] is already the bond vector from cation to neighbour
        # ng[2] is the bond length
        # only keep anion neighbours
        anion_nghs = [ng for ng in nghs if ng[0].type in anions]

        # ── bond lengths ─────────────────────────────────────────────────
        for i, ng in enumerate(anion_nghs):
            bond_rows.append({
                **meta,
                'cation':        atom.type,
                'cation_no':     cation_no,
                'anion':         ng[0].type,
                'bond_label':    f"{atom.type}-{ng[0].type}",
                'bond_index':    i + 1,
                'bond_length_A': round(ng[2], 6),
            })

        # ── bond angles: unique N-Cation-N pairs ─────────────────────────
        for i in range(len(anion_nghs)):
            for j in range(i + 1, len(anion_nghs)):
                v1    = anion_nghs[i][1]    # bond vector to 1st anion
                v2    = anion_nghs[j][1]    # bond vector to 2nd anion
                angle = calc_angle(v1, v2)
                angle_rows.append({
                    **meta,
                    'cation':         atom.type,
                    'cation_no':      cation_no,
                    'angle_label':    f"{anion_nghs[i][0].type}-{atom.type}-{anion_nghs[j][0].type}",
                    'angle_pair':     f"{i+1}-{j+1}",
                    'bond_angle_deg': round(angle, 4),
                })

    return bond_rows, angle_rows

def compile_analysis(glob_pattern, tol=0.25):
    all_bond_rows  = []
    all_angle_rows = []

    for name in sorted(iglob(glob_pattern)):
        if not os.path.isdir(name):
            continue

        poscar = os.path.join(name, 'CONTCAR')
        if not os.path.isfile(poscar):
            print(f"  [SKIP] no CONTCAR in {name}")
            continue

        meta = parse_folder_name(name)
        strc = read.poscar(poscar)

        bond_rows, angle_rows = analyse_structure(strc, meta, tol)
        all_bond_rows.extend(bond_rows)
        all_angle_rows.extend(angle_rows)
        print(f"Processed: {name.split('/')[-1]}  "
              f"({len(bond_rows)} bonds, {len(angle_rows)} angles)")

    if not all_bond_rows:
        print("No data found.")
        return None, None

    return pd.DataFrame(all_bond_rows), pd.DataFrame(all_angle_rows)

# ── main ─────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    GLOB_PATTERN = 'GaN_w/Ga_Al_GaN_w/relaxed_polar/*'
    OUTPUT_DIR   = 'GaN_w/Ga_Al_GaN_w/relaxed_polar/'

    bond_df, angle_df = compile_analysis(GLOB_PATTERN, tol=0.25)

    if bond_df is not None:
        bond_path  = os.path.join(OUTPUT_DIR, 'bond_length_analysis.csv')
        angle_path = os.path.join(OUTPUT_DIR, 'bond_angle_analysis.csv')

        bond_df.to_csv(bond_path,   index=False)
        angle_df.to_csv(angle_path, index=False)

        print(f"\nSaved: {bond_path}   ({len(bond_df)} rows)")
        print(f"Saved: {angle_path}  ({len(angle_df)} rows)")

        print("\n── Bond Lengths ──")
        print(bond_df.to_string(index=False))
        print("\n── Bond Angles ──")
        print(angle_df.to_string(index=False))
