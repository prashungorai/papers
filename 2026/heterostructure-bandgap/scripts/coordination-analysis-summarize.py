import pandas as pd
from glob import iglob
import os

def parse_folder_name(name):
    folder = name.split('/')[-1]
    parts  = folder.split('_')
    return {
        'base_strc': parts[-3],
        'alloy':     parts[2],
        'conc':      parts[0],
        'version':   parts[-1],
    }

def load_coordination_csv(folder_path):
    csv_path = os.path.join(folder_path, 'coordination_analysis.csv')
    if not os.path.isfile(csv_path):
        print(f"  [SKIP] no coordination_analysis.csv in {folder_path}")
        return None
    return pd.read_csv(csv_path)

def add_coordination_percent(df):
    """
    For each row, coordination_percent = count of atoms with same species AND same CN
    divided by total cations in the folder, as a percentage.
    """
    total_cations = len(df)

    # count how many atoms share the same (cation, coordination_number)
    cn_counts = df.groupby(['cation', 'coordination_number'])['cation_no'].transform('count')

    df['coordination_percent'] = round(100 * cn_counts / total_cations, 2)
    return df

def compile_all_folders(glob_pattern):
    all_rows = []

    for name in sorted(iglob(glob_pattern)):
        if not os.path.isdir(name):
            continue

        df = load_coordination_csv(name)
        if df is None:
            continue

        meta = parse_folder_name(name)
        print(f"Processing: {name.split('/')[-1]}  "
              f"({len(df)} cations, conc={meta['conc']}, {meta['version']})")

        df = add_coordination_percent(df)

        for key, val in reversed(meta.items()):
            df.insert(0, key, val)

        all_rows.append(df)

    if not all_rows:
        print("No data found. Check your glob pattern.")
        return None

    summary_df = pd.concat(all_rows, ignore_index=True)
    summary_df = summary_df[['base_strc', 'alloy', 'conc', 'version',
                              'cation', 'cation_no', 'coordination_number',
                              'coordination_percent']]
    return summary_df

# ── main ─────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    GLOB_PATTERN = 'GaN_w/band_gap_bowing/AlScN_HSE/static_HSE_222/0.*'
    OUTPUT_DIR   = 'GaN_w/band_gap_bowing/AlScN_HSE/static_HSE_222'

    summary_df = compile_all_folders(GLOB_PATTERN)

    if summary_df is not None:
        out_path = os.path.join(OUTPUT_DIR, 'summary_coordination.csv')
        summary_df.to_csv(out_path, index=False)
        print(f"\nSaved: {out_path}  ({len(summary_df)} rows)")
        print(summary_df.to_string(index=False))
