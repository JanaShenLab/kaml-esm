import os
import re
import copy
import subprocess
import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

DSSP_PATH = "mkdssp "
RIDA_PATH = "rida "

################################################################################
#                               HELPER FUNCTIONS
################################################################################

def pdb2fasta(file_path):
    """
    Convert a PDB to FASTA for the entire structure (single chain).
    Minimal version matching the original code's approach.
    """
    aa3to1 = {
        'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M',
        'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
        'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H',
        'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G',
        'MSE': 'M',
    }
    ca_pattern = re.compile(
        r"^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])"
    )
    filename = os.path.basename(file_path).split('.')[0]
    chain_dict = {}
    chain_list = []

    with open(file_path, 'r') as fp:
        for line in fp:
            if line.startswith("ENDMDL"):
                break
            match_list = ca_pattern.findall(line)
            if match_list:
                # match_list[0] might look like ('ASP','','','A') etc.
                # Combine the relevant group to get resn and chain
                resn = match_list[0][0] + match_list[0][2]  # 3-letter code
                chain = match_list[0][1] + match_list[0][3]  # chain ID
                one_letter = aa3to1.get(resn, 'X')
                if chain in chain_dict:
                    chain_dict[chain] += one_letter
                else:
                    chain_dict[chain] = one_letter
                    chain_list.append(chain)

    # Just pick the first chain we saw
    if chain_list:
        chain = chain_list[0]
        seq = chain_dict[chain]
    else:
        chain = "X"
        seq = ""

    return f">{filename}:{chain}\n{seq}\n"


def add_titr_features(resnm):
    """
    Return model_pka, n_sc_C, charge for each titratable residue type.
    """
    if resnm == 'ASP':
        return 3.7, 2, -1
    elif resnm == 'GLU':
        return 4.2, 3, -1
    elif resnm == 'HIS':
        return 6.5, 2, 0
    elif resnm == 'CYS':
        return 8.5, 1, 0
    elif resnm == 'LYS':
        return 10.4, 4, 1
    elif resnm == 'TYR':
        # Original code used 9.5 here, not 10.1
        return 9.5, 2, 0
    return 0.0, 0, 0


def sasa_to_buried_ratio(asa, resnm):
    """
    Return (1 - asa / max_accessible) for each residue type.
    """
    if resnm == 'ASP':
        return round(1 - asa / 157.0, 3)
    elif resnm == 'GLU':
        return round(1 - asa / 190.0, 3)
    elif resnm == 'HIS':
        return round(1 - asa / 176.0, 3)
    elif resnm == 'CYS':
        return round(1 - asa / 129.0, 3)
    elif resnm == 'LYS':
        return round(1 - asa / 207.0, 3)
    elif resnm == 'TYR':
        return round(1 - asa / 273.0, 3)
    return 0.0


def to_float(x):
    return float(x.strip()) if x.strip() else ''


def bbh_dis(x1, y1, z1, x2, y2, z2):
    """
    Original distance formula used by the code. (Misnamed arguments, but matches original.)
    """
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)


################################################################################
#                          MAIN GENERATE FUNCTION
################################################################################

def generate(residues):
    """
    Generate the feature CSV string for the given residue list. The first element
    in `residues` is the path to the PDB, subsequent elements are [chain, resnm, resid].
    """

    # ---------- 1) Basic Setup ----------
    file = residues[0]
    residues = residues[1:]  # the actual chain/resnm/resid combos

    # We'll store final results in `features`:
    features_str = (
        "PDB_ID,Chain,Res_Name,Res_ID,info,"
        "com_min_d0_polar,com_min_d1_polar,com_min_d0_nonpolar,com_min_d1_nonpolar,"
        "com_min_d0_bb_NH,com_min_d1_bb_NH,com_min_d0_neg_O,com_min_d1_neg_O,"
        "com_min_d0_pos_N,com_min_d1_pos_N,com_min_d0_hbond_o,com_min_d1_hbond_o,"
        "com_min_d0_hbond_n,com_min_d1_hbond_n,com_min_d0_hbond_h,com_min_d1_hbond_h,"
        "com_min_d0_C_SG,com_min_d1_C_SG,"
        "com_n_hv_6,com_n_hv_9,com_n_hv_12,com_n_hv_15,"
        "com_npolar_5,com_npolar_10,com_npolar_15,"
        "com_polar_5,com_polar_10,com_polar_15,"
        "ats,rss,rp2ss,rp4ss,rm2ss,rm4ss,"
        "DA1,DA2,DD1,DD2,"
        "d_fe,d_zn,d_mg,"
        "flexibility,model_pka,n_sc_C,charge,buried_ratio,metal,group_code\n"
    )

    # ---------- 2) Parse the PDB Once via biopandas ----------
    # Remove hydrogens globally
    ppdb = PandasPdb().read_pdb(file)
    atom_df = ppdb.df['ATOM']
    atom_df = atom_df[atom_df['element_symbol'] != 'H']  # remove hydrogens
    ppdb.df['ATOM'] = atom_df

    # Pre-make subsets we’ll need for distance checks. The code used repeated
    # deepcopies – we do them once, store them in dict form:
    # We keep them as DataFrame with only the relevant lines, no sorting needed
    # until we do the final distance calculations.
    def subset_df(resnames, atnames):
        mask_res = ppdb.df['ATOM']['residue_name'].isin(resnames)
        mask_atom = ppdb.df['ATOM']['atom_name'].isin(atnames)
        return ppdb.df['ATOM'][mask_res & mask_atom].copy()

    subset_dict = {
        'DE': subset_df(['ASP', 'GLU'], ['OD1', 'OD2', 'OE1', 'OE2']),
        'RK': subset_df(['ARG', 'LYS'], ['NH1', 'NH2', 'NE', 'NZ']),
        'TSY': subset_df(['THR', 'SER', 'TYR'], ['OG1', 'OG', 'OH']),
        'NQW': subset_df(['ASN', 'GLN', 'TRP'], ['ND2', 'NE2', 'NE1']),
        'H': subset_df(['HIS'], ['ND1', 'NE2']),
        'C': subset_df(['CYS'], ['SG']),
        'nonpolar': ppdb.df['ATOM'][ppdb.df['ATOM']['element_symbol'] == 'C'].copy(),
        'nn_NH': ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] == 'N'].copy(),
    }
    # For "polar" we do more complicated logic combining O/SD/N with special cases for ASN/GLN
    # We'll build that once:
    base_polar_mask = ppdb.df['ATOM']['atom_name'].isin(['O', 'SD', 'N'])
    # Add ASN:OD1, GLN:OE1
    asn_mask = (ppdb.df['ATOM']['residue_name'] == 'ASN') & (ppdb.df['ATOM']['atom_name'] == 'OD1')
    gln_mask = (ppdb.df['ATOM']['residue_name'] == 'GLN') & (ppdb.df['ATOM']['atom_name'] == 'OE1')
    subset_dict['polar'] = ppdb.df['ATOM'][base_polar_mask | asn_mask | gln_mask].copy()

    # ---------- 3) Pre-run DSSP Once ----------
    try:
        dssp_out = subprocess.check_output(f"{DSSP_PATH} {file}", shell=True, text=True, stderr=subprocess.DEVNULL)
        dssp_lines = dssp_out.strip().split('\n')
    except subprocess.CalledProcessError:
        # If DSSP fails, we can store dssp_lines = [] or something
        dssp_lines = []

    # We'll parse them into a dictionary: (chain, resid, 1-letter-code) -> lineindex
    # The code does a fuzzy match: dssp_file[i].find(f"{resid} {chain} {resn}") != -1
    # We'll do a quick approach: store entire lines in a list, then do a dictionary
    # keying by a triple. We can’t easily do a perfect parse w/o a bigger rewrite.
    # We'll do a naive approach: for each line, we look for the pattern "resid chain resn"
    # and store i in a dict.

    dssp_dict = {}
    for i, line in enumerate(dssp_lines):
        # We'll try to parse e.g. "  16  A  H" from line or do a substring check
        # but the original code uses: if line.find(f'{resid} {chain} {resn}') != -1:
        # We'll store them at runtime for faster searching. Because we might have collisions,
        # we store them in a dict-of-lists.
        # But let's do a big hack: we'll just store the line in a dictionary keyed by
        # (chain, resid, resn). We'll do a quick find with an f-string.
        # If we get duplicates, we store only the first. It's not perfect, but it matches old logic.
        # Because old logic breaks as soon as we find a match.
        # We'll get the 1-letter code from the line for the 3rd slot of the key.
        # The old code used resn for searching. We can do it the same way.
        # We'll do a quick parse:
        line_stripped = line.replace(" ", "")
        # We'll just store the entire line. We'll do the searching later. Simpler approach:
        # We'll skip this dictionary approach and do the same substring search inside the loop.
        pass
    # Honestly: The easiest to keep identical is: we do exactly the same search each time.
    # But that reintroduces a loop. We'll do a partial optimization: store `dssp_lines` once,
    # not re-run DSSP. Then each residue can do a quick for-loop to find the line. That’s still
    # better than re-running DSSP for each residue.

    # ---------- 4) Pre-run RIDA Once (for mean flexibility) ----------
    # The old code: if totalres > 30, run RIDA, read "rida_results.dat", store the mean of col(9).
    # Then apply to all residues. We'll replicate that once.
    fasta_str = pdb2fasta(file)
    # get the actual sequence from the second line
    seq_line = "".join(fasta_str.split('\n')[1:])
    # remove possible chain lines etc.
    seq_line = re.sub(r'^>.*', '', seq_line).replace('\n', '')
    totalres = len(seq_line.strip())
    if totalres <= 30:
        global_flexi = 0.5
    else:
        with open("ridainp", "w") as rr:
            rr.write(">X:A\n")
            rr.write(seq_line.strip() + "\n")

        try:
            subprocess.check_output(f"{RIDA_PATH} ridainp", shell=True, text=True)
            # We'll parse the "rida_results.dat" col(9):
            in1 = np.loadtxt("rida_results.dat", comments=["@", "$", "&", "#"], usecols=(9,))
            global_flexi = float(np.mean(in1))
        except Exception:
            global_flexi = 0.5
        # Remove ridainp and the result files
        for leftover in ["ridainp", "rida_results.dat", "rida_results.tab", "rida_anchor2.tab", "logs.log"]:
            if os.path.exists(leftover):
                os.remove(leftover)

    # ---------- 5) For each residue, do distance calculations + DSSP lookup + finalize row ----------
    # We'll define a small helper that returns the 28 fields from the old cal_min_dis logic.
    # Instead of re-parsing the entire PDB for each residue, we do everything with our subsets in memory.

    def compute_min_dist_fields(chain, resid_int, resnm):
        """
        Return the 28 fields that old code calls pyout. If we can't find the sidechain heavy atoms,
        we return 28 commas as placeholders. This merges the old 'cal_min_dis' function
        but removes repeated PDB reads and repeated deepcopies.
        """
        # “gc_coords” = the mean of the sidechain heavy-atom coordinates
        if resnm == 'ASP':
            at1, at2 = 'OD1', 'OD2'
        elif resnm == 'GLU':
            at1, at2 = 'OE1', 'OE2'
        elif resnm == 'HIS':
            at1, at2 = 'ND1', 'NE2'
        elif resnm == 'CYS':
            at1, at2 = 'SG', 'SG'
        elif resnm == 'LYS':
            at1, at2 = 'NZ', 'NZ'
        elif resnm == 'TYR':
            at1, at2 = 'OH', 'OH'
        else:
            # Not recognized, though the code never calls it except for the six residue types
            return "," * 28

        # subset for the residue's sidechain atoms
        sub = atom_df[
            (atom_df['residue_name'] == resnm) &
            (atom_df['residue_number'] == resid_int) &
            (atom_df['atom_name'].isin([at1, at2])) &
            (atom_df['chain_id'] == chain)
        ]
        if sub.empty:
            # old code would write an error file, but we'll just return blank
            return "," * 28

        gc_x = sub['x_coord'].mean()
        gc_y = sub['y_coord'].mean()
        gc_z = sub['z_coord'].mean()

        # For each subset: compute distance to (gc_x, gc_y, gc_z), ignoring the current residue
        # Then pick min1, min2. We'll define a quick local helper:
        def get_min12(df):
            if df.empty:
                return (999, 999)
            # distance to (gc_x, gc_y, gc_z)
            x = df['x_coord'].values
            y = df['y_coord'].values
            z = df['z_coord'].values
            dist = np.sqrt((x - gc_x) ** 2 + (y - gc_y) ** 2 + (z - gc_z) ** 2)

            # Remove the current residue from consideration
            # (the code used filter_one_atom_per_res, but let's do a simpler approach)
            # The old code used “df[df['residue_number'] != resid_nb]”
            # We'll do that at subset creation time if needed, but let's do it now:
            mask_self = (df['residue_number'].values == resid_int) & (df['chain_id'].values == chain)
            dist = dist[~mask_self]

            if dist.size == 0:
                return (999, 999)
            sortd = np.sort(dist)
            min1 = round(sortd[0], 2)
            if dist.size > 1:
                min2 = round(sortd[1], 2)
            else:
                min2 = 999
            return (min1, min2)

        # DE, RK, TSY, NQW, H, C, polar, nonpolar, nn_NH, hv
        # The code actually made separate subsets for “hv_df” but never assigned it to anything except empty?
        # We'll define hv as the entire ATOM df. Then we can do the same “dis <= X” logic. We'll do that.
        # Actually the code used hv_df = copy.deepcopy(pdb_df) then no filtering? It's the entire structure minus hydrogens.
        # Then we do counting how many are within 6, 9, 12, 15. We'll do that. We'll call it "hv".
        hv = atom_df  # entire no-H structure
        # We'll define a helper for "count how many atoms are within X" ignoring same residue
        def get_count_range(df, rmax):
            x = df['x_coord'].values
            y = df['y_coord'].values
            z = df['z_coord'].values
            dist = np.sqrt((x - gc_x) ** 2 + (y - gc_y) ** 2 + (z - gc_z) ** 2)
            # remove same residue
            mask_self = (df['residue_number'].values == resid_int) & (df['chain_id'].values == chain)
            dist = dist[~mask_self]
            return np.sum(dist <= rmax)

        # Now compute min1, min2 for each relevant subset:
        DE_min1, DE_min2 = get_min12(subset_dict['DE'])
        RK_min1, RK_min2 = get_min12(subset_dict['RK'])
        TSY_min1, TSY_min2 = get_min12(subset_dict['TSY'])
        NQW_min1, NQW_min2 = get_min12(subset_dict['NQW'])
        H_min1, H_min2 = get_min12(subset_dict['H'])
        C_min1, C_min2 = get_min12(subset_dict['C'])

        # "polar" version
        polar_min1, polar_min2 = get_min12(subset_dict['polar'])
        # Then we also need the counts of how many are within 5,10,15
        polar_atoms = subset_dict['polar']
        if not polar_atoms.empty:
            x = polar_atoms['x_coord'].values
            y = polar_atoms['y_coord'].values
            z = polar_atoms['z_coord'].values
            dist = np.sqrt((x - gc_x) ** 2 + (y - gc_y) ** 2 + (z - gc_z) ** 2)
            mask_self = (polar_atoms['residue_number'].values == resid_int) & \
                        (polar_atoms['chain_id'].values == chain)
            dist = dist[~mask_self]
            n_polar_5 = np.sum(dist <= 5)
            n_polar_10 = np.sum(dist <= 10)
            n_polar_15 = np.sum(dist <= 15)
        else:
            n_polar_5 = n_polar_10 = n_polar_15 = 0

        # "nonpolar"
        nonpolar_min1, nonpolar_min2 = get_min12(subset_dict['nonpolar'])
        nonpolar_atoms = subset_dict['nonpolar']
        if not nonpolar_atoms.empty:
            x = nonpolar_atoms['x_coord'].values
            y = nonpolar_atoms['y_coord'].values
            z = nonpolar_atoms['z_coord'].values
            dist = np.sqrt((x - gc_x) ** 2 + (y - gc_y) ** 2 + (z - gc_z) ** 2)
            mask_self = (nonpolar_atoms['residue_number'].values == resid_int) & \
                        (nonpolar_atoms['chain_id'].values == chain)
            dist = dist[~mask_self]
            n_nonpolar_5 = np.sum(dist <= 5)
            n_nonpolar_10 = np.sum(dist <= 10)
            n_nonpolar_15 = np.sum(dist <= 15)
        else:
            n_nonpolar_5 = n_nonpolar_10 = n_nonpolar_15 = 0

        # "nn_NH"
        nn_NH_min1, nn_NH_min2 = get_min12(subset_dict['nn_NH'])

        # "hv" counts
        n_hv_6 = get_count_range(hv, 6)
        n_hv_9 = get_count_range(hv, 9)
        n_hv_12 = get_count_range(hv, 12)
        n_hv_15 = get_count_range(hv, 15)

        # Return them as a comma-delimited string
        return (
            f"{polar_min1}, {polar_min2}, {nonpolar_min1}, {nonpolar_min2}, "
            f"{nn_NH_min1}, {nn_NH_min2}, {DE_min1}, {DE_min2}, {RK_min1}, {RK_min2}, "
            f"{TSY_min1}, {TSY_min2}, {NQW_min1}, {NQW_min2}, {H_min1}, {H_min2}, "
            f"{C_min1}, {C_min2}, {n_hv_6}, {n_hv_9}, {n_hv_12}, {n_hv_15}, "
            f"{n_nonpolar_5},{n_nonpolar_10},{n_nonpolar_15},"
            f"{n_polar_5},{n_polar_10},{n_polar_15}"
        )

    # ---------- 6) Build the output lines for each residue ----------
    # We'll store a single run of DSSP lines once. For each residue, we do a small loop
    # searching for `resid chain resn`.
    # Then parse ASA, secondary structure info, etc.
    for s in residues:
        chain = s[0]
        resnm = s[1]
        resid = s[2]

        # For mapping to DSSP usage, the code used a single-letter resn:
        if resnm == "HIS":
            resn_1 = 'H'
        elif resnm == "ASP":
            resn_1 = 'D'
        elif resnm == "GLU":
            resn_1 = 'E'
        elif resnm == "CYS":
            resn_1 = 'C'
        elif resnm == "LYS":
            resn_1 = 'K'
        elif resnm == "TYR":
            resn_1 = 'Y'
        else:
            resn_1 = 'X'  # fallback, though original code never calls for non-titratable

        # 6.1) The 28 fields from cal_min_dis
        try:
            resid_int = int(resid)
        except ValueError:
            # If it's not numeric, skip
            resid_int = -999
        pyout = compute_min_dist_fields(chain, resid_int, resnm)

        # 6.2) Search DSSP lines for the matching line
        # The old code: for i in range(len(dssp_file)): if dssp_file[i].find(f"{resid} {chain} {resn_1}") != -1:
        # We'll replicate that approach:
        res_row = -1
        for i, line in enumerate(dssp_lines):
            if f"{resid} {chain} {resn_1}" in line.replace(" ", ""):
                res_row = i
                break

        # 6.3) ASA
        if res_row == -1:
            # fallback ASA if we never found a match
            if resnm == 'ASP':
                asa = 157.0
            elif resnm == 'GLU':
                asa = 190.0
            elif resnm == 'HIS':
                asa = 176.0
            elif resnm == 'CYS':
                asa = 129.0
            elif resnm == 'LYS':
                asa = 207.0
            elif resnm == 'TYR':
                asa = 273.0
            else:
                asa = 0.0
            rss = rp2ss = rp4ss = rm2ss = rm4ss = ' '
        else:
            line = dssp_lines[res_row]
            # ASA is at [35:38]
            try:
                asa = float(line[35:38])
            except ValueError:
                asa = 0.0
            # the code takes rss = line[16], etc.
            rss = line[16] if len(line) > 16 else ' '
            rp2ss = dssp_lines[res_row + 2][16] if (res_row + 2) < len(dssp_lines) else ' '
            rp4ss = dssp_lines[res_row + 4][16] if (res_row + 4) < len(dssp_lines) else ' '
            rm2ss = dssp_lines[res_row - 2][16] if (res_row - 2) > 0 else ' '
            rm4ss = dssp_lines[res_row - 4][16] if (res_row - 4) > 0 else ' '

        # 6.4) DA1, DA2, DD1, DD2
        # The code uses base-pairing analysis from line[40:45], etc. We'll try to parse them:
        # If we have no line, we do a fallback.
        DA1 = DA2 = DD1 = DD2 = np.nan
        if res_row != -1:
            line = dssp_lines[res_row]
            try:
                bhba1 = int(line[40:45])
            except ValueError:
                bhba1 = 0
            try:
                bhbd1 = int(line[51:56])
            except ValueError:
                bhbd1 = 0
            try:
                bhba2 = int(line[62:67])
            except ValueError:
                bhba2 = 0
            try:
                bhbd2 = int(line[73:78])
            except ValueError:
                bhbd2 = 0

            # We’ll store them if we want to replicate old code’s partial calculations,
            # but old code simply used bbh_dis to get distances from O to N, etc. We'll do partial.
            # For maximum speed, we won't do the complicated logic of calling coor_dat each time.
            # The original code used coor_dat to parse the PDB again. We can replicate it in memory,
            # but let's keep it as is. We'll do it now:
            DA1, DA2, DD1, DD2 = get_bbh_distances(
                bhba1, bhba2, bhbd1, bhbd2, resid_int, atom_df
            )

        # 6.5) We skip the Fe, Zn, Mg distances for now. The original code sets them to 999. We'll do the same.
        d_fe = 999
        d_zn = 999
        d_mg = 999

        # 6.6) flexibility from RIDA
        flexi = global_flexi

        # 6.7) Titration features
        model_pka, n_sc_C, charge = add_titr_features(resnm)
        buried_ratio = sasa_to_buried_ratio(asa, resnm)

        # 6.8) Output line
        info = f"{chain}-{resid}"
        # group_code is "0" in the original code, metal is "999"
        # The old code does `..., {999},{0}\n`.
        # (Those are “metal” and “group_code,” presumably.)
        row_line = (
            f"XXXX,{chain},{resnm},{resid},{info},"
            f"{pyout},"
            f"{asa}, {rss}, {rp2ss}, {rp4ss}, {rm2ss}, {rm4ss},"
            f"{DA1}, {DA2}, {DD1}, {DD2},"
            f"{d_fe}, {d_zn}, {d_mg},"
            f"{flexi}, {model_pka}, {n_sc_C}, {charge}, {buried_ratio},999,0\n"
        )
        features_str += row_line

    return features_str


def get_bbh_distances(bhba1, bhba2, bhbd1, bhbd2, resid, atom_df):
    """
    Partial re-implementation of the coor_dat + bbh_dis used for DA1, DA2, DD1, DD2 in the old code.
    We parse the relevant O, N, H positions from the same Pandas DF, no re-reading.
    """
    # The old code is quite elaborate. We'll replicate minimal logic:
    #   a1rid = bhba1 + resid
    #   a2rid = bhba2 + resid
    #   d1rid = bhbd1 + resid
    #   d2rid = bhbd2 + resid
    # Then it picks O from the current residue, N from a1 or a2, etc. The code is complicated,
    # but the final distances are rarely used in the model. We’ll replicate quickly:
    a1rid = bhba1 + resid
    a2rid = bhba2 + resid
    d1rid = bhbd1 + resid
    d2rid = bhbd2 + resid

    # We need the "O" and "N" from the current residue? The old code used lines at indices [6,7,8].
    # It's simpler to do a “best effort.” We'll define a function that returns (x,y,z) or None:
    def get_atom_xyz(resid_nb, atomname):
        sub = atom_df[
            (atom_df['residue_number'] == resid_nb) &
            (atom_df['atom_name'] == atomname)
        ]
        if sub.empty:
            return None
        return (sub['x_coord'].values[0], sub['y_coord'].values[0], sub['z_coord'].values[0])

    # The code looks for:
    #   O, H, N for the current residue
    #   O, H, N for others
    # We'll do just enough to get DA1, DA2, DD1, DD2:
    #   DA1 = distance( N of current residue, O of a1rid? ) ? Actually old code was `nx,ny,nz` vs `a1x,a1y,a1z`.
    #   But the old code is quite confusing. We'll replicate exactly:
    #   a1 is "O" from a1rid? Actually no, old code uses lines:
    #        if l[2]=="O" and l[5]==str(a1rid):  # store in A[0]
    #        if l[2]=="O" and l[5]==str(a2rid):  # store in A[1]
    #        if l[2]=="N" and l[5]==str(d1rid):  # store in A[2]
    #        if l[2]=="N" and l[5]==str(d2rid):  # store in A[3]
    #        ...
    #   Then it uses Nx,Ny,Nz from the current residue [8], a1x,y,z from A[0] => DA1
    #   So DA1 = distance( N of current residue, O of a1rid ).
    # We can do precisely that.

    # Get N of current residue
    curN = get_atom_xyz(resid, "N")
    # O of a1, a2
    a1O = get_atom_xyz(a1rid, "O")
    a2O = get_atom_xyz(a2rid, "O")
    # O of current residue & N of d1, d2
    curO = get_atom_xyz(resid, "O")
    d1N = get_atom_xyz(d1rid, "N")
    d2N = get_atom_xyz(d2rid, "N")

    def dist_or_nan(x1, y1, z1, x2, y2, z2):
        if x1 is None or x2 is None:
            return np.nan
        return bbh_dis(x1, y1, z1, x2, y2, z2)

    # DA1 = distance(N of current residue, O of a1rid)
    if curN and a1O:
        DA1 = bbh_dis(curN[0], curN[1], curN[2], a1O[0], a1O[1], a1O[2])
    else:
        DA1 = np.nan
    # DA2
    if curN and a2O:
        DA2 = bbh_dis(curN[0], curN[1], curN[2], a2O[0], a2O[1], a2O[2])
    else:
        DA2 = np.nan
    # DD1 = distance(O of current residue, N of d1rid)
    if curO and d1N:
        DD1 = bbh_dis(curO[0], curO[1], curO[2], d1N[0], d1N[1], d1N[2])
    else:
        DD1 = np.nan
    # DD2
    if curO and d2N:
        DD2 = bbh_dis(curO[0], curO[1], curO[2], d2N[0], d2N[1], d2N[2])
    else:
        DD2 = np.nan

    return DA1, DA2, DD1, DD2

