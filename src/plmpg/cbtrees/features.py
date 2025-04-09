import pandas as pd
import numpy as np
import sys
from biopandas.pdb import PandasPdb
import copy
import sys,os
import re
import subprocess

DSSP_PATH="mkdssp "
RIDA_PATH="rida "

# cal_min_dis_py functions 
def filter_one_atom_per_res(df, dis, resnb):
    df['dis'] = dis
    df = df[df['residue_number'] != resnb]
    df = df.sort_values(by = ['residue_number', 'dis'])
    df = df.drop_duplicates(subset = ['residue_number'], keep='first')
    return df

def filter_polar(df, resnb):
    df = df[df['residue_number'] != resnb]
    df = df.sort_values(by = ['residue_number', 'dis'])
    return df

def dis_polar(df, gc_coords):
    df['gc_x'] = gc_coords[0]
    df['gc_y'] = gc_coords[1]
    df['gc_z'] = gc_coords[2]
    df['dis'] = np.sqrt((df['x_coord'] - df['gc_x'])**2 + 
                 (df['y_coord'] - df['gc_y'])**2 + 
                 (df['z_coord'] - df['gc_z'])**2) 
    return df

def get_min1_min2(df):
    min1 = round(df['dis'].min(), 2)
    if len(df) > 1:
        min2 = round(df['dis'].nsmallest(2).iloc[-1], 2)
    else:
        min2 = 999
    return min1, min2

def cal_min_dis(file, chain, resid, resnm, txt_id): 
    nan=' ,'
    resid = int(resid) 
    with open(f"error_record_{txt_id}.txt", "w") as err_txt:
        try:
            pth = file
            gc_coords = []
            at1 = ''
            at2 = ''
            short_nm = ''
            if resnm == 'ASP':
                at1 = 'OD1'
                at2 = 'OD2'
                short_nm = 'D'
            elif resnm == 'GLU':
                at1 = 'OE1'
                at2 = 'OE2'
                short_nm = 'E'
            elif resnm == 'HIS':
                at1 = 'ND1'
                at2 = 'NE2'
                short_nm = 'H'    
            elif resnm == 'CYS':
                at1 = 'SG'
                at2 = 'SG'
                short_nm = 'C'
            elif resnm == 'LYS':
                at1 = 'NZ'
                at2 = 'NZ'
                short_nm = 'K'
            elif resnm == 'TYR':
                at1 = 'OH'
                at2 = 'OH'
                short_nm = 'Y'
            pdb_df = PandasPdb().read_pdb(pth)
            pdb_df.df['ATOM'] = pdb_df.df['ATOM'][pdb_df.df['ATOM']['element_symbol'] != 'H']
            at_df = pdb_df.df['ATOM'][(pdb_df.df['ATOM']['residue_name'] == resnm) & 
                                      (pdb_df.df['ATOM']['residue_number'] == resid) &
                                      (pdb_df.df['ATOM']['atom_name'].isin([at1, at2])) &
                                      (pdb_df.df['ATOM']['chain_id'] == chain)]
            if not at_df.empty:
                gc_coords = [at_df['x_coord'].mean(), at_df['y_coord'].mean(), at_df['z_coord'].mean()]
                DE_df = copy.deepcopy(pdb_df)
                RK_df = copy.deepcopy(pdb_df)
                TSY_df = copy.deepcopy(pdb_df)
                NQW_df = copy.deepcopy(pdb_df)
                H_df = copy.deepcopy(pdb_df)
                C_df = copy.deepcopy(pdb_df)
                nonpolar_df = copy.deepcopy(pdb_df)
                polar_df = copy.deepcopy(pdb_df)
                nn_NH_df = copy.deepcopy(pdb_df)
                hv_df = copy.deepcopy(pdb_df)
                DE_df.df['ATOM'] = DE_df.df['ATOM'] = DE_df.df['ATOM'][DE_df.df['ATOM']['residue_name'].isin(['ASP', 'GLU'])]
                DE_df.df['ATOM'] = DE_df.df['ATOM'][DE_df.df['ATOM']['atom_name'].isin(['OD1', 'OD2', 'OE1', 'OE2'])]
                RK_df.df['ATOM'] = RK_df.df['ATOM'] = RK_df.df['ATOM'][RK_df.df['ATOM']['residue_name'].isin(['ARG', 'LYS'])]
                RK_df.df['ATOM'] = RK_df.df['ATOM'][RK_df.df['ATOM']['atom_name'].isin(['NH1', 'NH2', 'NE', 'NZ'])]
                TSY_df.df['ATOM'] = TSY_df.df['ATOM'][TSY_df.df['ATOM']['residue_name'].isin(['THR', 'SER', 'TYR'])]
                TSY_df.df['ATOM'] = TSY_df.df['ATOM'][TSY_df.df['ATOM']['atom_name'].isin(['OG1', 'OG', 'OH'])]
                NQW_df.df['ATOM'] = NQW_df.df['ATOM'][NQW_df.df['ATOM']['residue_name'].isin(['ASN', 'GLN', 'TRP'])]
                NQW_df.df['ATOM'] = NQW_df.df['ATOM'][NQW_df.df['ATOM']['atom_name'].isin(['ND2', 'NE2', 'NE1'])]
                H_df.df['ATOM'] = H_df.df['ATOM'][H_df.df['ATOM']['residue_name'].isin(['HIS'])]
                H_df.df['ATOM'] = H_df.df['ATOM'][H_df.df['ATOM']['atom_name'].isin(['ND1', 'NE2'])]
                C_df.df['ATOM'] = C_df.df['ATOM'][C_df.df['ATOM']['residue_name'].isin(['CYS'])]
                C_df.df['ATOM'] = C_df.df['ATOM'][C_df.df['ATOM']['atom_name'].isin(['SG'])]
                nonpolar_df.df['ATOM'] = nonpolar_df.df['ATOM'][nonpolar_df.df['ATOM']['element_symbol'] == 'C']
                nn_NH_df.df['ATOM'] = nn_NH_df.df['ATOM'][nn_NH_df.df['ATOM']['atom_name'].isin(['N'])]
                polar_df = pdb_df.df['ATOM'][pdb_df.df['ATOM']['atom_name'].isin(['O', 'SD', 'N'])]
                df_tmp1 = pdb_df.df['ATOM'][(pdb_df.df['ATOM']['residue_name'] == 'ASN') &
                                            (pdb_df.df['ATOM']['atom_name'] == 'OD1')]
                df_tmp2 = pdb_df.df['ATOM'][(pdb_df.df['ATOM']['residue_name'] == 'GLN') &
                                            (pdb_df.df['ATOM']['atom_name'] == 'OE1')]
                polar_df = pd.concat([polar_df, df_tmp1, df_tmp2])
                if not DE_df.df['ATOM'].empty:
                    DE_dis = DE_df.distance(xyz = gc_coords, records = ('ATOM',))
                    DE_res_dis = filter_one_atom_per_res(DE_df.df['ATOM'], DE_dis, resid)
                    DE_min1, DE_min2 = get_min1_min2(DE_res_dis)
                else:
                    DE_min1 = DE_min2 = 999
                if not RK_df.df['ATOM'].empty:
                    RK_dis = RK_df.distance(xyz = gc_coords, records = ('ATOM',))
                    RK_res_dis = filter_one_atom_per_res(RK_df.df['ATOM'], RK_dis, resid)
                    RK_min1, RK_min2 = get_min1_min2(RK_res_dis)
                else:
                    RK_min1 = RK_min2 = 999
                if not TSY_df.df['ATOM'].empty:
                    TSY_dis = TSY_df.distance(xyz = gc_coords, records = ('ATOM',))
                    TSY_res_dis = filter_one_atom_per_res(TSY_df.df['ATOM'], TSY_dis, resid)
                    TSY_min1, TSY_min2 = get_min1_min2(TSY_res_dis)
                else:
                    TSY_min1 = TSY_min2 = 999
                if not NQW_df.df['ATOM'].empty:
                    NQW_dis = NQW_df.distance(xyz = gc_coords, records = ('ATOM',))
                    NQW_res_dis = filter_one_atom_per_res(NQW_df.df['ATOM'], NQW_dis, resid)
                    NQW_min1, NQW_min2 = get_min1_min2(NQW_res_dis)
                else:
                    NQW_min1 = NQW_min2 = 999
                if not H_df.df['ATOM'].empty:
                    H_dis = H_df.distance(xyz = gc_coords, records = ('ATOM',))
                    H_res_dis = filter_one_atom_per_res(H_df.df['ATOM'], H_dis, resid)
                    H_min1, H_min2 = get_min1_min2(H_res_dis)
                else:
                    H_min1 = H_min2 = 999
                if not C_df.df['ATOM'].empty:
                    C_dis = C_df.distance(xyz = gc_coords, records = ('ATOM',))
                    C_res_dis = filter_one_atom_per_res(C_df.df['ATOM'], C_dis, resid)
                    C_min1, C_min2 = get_min1_min2(C_res_dis)
                else:
                    C_min1 = C_min2 = 999
                if not polar_df.empty:
                    polar_df = dis_polar(polar_df, gc_coords)
                    polar_res_dis = filter_polar(polar_df, resid)
                    polar_min1, polar_min2 = get_min1_min2(polar_res_dis)
                    n_polar_5 = len(polar_res_dis[polar_res_dis['dis'] <= 5])
                    n_polar_10 = len(polar_res_dis[polar_res_dis['dis'] <= 10])
                    n_polar_15 = len(polar_res_dis[polar_res_dis['dis'] <= 15])

                if not nonpolar_df.df['ATOM'].empty:
                    nonpolar_dis = nonpolar_df.distance(xyz = gc_coords, records = ('ATOM',))

                    nonpolar_df = nonpolar_df.df['ATOM']
                    nonpolar_df['dis'] = nonpolar_dis
                    nonpolar_res_dis = filter_polar(nonpolar_df, resid)
                    nonpolar_min1, nonpolar_min2 = get_min1_min2(nonpolar_res_dis)
                    n_nonpolar_5 = len(nonpolar_res_dis[nonpolar_res_dis['dis'] <= 5])
                    n_nonpolar_10 = len(nonpolar_res_dis[nonpolar_res_dis['dis'] <= 10])
                    n_nonpolar_15 = len(nonpolar_res_dis[nonpolar_res_dis['dis'] <= 15])
                if not nn_NH_df.df['ATOM'].empty:
                    nn_NH_dis = nn_NH_df.distance(xyz = gc_coords, records = ('ATOM',))
                    nn_NH_res_dis = filter_one_atom_per_res(nn_NH_df.df['ATOM'], nn_NH_dis, resid)
                    nn_NH_min1, nn_NH_min2 = get_min1_min2(nn_NH_res_dis)
                if not hv_df.df['ATOM'].empty:
                    hv_dis = hv_df.distance(xyz = gc_coords, records = ('ATOM',))
                    hv_res_df = hv_df.df['ATOM']
                    hv_res_df['dis'] = hv_dis
                    hv_res_df = hv_res_df[hv_res_df['residue_number'] != resid]
                    n_hv_6 = len(hv_res_df[hv_res_df['dis'] <= 6])
                    n_hv_9 = len(hv_res_df[hv_res_df['dis'] <= 9])
                    n_hv_12 = len(hv_res_df[hv_res_df['dis'] <= 12])
                    n_hv_15 = len(hv_res_df[hv_res_df['dis'] <= 15])
                return f'{polar_min1}, {polar_min2}, {nonpolar_min1}, {nonpolar_min2}, {nn_NH_min1}, {nn_NH_min2}, {DE_min1}, {DE_min2}, {RK_min1}, {RK_min2}, {TSY_min1}, {TSY_min2}, {NQW_min1}, {NQW_min2}, {H_min1}, {H_min2}, {C_min1}, {C_min2}, {n_hv_6}, {n_hv_9}, {n_hv_12}, {n_hv_15}, {n_nonpolar_5},{n_nonpolar_10},{n_nonpolar_15},{n_polar_5},{n_polar_10},{n_polar_15}'
            else:
                err_txt.write(f'{chain}-{resid}-{resnm}, atom not in PDB file')
                return f'{nan}' * 28 
        except FileNotFoundError:
            err_txt.write(f'{chain}-{resid}-{resnm}, PDB file not exist')
            return f'{nan}' * 28
            pass
        except ValueError:
            err_txt.write(f'{chain}-{resid}-{resnm}, Value error')
            return f'{nan}' * 28 
            pass
    err_txt.close()

#pdb2fasta and flex scripts as functions 
def pdb2fasta(file):
    aa3to1={
   'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
   'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
   'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
   'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
   'MSE':'M',
    }

    ca_pattern=re.compile(r"^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
    filename=os.path.basename(file).split('.')[0]
    chain_dict=dict()
    chain_list=[]

    fp=open(file,'r')
    for line in fp.read().splitlines():
        if line.startswith("ENDMDL"):
            break
        match_list=ca_pattern.findall(line)
        if match_list:
            resn=match_list[0][0]+match_list[0][2]
            chain=match_list[0][1]+match_list[0][3]
            if chain in chain_dict:
                chain_dict[chain]+=aa3to1[resn]
            else:
                chain_dict[chain]=aa3to1[resn]
                chain_list.append(chain)
    fp.close()
    return '>%s:%s\n%s\n'%(filename,chain,chain_dict[chain])

def coor_dat(a1rid, a2rid, d1rid, d2rid, resid, file): 
    A = [''] * 9 
    in1=tuple(open(file,"r"))
    for i in range(len(in1)):
        if in1[i][0:4]=="ATOM":        
            l=in1[i].split()
            if l[2]=="O" and l[5]==str(a1rid):
                A[0] = in1[i]
            if l[2]=="O" and l[5]==str(a2rid):
                A[1] = in1[i]
            if l[2]=="N" and l[5]==str(d1rid):
                A[2] = in1[i]
            if l[2]=="N" and l[5]==str(d2rid):
                A[3] = in1[i]
            if l[2]=="H" and l[5]==str(d1rid):
                A[4] = in1[i] 
            if l[2]=="H" and l[5]==str(d2rid):
                A[5] = in1[i] 
            if l[2]=="O" and l[5]==str(resid):
                A[6] = in1[i] 
            if l[2]=="H" and l[5]==str(resid):
                A[7] = in1[i] 
            if l[2]=="N" and l[5]==str(resid):
                A[8] = in1[i] 
    return A

def bbh_dis(x1, y1, z1, x2, y2, z2):
    xdis = x1 - y1 
    ydis = y1 - y2  
    zdis = z1 - z2
    return np.sqrt(xdis**2 + ydis**2 + zdis**2) 

def to_float(x): 
    if x.strip() == '':
        return ''
    else: 
        return float(x) 

# return: model_pka, n_sc_C, charge 
def add_titr_features(resnm):
    """Add titratable residue related features    Signature: df --> df
    """
    if resnm == 'ASP':
        model_pka = float(3.7)
        n_sc_C = 2
        charge = -1
#             data_df.loc[index, 'charge_change'] = 'neu_to_neg'
    elif resnm == 'GLU':
        model_pka = float(4.2)
        n_sc_C = 3
        charge = -1
#             data_df.loc[index, 'charge_change'] = 'neu_to_neg'
    elif resnm == 'HIS':
        model_pka = float(6.5)
        n_sc_C = 2
        charge = 0
#             data_df.loc[index, 'charge_change'] = 'pos_to_neu'
    elif resnm == 'CYS':
        model_pka = float(8.5)
        n_sc_C = 1
        charge = 0
#             data_df.loc[index, 'charge_change'] = 'neu_to_neg'
    elif resnm == 'LYS':
        model_pka = float(10.4)
        n_sc_C = 4
        charge = 1
#             data_df.loc[index, 'charge_change'] = 'pos_to_neu'
    elif resnm == 'TYR':
        model_pka = float(9.5)
        n_sc_C = 2
        charge = 0
#             data_df.loc[index, 'charge_change'] = 'neu_to_neg'
    return model_pka, n_sc_C, charge 


def sasa_to_buried_ratio(ats,resnm):

    
        if resnm == 'ASP':
            buried_ratio  = round(1-ats/157.0, 3)
            
        elif resnm == 'GLU':
            buried_ratio  = round(1-ats/190.0, 3)
            
        elif resnm == 'HIS':
            buried_ratio  = round(1-ats/176.0, 3)
            
        elif resnm == 'CYS':
            buried_ratio  = round(1-ats/129.0, 3)

        elif resnm == 'LYS':
            buried_ratio = round(1-ats/207.0, 3)
            
        elif resnm == 'TYR':
            buried_ratio = round(1-ats/273.0, 3)
            
        return buried_ratio

# generate features
def generate(residues): 
    features = "PDB_ID,Chain,Res_Name,Res_ID,info,com_min_d0_polar,com_min_d1_polar,com_min_d0_nonpolar,com_min_d1_nonpolar,com_min_d0_bb_NH,com_min_d1_bb_NH,com_min_d0_neg_O,com_min_d1_neg_O,com_min_d0_pos_N,com_min_d1_pos_N,com_min_d0_hbond_o,com_min_d1_hbond_o,com_min_d0_hbond_n,com_min_d1_hbond_n,com_min_d0_hbond_h,com_min_d1_hbond_h,com_min_d0_C_SG,com_min_d1_C_SG,com_n_hv_6,com_n_hv_9,com_n_hv_12,com_n_hv_15,com_npolar_5,com_npolar_10,com_npolar_15,com_polar_5,com_polar_10,com_polar_15,ats,rss,rp2ss,rp4ss,rm2ss,rm4ss,DA1,DA2,DD1,DD2,d_fe,d_zn,d_mg,flexibility,model_pka,n_sc_C,charge,buried_ratio,metal,group_code\n"
    file = residues[0] 
    residues = residues[1:] 
    for s in residues: 
        chain = s[0] 
        resnm = s[1] 
        resid = s[2]  
        
        if resnm== "HIS": 
            resn='H'
        elif resnm== "ASP": 
            resn='D'
        elif resnm == "GLU":
            resn='E'
        elif resnm == "CYS":
            resn='C'
        elif resnm == "LYS":
            resn='K'
        elif resnm == "TYR":
            resn='Y'
            
        ## calculate distance to various groups ##
        pyout = cal_min_dis(file, chain, resid, resnm, '100') 

        # run DSSP and read results 
        dssp_out = subprocess.check_output(f'{DSSP_PATH} {file}', shell=True, text=True,stderr=subprocess.DEVNULL) 
        dssp_file= dssp_out.strip().split('\n')

        # run through dssp file for residue row
        res_row = 0
        for i in range(len(dssp_file)):
            if dssp_file[i].find(f'{resid} {chain} {resn}') != -1: 
                res_row = i 
                break 
        if res_row==0:
            if resnm == 'ASP':
                asa = 157.0
            
            elif resnm == 'GLU':
                asa  = 190.0
            
            elif resnm == 'HIS':
                asa  = 176.0
            
            elif resnm == 'CYS':
                asa = 129.0

            elif resnm == 'LYS':
                asa = 207.0
                
            elif resnm == 'TYR':
                asa = 273.0
        else:
            res_dssp = dssp_file[res_row]
            asa = float(res_dssp[35:38]) #accessible surface area 
        
        #protein secondary structure of res and residues that are 2 and 4 AAs away from res 
        rss = res_dssp[16] 
        rp2ss = dssp_file[res_row+2][16] if res_row+2 < len(dssp_file) else ' '
        rp4ss = dssp_file[res_row+4][16] if res_row+4 < len(dssp_file) else ' ' 
        rm2ss = dssp_file[res_row-2][16] 
        rm4ss = dssp_file[res_row-4][16] 
        
        bhba1 = int(res_dssp[40:45])
        bhbd1 = int(res_dssp[51:56])
        bhba2 = int(res_dssp[62:67])
        bhbd2 = int(res_dssp[73:78])
        
        a1rid = bhba1 + int(resid) 
        d1rid = bhbd1 + int(resid) 
        a2rid = bhba2 + int(resid) 
        d2rid = bhbd2 + int(resid) 
        #run coor_dat 
        coor_dis = coor_dat(a1rid, a2rid, d1rid, d2rid, resid, file) #bbha1, bbha2, bbhd1, bbhd2, bbhd1H, bbhd2H, O, H, N 
        ocor = coor_dis[6] 
        ox = to_float(ocor[30:38])
        oy = to_float(ocor[38:46])
        oz = to_float(ocor[46:54]) 
        ncor = coor_dis[8] 
        nx = to_float(ncor[30:38]) 
        ny = to_float(ncor[38:46]) 
        nz = to_float(ncor[46:54]) 
        hcor = coor_dis[7] 
        hx = to_float(hcor[30:38]) 
        hy = to_float(hcor[38:46]) 
        hz = to_float(hcor[46:54]) 
        
        # distance features from coor_dat and bbh_dis (DA1, DA2, DD1, DD2) 
        a1cor = coor_dis[0] 
        a1x = to_float(a1cor[30:38]) 
        a1y = to_float(a1cor[38:46]) 
        a1z = to_float(a1cor[46:54]) 
        if nx != '' and ny != '' and nz != '' and a1x != '' and a1y != '' and a1z != '': 
            DA1 = bbh_dis(nx, ny, nz, a1x, a1y, a1z) 
        
        a2cor = coor_dis[1] 
        a2x = to_float(a2cor[30:38])
        a2y = to_float(a2cor[38:46]) 
        a2z = to_float(a2cor[46:54])
        if nx != '' and ny != '' and nz != '' and a2x != '' and a2y != '' and a2z != '': 
            DA2 = bbh_dis(nx, ny, nz, a2x, a2y, a2z) 
        
        d1cor = coor_dis[2] 
        d1x = to_float(d1cor[30:38])
        d1y = to_float(d1cor[38:46]) 
        d1z = to_float(d1cor[46:54]) 
        if ox != '' and oy != '' and oz != '' and d1x != '' and d1y != '' and d1z != '':
            DD1 = bbh_dis(ox, oy, oz, d1x, d1y, d1z) 
        
        d2cor = coor_dis[3] 
        d2x = to_float(d2cor[30:38]) 
        d2y = to_float(d2cor[38:46]) 
        d2z = to_float(d2cor[46:54]) 
        if ox != '' and oy != '' and oz != '' and d2x != '' and d2y != '' and d2z != '':
            DD2 = bbh_dis(ox, oy, oz, d2x, d2y, d2z)
        
        com_min_d0_fe=999
        com_min_d0_zn=999
        com_min_d0_mg=999
            
        ## process residue numbering for RIDA 
        pdbfile = open(file, 'r')
        orig_resnm = [] 
        orig_resid = [] 
        for line in pdbfile: 
            x = line.split() 
            if x[0] == 'ATOM': 
                orig_resnm.append(x[3]) 
                if len(x[4]) > 4: 
                    orig_resid.append(x[4][:-4])
                else: 
                    orig_resid.append(x[5]) 
        pdbfile.close()
        i = 0
        # remove grouped sequences 
        while i < len(orig_resid) - 1: 
            if orig_resid[i] == orig_resid[i+1]:
                orig_resid.pop(i) 
                orig_resnm.pop(i) 
            else: 
                i += 1 
        res_index = orig_resid.index(resid) # index of resid for flexibility 
        
        #flexibility 
        fastastr = pdb2fasta(file)
        fasta = fastastr.split()[1:] 
        text = ''
        for line in fasta:
            text = text + line
        text = text.replace('\n','')
        rr=open("ridainp","w")
        rr.write(">X:A"+"\n")
        rr.write(text+"\n")
        
        rr.close()
        totalres = len(text) 
        text="ridainp"
        

        flexi = 0.5
        if totalres > 30: 
            rida_out = subprocess.check_output(f'{RIDA_PATH} {text}', shell = True, text = True) 
            in1=np.loadtxt("rida_results.dat",comments=["@","$","&","#"],usecols = (9))

            flexi = np.mean(in1)#rida[res_index][9] #MDP col. 
        
        #titr_features and buried_ratio 
        model_pka, n_sc_C, charge = add_titr_features(resnm) 
        buried_ratio = sasa_to_buried_ratio(asa,resnm) 
        
        #final outfile print 
        L = f'{s[0]}, {s[1]}, {s[2]}'
        info = f'{chain}-{resid}'
        features += f'XXXX,{L}, {info}, {pyout}, {asa}, {rss}, {rp2ss}, {rp4ss}, {rm2ss}, {rm4ss}, {DA1}, {DA2}, {DD1}, {DD2}, {com_min_d0_fe}, {com_min_d0_zn}, {com_min_d0_mg}, {flexi}, {model_pka}, {n_sc_C}, {charge}, {buried_ratio},{999},{0}\n'
    return features 
