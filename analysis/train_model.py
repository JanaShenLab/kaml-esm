import torch
import os
import re
import os.path
from glob import glob
import multiprocessing as mp
import numpy as np
import pandas as pd
from torch.utils.data import Dataset, DataLoader, ConcatDataset, SubsetRandomSampler
from torch_geometric.data import DataLoader as GeometricDataLoader
from tqdm import tqdm
from torch_geometric.utils import dense_to_sparse, to_networkx
from torch_geometric.data import Data, Batch, DataListLoader
from torch_geometric.loader import DataLoader
#from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances
import networkx as nx
from torch.nn import Linear, L1Loss
from torch_geometric.nn import GATConv
from torch_geometric.nn import DataParallel as GeometricDataParallel
from torch_geometric.nn import global_mean_pool
from torch.nn import Linear, Parameter
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops, degree
import torch.nn.functional as F
import torch.nn as nn
#import h5py
from sklearn.metrics import r2_score
import warnings
from typing import Optional, Tuple, Union

import torch
import torch.nn.functional as F
from torch import Tensor
from torch.nn import Parameter

from torch_geometric.nn.conv import MessagePassing
from torch_geometric.nn.dense.linear import Linear
from torch_geometric.nn.inits import glorot, zeros, uniform
from torch_geometric.typing import NoneType  # noqa
from torch_geometric.typing import (
    Adj,
    OptPairTensor,
    OptTensor,
    Size,
    SparseTensor,
    torch_sparse,
)
from torch_geometric.utils import (
    add_self_loops,
    is_torch_sparse_tensor,
    remove_self_loops,
    softmax,
)
from torch_geometric.utils.sparse import set_sparse_value
import sys
sys.path.append('/home/mings/his_pka/pkPDB/')
from model_class import get_model

device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

model_pka_dict = {'D': 3.7, 'E': 4.2, 'H': 6.5, 'C': 8.5, 'K': 10.4, 'Y': 9.8}
def cal_dpka_std(ID, Expt_pKa):
    resnm = ID.split('_')[-3]
    dpka_std = ((Expt_pKa - model_pka_dict[resnm]) - 0.24) / 2.3   ###split_restype_v3
    #dpka_std = ((Expt_pKa - model_pka_dict[resnm]) - 0.04) / 1.46   ###split_cdhit_v1
    #dpka_std = ((Expt_pKa - model_pka_dict[resnm]) - 0.21) / 2.3   ####split_restype_EXP67S
    ###dpka_std = ((Expt_pKa - model_pka_dict[resnm]) - (-0.2)) / 2.12   ####split_cossim_v1
    return dpka_std

mc = 'esm2V1layerV1p2'
lr = 0.0005
LR = 'lr0005'
data_pth = "/home/mings/his_pka/pkPDB/cleanV1_PKAD3_propsample/20times_entrybased_replace"
model = get_model(mc).to(device)
batch_size = 64
layer = 31

PKAD3_del0 = ['P00829-0-S-1BMF_D_Y_368_418', 'P00829-0-S-1BMF_D_Y_345_395', 'P09147-0-S-2UDP_A_Y_149_149']

multi_conf_ID = ['P10599-2-S-1ERT_A_D_26_26', 'P10599-2-S-1TRW_A_D_26_26', 'P10599-1-S-1ERU_A_D_26_26', 
                 'P10599-1-S-1TRS_A_D_26_26', 'Q72F37-1-S-2L6C_A_H_33_33', 'Q72F37-2-S-2L6D_A_H_33_33', 
                 'P69905-0-S-1HHO_A_H_89_90', 'P69905-0-S-4HHB_A_H_89_90', 'P69905-0-S-3KMF_A_H_89_90', 
                 'P68871-0-S-1HHO_B_H_146_147', 'P68871-0-S-4HHB_B_H_146_147']

### load the pre-trained weight
model.load_state_dict(torch.load(f'../E0f0B.pth'))

optimizer = torch.optim.Adam(model.parameters(), lr=lr)
criterion = nn.MSELoss()

def train(r_tl, train_res_li, train_pred_li, train_expt_li, device):
    model.train()
    
    for data in train_loader:  # Iterate in batches over the training dataset.
        data = data.to(device)
        out = model(data.x, data.batch)  # Perform a single forward pass.
        train_loss = criterion(out, data.y)
        train_r2 = r2_score(data.y.detach().cpu().numpy(), out.detach().cpu().numpy())
        r_tl += train_loss.item()*len(data.resname)
        res = data.resname
        pka = data.y.cpu().numpy()
        OUT = out.detach().cpu().numpy()
        train_res_li.extend(res)
        train_pred_li.extend(OUT)
        train_expt_li.extend(pka)
        #train_r = torch.corrcoef(torch.tensor([data.y.detach().numpy(), out.detach().numpy()]))
        #loss = criterion(out.squeeze(), data.y.squeeze())
        train_loss.backward()  # Derive gradients.
        optimizer.step()  # Update parameters based on gradients.
        optimizer.zero_grad()  # Clear gradients.
    r_tl /= len(train_dataset_list)
    r_tl = r_tl ** 0.5
    return train_loss, train_r2, r_tl, train_res_li, train_pred_li, train_expt_li


def val(loader, r_vl, val_res_li, val_pred_li, val_expt_li, device):
    model.eval()

    for data in loader:  # Iterate in batches over the training/test dataset.
        data = data.to(device)
        out = model(data.x, data.batch)
        res = data.resname
        pka = data.y.cpu().numpy()
        OUT = out.detach().cpu().numpy()
        val_res_li.extend(res)
        val_pred_li.extend(OUT)
        val_expt_li.extend(pka)
        #mae = torch.mean(torch.abs(torch.sub(out, data.y)))
        loss = criterion(out, data.y)
        r2 = r2_score(data.y.detach().cpu().numpy(), out.detach().cpu().numpy())
        r_vl += loss.item()*len(data.resname)
        #r = torch.corrcoef([data.y.detach().numpy(), out.detach().numpy()])
        #r = torch.corrcoef(torch.tensor([data.y.detach().numpy(), out.detach().numpy()]))
    r_vl /= len(val_dataset_list)
    r_vl = r_vl ** 0.5
    return loss, r2, r_vl, val_res_li, val_pred_li, val_expt_li

def test(loader, test_res_li, test_pred_li, test_expt_li, device):
    model.eval()

    for data in loader:  # Iterate in batches over the training/test dataset.
        data = data.to(device)
        out = model(data.x, data.batch)
        res = data.resname
        pka = data.y.cpu().numpy()
        OUT = out.detach().cpu().numpy()
        test_res_li.extend(res)
        test_pred_li.extend(OUT)
        test_expt_li.extend(pka)
    return test_res_li, test_pred_li, test_expt_li


class EarlyStopper:
    def __init__(self, patience = 5, min_delta = 0.1):
        self.patience = patience
        self.min_delta = min_delta
        self.counter = 0
        self.min_val_loss = float('inf')

    def early_stop(self, train_loss, val_loss, train_res_li, train_pred_li, train_expt_li, val_res_li, val_pred_li, val_expt_li, test_res_li, test_pred_li, test_expt_li):
        if val_loss < self.min_val_loss:
            if abs(train_loss - val_loss) <= 0.5:
                #if train_r2 > 0.3:
                self.min_val_loss = val_loss
                torch.save(model.state_dict(), f'base/E{n_expt}f{fold}.pth')
                df_record.loc[len(df_record.index)] = [epoch, epoch, epoch, epoch, epoch]
                df_train_detail = pd.DataFrame(train_res_li)
                df_train_detail['expt'] = np.concatenate(train_expt_li)
                df_train_detail['pred'] = np.concatenate(train_pred_li)
                df_train_detail.to_csv(f'base/E{n_expt}f{fold}_traindtl.csv', index = False)
                df_val_detail = pd.DataFrame(val_res_li)
                df_val_detail['expt'] = np.concatenate(val_expt_li)
                df_val_detail['pred'] = np.concatenate(val_pred_li)
                df_val_detail.to_csv(f'base/E{n_expt}f{fold}_valdtl.csv', index = False)
                df_test_detail = pd.DataFrame(test_res_li)
                df_test_detail['expt'] = np.concatenate(test_expt_li)
                df_test_detail['pred'] = np.concatenate(test_pred_li)
                df_test_detail.to_csv(f'base/E{n_expt}f{fold}_testdtl.csv', index = False)
                self.counter = 0
        elif val_loss > (self.min_val_loss + self.min_delta):
            self.counter += 1
            if self.counter >= self.patience:
                return True
        return False

#n_aug_in_test = 0
for n_expt in [%n_expt%]:
    for fold in [%n_fold%]:
        df_val = pd.read_csv(f'/home/mings/his_pka/split_restype_v3/expt{n_expt}f{fold}_validation_10.csv')
        df_val = df_val[df_val['tag'] != 'CpHMD']   ### drop CpHMD auged CYS (save to split_restype_v3_finetune)
        df_val = df_val[~df_val['ID'].isin(PKAD3_del0)]
        df_val = df_val[~df_val['ID'].isin(multi_conf_ID)]
        df_val['dpka_stded'] = df_val.apply(lambda x: cal_dpka_std(x.ID, x.Expt_pKa), axis=1)
        #df_val = df_val[(df_val['Res_Name'] != 'HIS') & (df_val['Res_Name'] != 'LYS')]  ### acid
        df_val = df_val[(df_val['Res_Name'] == 'HIS') | (df_val['Res_Name'] == 'LYS')]  ### base
        df_train = pd.read_csv(f'/home/mings/his_pka/split_restype_v3/expt{n_expt}f{fold}_train_10.csv')
        df_train = df_train[df_train['tag'] != 'CpHMD']   ### drop CpHMD auged CYS (save to split_restype_v3_finetune)
        df_train = add_id(df_train)
        #df_train = df_train[~df_train['ID'].isin(PKAD3_del0)]
        df_train = df_train[~df_train['ID'].isin(multi_conf_ID)]
        df_train['dpka_stded'] = df_train.apply(lambda x: cal_dpka_std(x.ID, x.Expt_pKa), axis=1)
        #df_train = df_train[(df_train['Res_Name'] != 'HIS') & (df_train['Res_Name'] != 'LYS')]   ### acid
        df_train = df_train[(df_train['Res_Name'] == 'HIS') | (df_train['Res_Name'] == 'LYS')]   ### base
        df_test = pd.read_csv(f'/home/mings/his_pka/split_restype_v3/expt{n_expt}_test.csv')
        df_test = df_test[~df_test['ID'].isin(PKAD3_del0)]
        df_test['dpka_stded'] = df_test.apply(lambda x: cal_dpka_std(x.ID, x.Expt_pKa), axis=1)
        #df_test = df_test[(df_test['Res_Name'] != 'HIS') & (df_test['Res_Name'] != 'LYS')]   ### acid
        df_test = df_test[(df_test['Res_Name'] == 'HIS') | (df_test['Res_Name'] == 'LYS')]   ### base
        train_dataset_list = []
        val_dataset_list = []
        test_dataset_list = []
        for index, row in df_train.iterrows():
            ID = row['ID']
            dpka_stded = row['dpka_stded']
            if os.path.isfile(f'{data_pth}/esm2V1layerV1p2/lr0005/split_restype_v3_finetune/esm2_t33_650M_UR50D_feat_v1_split_restype_v3_layer31/{ID}.pt'):
                data = torch.load(f'{data_pth}/esm2V1layerV1p2/lr0005/split_restype_v3_finetune/esm2_t33_650M_UR50D_feat_v1_split_restype_v3_layer31/{ID}.pt')
                data.x = torch.tensor(data.x,dtype=torch.float32)
                data.y = torch.FloatTensor(np.array(dpka_stded)).view(-1, 1)
                train_dataset_list.append(data)
        for index, row in df_val.iterrows():
            ID = row['ID']
            dpka_stded = row['dpka_stded']
            if os.path.isfile(f'{data_pth}/esm2V1layerV1p2/lr0005/split_restype_v3_finetune/esm2_t33_650M_UR50D_feat_v1_split_restype_v3_layer31/{ID}.pt'):
                data = torch.load(f'{data_pth}/esm2V1layerV1p2/lr0005/split_restype_v3_finetune/esm2_t33_650M_UR50D_feat_v1_split_restype_v3_layer31/{ID}.pt')
                data.x = torch.tensor(data.x,dtype=torch.float32)
                data.y = torch.FloatTensor(np.array(dpka_stded)).view(-1, 1)
                val_dataset_list.append(data)
        for index, row in df_test.iterrows():
            ID = row['ID']
            dpka_stded = row['dpka_stded']
            if os.path.isfile(f'{data_pth}/esm2V1layerV1p2/lr0005/split_restype_v3_finetune/esm2_t33_650M_UR50D_feat_v1_split_restype_v3_layer31/{ID}.pt'):
                data = torch.load(f'{data_pth}/esm2V1layerV1p2/lr0005/split_restype_v3_finetune/esm2_t33_650M_UR50D_feat_v1_split_restype_v3_layer31/{ID}.pt')
                data.x = torch.tensor(data.x,dtype=torch.float32)
                data.y = torch.FloatTensor(np.array(dpka_stded)).view(-1, 1)
                test_dataset_list.append(data)
        print(len(train_dataset_list), len(val_dataset_list), len(test_dataset_list))
        train_loader = DataLoader(train_dataset_list, batch_size = batch_size, shuffle = True)
        val_loader = DataLoader(val_dataset_list, batch_size = batch_size, shuffle = True)
        test_loader = DataLoader(test_dataset_list, batch_size = 999, shuffle = False)
        df_record = pd.DataFrame(columns = ['epoch','train_RMSE_epoch','train_R2','val_RMSE_epoch','val_R2'])
        early_stopper = EarlyStopper()
        for epoch in range(1, 150):
            r_tl = 0.0
            r_vl = 0.0
            train_res_li = []
            train_pred_li = []
            train_expt_li = []
            val_res_li = []
            val_pred_li = []
            val_expt_li = []
            test_res_li = []
            test_pred_li = []
            test_expt_li = []
            train_loss, train_r2, r_tl, train_res_li, train_pred_li, train_expt_li = train(r_tl, train_res_li, train_pred_li, train_expt_li, device = device)
            val_loss, val_r2, r_vl, val_res_li, val_pred_li, val_expt_li = val(val_loader, r_vl, val_res_li, val_pred_li, val_expt_li, device = device)
            test_res_li, test_pred_li, test_expt_li = test(test_loader, test_res_li, test_pred_li, test_expt_li, device = device)
            df_record.loc[len(df_record.index)] = [epoch, r_tl, train_r2, r_vl, val_r2]
            print(f'Epoch: {epoch},train_loss:{round(float(r_tl),5)}, val_loss:{round(float(r_vl),5)}')
            if epoch > 0:
                if early_stopper.early_stop(r_tl, r_vl, train_res_li, train_pred_li, train_expt_li, val_res_li, val_pred_li, val_expt_li, test_res_li, test_pred_li, test_expt_li):
                    break
        df_record.to_csv(f'base/E{n_expt}f{fold}_training.csv', index = False)