import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv
from torch_geometric.nn import global_mean_pool
from torch_geometric.nn import global_max_pool
from typing import Union, Optional

# ptsV1layerV1nnV1p5p5, ptsV1layerV1nnV1p2p2, ptsV2layerV1nnV1p2p5 is not trained with this script but manually change the naming to keep consistent

class KaMLV2layerV1nnV1p2p2(torch.nn.Module):
    """
    KaMLV2: residue based node feature based on KaML-GAT but drop the feature for heteroatoms
    layerV1: 3 GAT + 3 FC
    nnV1: number of neurons: GAT 41, FC 41-32-16-1
    p2p2: dropout rate 0.2 for GAT & 0.2 for FC
    """
    def __init__(self):
        super().__init__()
        #torch.manual_seed(1)
        
        self.conv1 = GATConv(41, 41, dropout = 0.2) 
        self.conv2 = GATConv(41, 41, dropout = 0.2)
        self.conv3 = GATConv(41, 41, dropout = 0.2)
        
        self.fc1 = torch.nn.Linear(41, 32)
        self.bn1 = nn.BatchNorm1d(num_features=32)
        self.fc2 = torch.nn.Linear(32, 16)
        self.bn2 = nn.BatchNorm1d(num_features=16)
        self.fc3 = torch.nn.Linear(16, 1)

    def forward(self, x, edge_index, edge_attr, batch):
                
        #act = torch.nn.Tanhshrink()
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.2)
        
        h = self.conv1(x, edge_index, edge_attr)
        h = h.tanh()
        h = self.conv2(h, edge_index, edge_attr)
        h = h.tanh()
        h = self.conv3(h, edge_index, edge_attr)
        h = h.tanh()
        h = global_mean_pool(h, batch)
        
        h = act(self.bn1(self.fc1(h)))
        h = dropout(h)
        h = act(self.bn2(self.fc2(h)))
        h = dropout(h)
        h = self.fc3(h)
        
        out = h
        
        return out
    
class KaMLV2prottransV4layerV3nnV1p2p2(torch.nn.Module):
    """
    KaMLV2: residue based node feature based on KaML-GAT but drop the feature for heteroatoms
    prottransV4: painter_ProtTrans_v4_CBTree_cleanV1_split.ipynb
    layerV3: 3 GAT(KaML+ProtTrans) + 3 FC
    nnV1: number of neurons: KaML-GAT 41, ProtTransGAT 1058-512-256-32, FC 73-32-16-1
    p2p2: dropout rate 0.2 for GAT & 0.2 for FC
    """
    def __init__(self):
        super().__init__()
        #torch.manual_seed(1)
        ### GAT for KaML
        self.conv11 = GATConv(41, 41, dropout = 0.2)
        self.bn11 = nn.BatchNorm1d(41)
        self.conv12 = GATConv(41, 41, dropout = 0.2)
        self.bn12 = nn.BatchNorm1d(41)
        self.conv13 = GATConv(41, 41, dropout = 0.2)
        self.bn13 = nn.BatchNorm1d(41)
        ### GAT for ProtTrans
        self.conv21 = GATConv(1058, 512, dropout = 0.2)
        self.bn21 = nn.BatchNorm1d(512)
        self.conv22 = GATConv(512, 256, dropout = 0.2)
        self.bn22 = nn.BatchNorm1d(256)
        self.conv23 = GATConv(256, 32, dropout = 0.2)
        self.bn23 = nn.BatchNorm1d(32)
        
        self.fc1 = torch.nn.Linear(73, 32)
        self.bn1 = nn.BatchNorm1d(num_features=32)
        self.fc2 = torch.nn.Linear(32, 16)
        self.bn2 = nn.BatchNorm1d(num_features=16)
        self.fc3 = torch.nn.Linear(16, 1)

    def forward(self, x1, x2, edge_index1, edge_index2, edge_attr1, edge_attr2, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.2)
        
        ### GAT for KaML
        h = self.conv11(x1, edge_index1, edge_attr1)
        h = self.bn11(h)
        h = self.conv12(h, edge_index1, edge_attr1)
        h = self.bn12(h)
        h = self.conv13(h, edge_index1, edge_attr1)
        h = self.bn13(h)
        #print(h)
        #print(len(h))
        h = global_mean_pool(h, batch)
        #print('++++++++++++++++++++')
        ### GAT for ProtTrans
        H = self.conv21(x2, edge_index2, edge_attr2)
        H = self.bn21(H)
        H = self.conv22(H, edge_index2, edge_attr2)
        H = self.bn22(H)
        H = self.conv23(H, edge_index2, edge_attr2)
        H = self.bn23(H)
        #print(H)
        #print(len(H))
        H = global_max_pool(H, batch)
       
        ### FC layer
        ch = torch.cat((h, H), 1)
        #print(ch)
        ch = act(self.fc1(ch))
        ch = dropout(ch)
        ch = act(self.fc2(ch))
        ch = dropout(ch)
        ch = self.fc3(ch)
        
        out = ch
        
        return out
    
class KaMLV2prottransV5layerV4nnV1p2p2(torch.nn.Module):
    """
    KaMLV2: residue based node feature based on KaML-GAT but drop the feature for heteroatoms
    prottransV5: painter_ProtTrans_v5_CBTree_cleanV1_split.ipynb
    layerV4: 3 GAT(KaML) & 3 FC (ProtTrans) + 3 FC
    nnV1: number of neurons: KaML-GAT 41, ProtTransGAT 1024-512-256-32, FC 73-32-16-1
    p2p2: dropout rate 0.2 for GAT & 0.2 for FC
    """
    def __init__(self):
        super().__init__()
        #torch.manual_seed(1)
        ### GAT for KaML
        self.conv11 = GATConv(41, 41, dropout = 0.2)
        self.bn11 = nn.BatchNorm1d(41)
        self.conv12 = GATConv(41, 41, dropout = 0.2)
        self.bn12 = nn.BatchNorm1d(41)
        self.conv13 = GATConv(41, 41, dropout = 0.2)
        self.bn13 = nn.BatchNorm1d(41)
        ### FC for ProtTrans
        self.conv21 = torch.nn.Linear(1024, 512)
        self.bn21 = nn.BatchNorm1d(num_features=512)
        self.conv22 = torch.nn.Linear(512, 256)
        self.bn22 = nn.BatchNorm1d(num_features=256)
        self.conv23 = torch.nn.Linear(256, 32)
        self.bn23 = nn.BatchNorm1d(num_features=32)
        
        self.fc1 = torch.nn.Linear(73, 32)
        self.bn1 = nn.BatchNorm1d(num_features=32)
        self.fc2 = torch.nn.Linear(32, 16)
        self.bn2 = nn.BatchNorm1d(num_features=16)
        self.fc3 = torch.nn.Linear(16, 1)

    def forward(self, x1, x2, edge_index, edge_attr, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.2)
        
        ### GAT for KaML
        h = self.conv11(x1, edge_index, edge_attr)
        h = self.bn11(h)
        h = self.conv12(h, edge_index, edge_attr)
        h = self.bn12(h)
        h = self.conv13(h, edge_index, edge_attr)
        h = self.bn13(h)
        #print(h)
        #print(len(h))
        h = global_mean_pool(h, batch)
        #print('++++++++++++++++++++')
        ### FC for ProtTrans
        H = act(self.conv21(x2))
        H = dropout(H)
        H = act(self.conv22(H))
        H = dropout(H)
        H = act(self.conv23(H))
        #H = dropout(H)
        #print(H)
        #print(len(H))
       
        ### FC layer
        ch = torch.cat((h, H), 1)
        ch = act(self.fc1(ch))
        ch = dropout(ch)
        ch = act(self.fc2(ch))
        ch = dropout(ch)
        ch = self.fc3(ch)
        
        out = ch
        
        return out
    
class esm2V1layerV1p2(torch.nn.Module):
    """
    esm2V1: residue based node feature based on painter_ESM_V1.ipynb
    layerV1: 4 FC (1280-512-256-32-1)
    p2: dropout rate 0.2 for FC
    """
    def __init__(self):
        super().__init__()
        ### 3 FC layers
        self.fc1 = torch.nn.Linear(1280, 512)
        #self.bn1 = nn.BatchNorm1d(num_features=512)
        self.fc2 = torch.nn.Linear(512, 256)
        self.bn2 = nn.BatchNorm1d(num_features=256)
        self.fc3 = torch.nn.Linear(256, 32)
        self.bn3 = nn.BatchNorm1d(num_features=32)
        self.fc4 = torch.nn.Linear(32, 1)
        
    def forward(self, x, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.2)
        
        h = self.fc1(x)
        h = dropout(h)
        h = act(self.bn2(self.fc2(h)))
        h = dropout(h)
        h = act(self.bn3(self.fc3(h)))
        h = dropout(h)
        h = self.fc4(h)
        
        out = h
        
        return out
    
class esm2V1layerV3p2(torch.nn.Module):
    """
    esm2V1: residue based node feature based on painter_ESM_V1.ipynb
    layerV1: 1 FC (1280-1)
    p2: dropout rate 0.2 for FC
    """
    def __init__(self):
        super().__init__()
        ### 3 FC layers
        self.fc1 = torch.nn.Linear(1280, 1)
        
    def forward(self, x, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.2)
        
        h = self.fc1(x)
        
        out = h
        
        return out
    
class esm2V1layerV1p2_hidden(torch.nn.Module):
    """
    esm2V1: residue based node feature based on painter_ESM_V1.ipynb
    layerV1: 4 FC (1280-512-256-32-1)
    p2: dropout rate 0.2 for FC
    """
    def __init__(self):
        super().__init__()
        ### 3 FC layers
        self.fc1 = torch.nn.Linear(1280, 512)
        #self.bn1 = nn.BatchNorm1d(num_features=512)
        self.fc2 = torch.nn.Linear(512, 256)
        self.bn2 = nn.BatchNorm1d(num_features=256)
        self.fc3 = torch.nn.Linear(256, 32)
        self.bn3 = nn.BatchNorm1d(num_features=32)
        self.fc4 = torch.nn.Linear(32, 1)
        
    def forward(self, x, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0)
        
        h = self.fc1(x)
        h = dropout(h)
        h = act(self.bn2(self.fc2(h)))
        h = dropout(h)
        h = act(self.bn3(self.fc3(h)))
        hl = h
        h = dropout(h)
        h = self.fc4(h)
        
        out = h
        
        return out, hl
    
class esmCV1layerV1p2(torch.nn.Module):
    """
    esmCV1: residue based node feature based on esmC
    layerV1: 4 FC (2560-1024-512-64-1)
    p2: dropout rate 0.2 for FC
    """
    def __init__(self):
        super().__init__()
        ### 3 FC layers
        self.fc1 = torch.nn.Linear(2560, 1024)
        self.bn1 = nn.BatchNorm1d(num_features=1024)
        self.fc2 = torch.nn.Linear(1024, 512)
        self.bn2 = nn.BatchNorm1d(num_features=512)
        self.fc3 = torch.nn.Linear(512, 64)
        self.bn3 = nn.BatchNorm1d(num_features=64)
        self.fc4 = torch.nn.Linear(64, 1)
        
    def forward(self, x, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.2)
        
        h = act(self.bn1(self.fc1(x)))
        #h = self.fc1(x)
        h = dropout(h)
        h = act(self.bn2(self.fc2(h)))
        h = dropout(h)
        h = act(self.bn3(self.fc3(h)))
        h = dropout(h)
        h = self.fc4(h)
        
        out = h
        
        return out
    
class esmCV1layerV1p2_hidden(torch.nn.Module):
    """
    esmCV1: residue based node feature based on esmC
    layerV1: 4 FC (2560-1024-512-64-1)
    p2: dropout rate 0.2 for FC
    """
    def __init__(self):
        super().__init__()
        ### 3 FC layers
        self.fc1 = torch.nn.Linear(2560, 1024)
        self.bn1 = nn.BatchNorm1d(num_features=1024)
        self.fc2 = torch.nn.Linear(1024, 512)
        self.bn2 = nn.BatchNorm1d(num_features=512)
        self.fc3 = torch.nn.Linear(512, 64)
        self.bn3 = nn.BatchNorm1d(num_features=64)
        self.fc4 = torch.nn.Linear(64, 1)
        
    def forward(self, x, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0)
        
        h = act(self.bn1(self.fc1(x)))
        h = dropout(h)
        h = act(self.bn2(self.fc2(h)))
        h = dropout(h)
        h = act(self.bn3(self.fc3(h)))
        hl = h
        h = dropout(h)
        h = self.fc4(h)
        
        out = h
        
        return out, hl
     
class esm2V1pestoV0layerV1p2(torch.nn.Module):
    """
    esm2V1: residue based node feature based on painter_ESM_V1.ipynb
    pestoV0: residue embeddings based on pestoV0/cat_pesto.py (from E{n_expt}f0.pt)
    layerV1: 3 FC for esm2 (1280-512-256-32), 2 FC for pesto emb (64-32), 1 FC for fusion (64(32&32)-1) 
    p2: dropout rate 0.2 for all FC
    """
    def __init__(self):
        super().__init__()
        ### 3 FC layers
        self.fc1 = torch.nn.Linear(1280, 512)
        self.bn1 = nn.BatchNorm1d(num_features=512)
        self.fc2 = torch.nn.Linear(512, 256)
        self.bn2 = nn.BatchNorm1d(num_features=256)
        self.fc3 = torch.nn.Linear(256, 32)
        self.bn3 = nn.BatchNorm1d(num_features=32)
        
        self.fc4 = torch.nn.Linear(64, 32)
        self.bn4 = nn.BatchNorm1d(num_features=32)
        
        self.fc5 = torch.nn.Linear(64, 1)
        
    def forward(self, x, x2, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.2)
        ### FC for ems2V1
        h = act(self.bn1(self.fc1(x)))
        #h = self.fc1(x)
        h = dropout(h)
        h = act(self.bn2(self.fc2(h)))
        h = dropout(h)
        h = act(self.bn3(self.fc3(h)))
        h = dropout(h)
        
        ### FC for pestoV0
        H = act(self.bn4(self.fc4(x2)))
        H = dropout(H)
        
        ### FC for fusion
        ch = torch.cat((h, H), 1)
        ch = self.fc5(ch)
        
        out = ch
        
        return out
    
class esm2V1pestoV0layerV2p2(torch.nn.Module):
    """
    esm2V1: residue based node feature based on painter_ESM_V1.ipynb
    pestoV0: residue embeddings based on pestoV0/cat_pesto.py (by concat all the 10 fold models for one expt)
    layerV2: 1 FC for esm2 (1280-640), 0 FC for pesto emb (640), 4 FC for fusion (1280(640&640)-512-256-32-1) 
    p2: dropout rate 0.2 for all FC
    """
    def __init__(self):
        super().__init__()
        self.fc1 = torch.nn.Linear(1280, 640)
        self.bn1 = nn.BatchNorm1d(num_features=640)
        
        self.fc2 = torch.nn.Linear(1280, 512)
        self.bn2 = nn.BatchNorm1d(num_features=512)
        self.fc3 = torch.nn.Linear(512, 256)
        self.bn3 = nn.BatchNorm1d(num_features=256) 
        self.fc4 = torch.nn.Linear(256, 32)
        self.bn4 = nn.BatchNorm1d(num_features=32)
        self.fc5 = torch.nn.Linear(32, 1)
        
    def forward(self, x, x2, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.2)
        ### FC for ems2V1
        h = act(self.bn1(self.fc1(x)))
        h = dropout(h)
        
        ### FC for fusion
        ch = torch.cat((h, x2), 1)
        ch = act(self.bn2(self.fc2(ch)))
        ch = dropout(ch)
        ch = act(self.bn3(self.fc3(ch)))
        ch = dropout(ch)
        ch = act(self.bn4(self.fc4(ch)))
        ch = dropout(ch)
        ch = self.fc5(ch)
        
        out = ch
        
        return out
    
class esm2V1ctcV0layerV1p4(torch.nn.Module):
    """
    esm2V1: residue based node feature based on painter_ESM_V1_ctc.ipynb
    ctcV0: the contact map is extracted, residues w/ top5 contacts will be added to the MLP with contact values 
           multiplied to the original residue embedding
    layerV1: 4 FC (7680-3840-1024-256-1)
    p2: dropout rate 0.4 for FC
    """
    def __init__(self):
        super().__init__()
        ### 3 FC layers
        self.fc1 = torch.nn.Linear(7680, 3840)
        self.bn1 = nn.BatchNorm1d(num_features=3840)
        self.fc2 = torch.nn.Linear(3840, 1024)
        self.bn2 = nn.BatchNorm1d(num_features=1024)
        self.fc3 = torch.nn.Linear(1024, 256)
        self.bn3 = nn.BatchNorm1d(num_features=256)
        self.fc4 = torch.nn.Linear(256, 1)
        
    def forward(self, x, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.4)
        
        h = self.fc1(x)
        h = dropout(h)
        h = act(self.bn2(self.fc2(h)))
        h = dropout(h)
        h = act(self.bn3(self.fc3(h)))
        h = dropout(h)
        h = self.fc4(h)
        
        out = h
        
        return out


class esm2V2layerV1p2(torch.nn.Module):
    """
    esm2V2: residue based node feature based on painter_ESM_V2.ipynb
    layerV1: 4 FC (5120-2048-1024-256-1)
    p2: dropout rate 0.2 for FC
    """
    def __init__(self):
        super().__init__()
        ### 3 FC layers
        self.fc1 = torch.nn.Linear(5120, 2048)
        self.bn1 = nn.BatchNorm1d(num_features=2048)
        self.fc2 = torch.nn.Linear(2048, 1024)
        self.bn2 = nn.BatchNorm1d(num_features=1024)
        self.fc3 = torch.nn.Linear(1024, 256)
        self.bn3 = nn.BatchNorm1d(num_features=256)
        self.fc4 = torch.nn.Linear(256, 1)
        
    def forward(self, x, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.2)
        
        h = self.fc1(x)
        h = dropout(h)
        h = act(self.bn2(self.fc2(h)))
        h = dropout(h)
        h = act(self.bn3(self.fc3(h)))
        h = dropout(h)
        h = self.fc4(h)
        
        out = h
        
        return out
    

class KaMLV2esm2V1layerV1nnV1p2p2(torch.nn.Module):
    """
    KaMLV2: residue based node feature based on KaML-GAT but drop the feature for heteroatoms
    esm2V1: residue based node feature based on painter_ESM_V1.ipynb
    layerV1: 3 GAT(KaML) & 3 FC (esm2) + 3 FC
    nnV1: number of neurons: KaML-GAT 41, ProtTransGAT 1280-512-256-32, FC 73-32-16-1
    p2p2: dropout rate 0.2 for GAT & 0.2 for FC
    """
    def __init__(self):
        super().__init__()
        #torch.manual_seed(1)
        ### GAT for KaML
        self.conv11 = GATConv(41, 41, dropout = 0.2)
        self.bn11 = nn.BatchNorm1d(41)
        self.conv12 = GATConv(41, 41, dropout = 0.2)
        self.bn12 = nn.BatchNorm1d(41)
        self.conv13 = GATConv(41, 41, dropout = 0.2)
        self.bn13 = nn.BatchNorm1d(41)
        ### FC for ProtTrans
        self.conv21 = torch.nn.Linear(1280, 512)
        self.bn21 = nn.BatchNorm1d(num_features=512)
        self.conv22 = torch.nn.Linear(512, 256)
        self.bn22 = nn.BatchNorm1d(num_features=256)
        self.conv23 = torch.nn.Linear(256, 32)
        self.bn23 = nn.BatchNorm1d(num_features=32)
        
        self.fc1 = torch.nn.Linear(73, 32)
        self.bn1 = nn.BatchNorm1d(num_features=32)
        self.fc2 = torch.nn.Linear(32, 16)
        self.bn2 = nn.BatchNorm1d(num_features=16)
        self.fc3 = torch.nn.Linear(16, 1)

    def forward(self, x1, x2, edge_index, edge_attr, batch):
        act = F.relu
        dropout = torch.nn.Dropout(p = 0.2)
        
        ### GAT for KaML
        h = self.conv11(x1, edge_index, edge_attr)
        h = self.bn11(h)
        h = self.conv12(h, edge_index, edge_attr)
        h = self.bn12(h)
        h = self.conv13(h, edge_index, edge_attr)
        h = self.bn13(h)
        h = global_mean_pool(h, batch)

        ### FC for ProtTrans
        H = act(self.conv21(x2))
        H = dropout(H)
        H = act(self.conv22(H))
        H = dropout(H)
        H = act(self.conv23(H))
        H = dropout(H)
       
        ### FC layer
        ch = torch.cat((h, H), 1)
        ch = act(self.bn1(self.fc1(ch)))
        ch = dropout(ch)
        ch = act(self.bn2(self.fc2(ch)))
        ch = dropout(ch)
        ch = self.fc3(ch)
        
        out = ch
        
        return out
    
class esm2V1layerV1p2_DAGU(nn.Module):
    """
    esm2V1layerV1p2 + DAGU (Disjoined Adaptive Gating Unit)

    GOAL:
      - Take concatenated embeddings: [ESM | structure-derived features].
      - Normalize each block, compute per-sample/graph gates, rejoin.
      - Feed through the original 4-FC regression head.
      - Allow identity-like gate init and optional joint LayerNorm toggle.

    INPUT / OUTPUT:
      x:     [N, emb_dim]  (PyG/node-wise view; N = number of nodes/tokens total)
      batch: [N] long (optional). batch[i] == graph/sample id for node i.
             If None, all nodes are treated as a single sample.
      returns: [N, out_dim] (default out_dim=1 for regression)

    MAPPING from MLPHead8 (which used [B, L, D]):
      - MLPHead8 forward took x_BLD with shape [B, L, D].
      - This head expects x_ND with N = sum over B of L (i.e., flatten over L).
      - The 'batch' vector maps each of those N tokens back to its original sample b in [0..B-1].
      - MLPHead8's per-sample mean over L → here becomes a per-graph mean over nodes grouped by 'batch'.
      - Everything else (LN → gate nets → rejoin → joint LN) is behavior-equivalent.

    ARGS:
      emb_dim: total width of x (ESM + structure features).
      esm_dim: width of the leading ESM slice inside x (x[:, :esm_dim]).
      out_dim: regression output dim (1 by default).
      dropout_p: dropout prob inside the FC stack.
      gate_reduce: optional reduction ratio for gate nets (Squeeze-Excite style).
                   None means Linear(d→d). Example: 8 yields d→(d/8)→d.
      use_joint_norm: if True, apply LayerNorm(emb_dim) after gating; if False, skip it.
      freeze_epoch: (optional) call on_epoch_end(epoch) to freeze gates + piecewise LNs.

    SWAP-IN:
      - Same forward(x, batch) signature as your original head.
      - Same tower & dropout behavior; eval() disables dropout correctly.
    """
    def __init__(self,
                 gate_reduce: Optional[int] = None,
                 freeze_epoch: Optional[int] = None,
                 emb_dim: int = 1280 + 82,
                 esm_dim: int = 1280,
                 out_dim: int = 1,
                 dropout_p: float = 0.2,
                 use_joint_norm: bool = True
                 #gate_reduce: int | None = None,
                 
                 #freeze_epoch: int | None = None):
                 ):
        super().__init__()
        assert 0 <= esm_dim <= emb_dim, "esm_dim must be within [0, emb_dim]"
        self.emb_dim = emb_dim
        self.esm_dim = esm_dim
        self.feat_dim = emb_dim - esm_dim

        # --- Piecewise LayerNorms (pre-gating) ---
        self.norm_esm  = nn.LayerNorm(self.esm_dim)  if self.esm_dim  > 0 else None
        self.norm_feat = nn.LayerNorm(self.feat_dim) if self.feat_dim > 0 else None

        # --- Joint LayerNorm (post-rejoin), can be toggled off ---
        self.use_joint_norm = use_joint_norm
        self.norm_joint = nn.LayerNorm(emb_dim) if use_joint_norm else nn.Identity()

        # --- Gate nets: summarize normalized block → per-channel mask in (0,1) ---
        def make_gate(d: int):
            if d == 0:
                return None
            if gate_reduce and (d // gate_reduce) > 0:
                return nn.Sequential(
                    nn.Linear(d, d // gate_reduce),
                    nn.GELU(),
                    nn.Linear(d // gate_reduce, d),
                    nn.Sigmoid()
                )
            return nn.Sequential(nn.Linear(d, d), nn.Sigmoid())

        self.gate_net_esm  = make_gate(self.esm_dim)
        self.gate_net_feat = make_gate(self.feat_dim)

        # --- Original 4-layer regression tower (behavior preserved) ---
        self.fc1 = nn.Linear(emb_dim, 512)
        self.fc2 = nn.Linear(512, 256)
        self.bn2 = nn.BatchNorm1d(256)
        self.fc3 = nn.Linear(256, 32)
        self.bn3 = nn.BatchNorm1d(32)
        self.fc4 = nn.Linear(32, out_dim)
        self.dropout = nn.Dropout(p=dropout_p)

        # --- Optional freeze schedule (matches your MLPHead8 behavior) ---
        self.freeze_epoch = freeze_epoch
        self.gates_frozen = False

    # -------- Convenience: identity-ish gate initialization --------
    def init_gates_identity(self, bias_value: float = 4.0):
        """
        Initialize the *final* linear of each gate to zeros + positive bias so that
        sigmoid(bias) ~ 1.0 → near-pass-through at step 0.

        bias=4 → σ(4)=0.982; bias=5 → σ(5)=0.993. Use 0 if you want σ(0)=0.5.
        Call this right after constructing the module (before training).
        """
        for gate in (self.gate_net_esm, self.gate_net_feat):
            if gate is None:
                continue
            # find the last Linear layer in the Sequential
            last_linear = None
            for m in reversed(list(gate)):
                if isinstance(m, nn.Linear):
                    last_linear = m
                    break
            if last_linear is not None:
                nn.init.zeros_(last_linear.weight)
                nn.init.constant_(last_linear.bias, bias_value)

    # -------- Optional: freeze gating & piecewise norms after N epochs --------
    @torch.no_grad()
    def on_epoch_end(self, epoch: int):
        if self.freeze_epoch is None or self.gates_frozen:
            return
        if epoch >= self.freeze_epoch:
            for net in (self.gate_net_esm, self.gate_net_feat):
                if net is not None:
                    for p in net.parameters():
                        p.requires_grad = False
            for norm in (self.norm_esm, self.norm_feat):
                if norm is not None:
                    for p in norm.parameters():
                        p.requires_grad = False
            self.gates_frozen = True
            print(f"[DAGU] gates + piecewise LNs frozen at epoch {epoch}")

    # -------- Utility: mean over nodes per graph (PyG-style) --------
    #def _segment_mean(self, x: torch.Tensor, batch: torch.Tensor | None) -> torch.Tensor:
    def _segment_mean(self, x: torch.Tensor, batch: Union[torch.Tensor, None]) -> torch.Tensor:
        """
        x: [N, C], batch: [N] with values in [0..G-1].
        returns: [G, C]; if batch is None, returns [1, C].
        """
        if batch is None:
            return x.mean(dim=0, keepdim=True)
        batch = batch.long()
        G = int(batch.max().item()) + 1 if x.numel() > 0 else 0
        sums = x.new_zeros((G, x.size(1)))
        sums.index_add_(0, batch, x)
        counts = x.new_zeros((G,), dtype=x.dtype)
        counts.index_add_(0, batch, torch.ones_like(batch, dtype=x.dtype))
        return sums / counts.clamp_min(1).unsqueeze(1)

    # ------------------------------ Forward ------------------------------
    def forward(self, x: torch.Tensor, batch: Union[torch.Tensor, None]) -> torch.Tensor:
        """
        x: [N, emb_dim]; batch: [N] or None.
        """
        assert x.dim() == 2 and x.size(1) == self.emb_dim, \
            f"Expected x [N,{self.emb_dim}], got {tuple(x.shape)}"

        N = x.size(0)
        parts = []

        # ---- Split: [ESM | FEAT] ----
        if self.esm_dim > 0:
            esm = x[:, :self.esm_dim]
            esm_norm = self.norm_esm(esm) if self.norm_esm is not None else esm
            # summarize per graph/sample
            esm_summary = self._segment_mean(esm_norm, batch)  # [G, esm_dim] or [1, esm_dim]
            if self.gate_net_esm is not None:
                esm_mask = self.gate_net_esm(esm_summary)       # [G, esm_dim] / [1, esm_dim]
                esm_mask_nodes = esm_mask.expand(N, -1) if batch is None else esm_mask[batch]
                esm_g = esm_norm * esm_mask_nodes               # gate per-channel
            else:
                esm_g = esm_norm
            parts.append(esm_g)

        if self.feat_dim > 0:
            feat = x[:, self.esm_dim:]
            feat_norm = self.norm_feat(feat) if self.norm_feat is not None else feat
            feat_summary = self._segment_mean(feat_norm, batch) # [G, feat_dim] / [1, feat_dim]
            if self.gate_net_feat is not None:
                feat_mask = self.gate_net_feat(feat_summary)
                feat_mask_nodes = feat_mask.expand(N, -1) if batch is None else feat_mask[batch]
                feat_g = feat_norm * feat_mask_nodes
            else:
                feat_g = feat_norm
            parts.append(feat_g)

        # ---- Rejoin (+ optional joint LayerNorm) ----
        h = torch.cat(parts, dim=-1) if len(parts) > 1 else parts[0]  # [N, emb_dim]
        h = self.norm_joint(h)                                        # LN or Identity

        # ---- Original 4-FC regression tower ----
        h = self.fc1(h)
        h = self.dropout(h)
        h = F.relu(self.bn2(self.fc2(h)))
        h = self.dropout(h)
        h = F.relu(self.bn3(self.fc3(h)))
        h = self.dropout(h)
        out = self.fc4(h)  # [N, out_dim] (use MSE/MAE for regression)
        return out

def get_model(model_name):
    """
    Return module class based on model name.
    """
    model_dict = {
        'KaMLV2layerV1nnV1p2p2': KaMLV2layerV1nnV1p2p2(),
        'KaMLV2prottransV4layerV3nnV1p2p2': KaMLV2prottransV4layerV3nnV1p2p2(),
        'KaMLV2prottransV5layerV4nnV1p2p2': KaMLV2prottransV5layerV4nnV1p2p2(),
        'esm2V1layerV1p2': esm2V1layerV1p2(),
        'esm2V2layerV1p2': esm2V2layerV1p2(),
        'esm2V1ctcV0layerV1p4': esm2V1ctcV0layerV1p4(),
        'KaMLV2esm2V1layerV1nnV1p2p2': KaMLV2esm2V1layerV1nnV1p2p2(),
        'esm2V1pestoV0layerV1p2': esm2V1pestoV0layerV1p2(),
        'esm2V1pestoV0layerV2p2': esm2V1pestoV0layerV2p2(),
        'esmCV1layerV1p2': esmCV1layerV1p2(),
        'esm2V1layerV1p2_DAGU': esm2V1layerV1p2_DAGU(),
        'esm2V1layerV1p2_hidden': esm2V1layerV1p2_hidden(),
        'esmCV1layerV1p2_hidden': esmCV1layerV1p2_hidden(),
        'esm2V1layerV3p2': esm2V1layerV3p2(),
    }
    return model_dict[model_name]