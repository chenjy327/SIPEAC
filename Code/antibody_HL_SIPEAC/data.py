import torch
import torch.utils.data as data

class ProteinGraphDataset(data.Dataset):
    def __init__(self, dataset, index):
        super(ProteinGraphDataset, self).__init__()       
        index = list(set(index))
        self.dataset = [item for i,item in enumerate(dataset) if i in index] 
        self.init_pair_dict()        
        
    def __len__(self): return len(self.dataset)

    def get_csv(self): return self.dataset
    
    def init_pair_dict(self): 
        self.H2L_pair_dict = {hid:[] for hid in [data_dict['heavy_chain_id'] for data_dict in self.dataset]}
        self.L2H_pair_dict = {lid:[] for lid in [data_dict['light_chain_id'] for data_dict in self.dataset]}
        for data_dict in self.dataset:
            self.H2L_pair_dict[data_dict['heavy_chain_id']].append(data_dict['light_chain_id'] )
            self.L2H_pair_dict[data_dict['light_chain_id']].append(data_dict['heavy_chain_id'] )
    
    def __getitem__(self, idx): return self._featurize_graph(idx)
    
    def _featurize_graph(self, idx):
        antibody_data = self.dataset[idx]
        pdb_id = str(antibody_data['pairname']) # 唯一标识
        
        return {  
            'name': pdb_id, 
            'heavy_chain_seq':antibody_data['heavy_chain_seq'],
            'light_chain_seq':antibody_data['light_chain_seq'],
            'heavy_chain_id':antibody_data['heavy_chain_id'],
            'light_chain_id':antibody_data['light_chain_id'],
            }

