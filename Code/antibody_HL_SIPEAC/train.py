import numpy as np
import os, random, pickle
import datetime
from tqdm import tqdm

import torch
from torch.utils.data import RandomSampler,DataLoader
from model import *
from data import *
from metrics_utils import *
import argparse
from esm import pretrained

def Seed_everything(seed=42):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True

def Write_log(logFile, text, isPrint=True):
    if isPrint:
        print(text)
    logFile.write(text)
    logFile.write('\n')


CLS_IDX=0
PADDING_idx=1
EOS_IDX=2
UNK_IDX=3
MASK_IDX=4
def esm_featurize(data,esm_alphabet,device):
    batch_converter = esm_alphabet.get_batch_converter()
    IGH_seqs = [(pair_idx,wt_seq) for pair_idx,wt_seq in zip(data['name'],data['heavy_chain_seq'])]
    IGL_seqs = [(pair_idx,wt_seq) for pair_idx,wt_seq in zip(data['name'],data['light_chain_seq'])]
    _,_,IGH_batch_tokens=batch_converter(IGH_seqs)
    IGH_mask = (IGH_batch_tokens != PADDING_idx) & (IGH_batch_tokens != EOS_IDX) & (IGH_batch_tokens != CLS_IDX)
    _,_,IGL_batch_tokens=batch_converter(IGL_seqs)
    IGL_mask = (IGL_batch_tokens != PADDING_idx) & (IGL_batch_tokens != EOS_IDX) & (IGL_batch_tokens != CLS_IDX)
    
    return IGH_batch_tokens.to(device),IGL_batch_tokens.to(device),IGH_mask.to(device),IGL_mask.to(device)

def cal_pair_label(data,H2L_pair_dict,device):
    heavy_chain_id = [hid for hid in data['heavy_chain_id']]
    light_chain_id = [lid for lid in data['light_chain_id']]
    
    B = len(heavy_chain_id)
    pair_label = torch.zeros((B,B), dtype=torch.bool, device=device)

    # Traverse the indexes of heavy_chain_id and light_chain_id to fill in the positive sample relationship
    for i, hid in enumerate(heavy_chain_id):
        valid_light_ids = H2L_pair_dict[hid]  
        for j, lid in enumerate(light_chain_id):
            if lid in valid_light_ids: 
                pair_label[i, j] = True 
    
    return pair_label

def evaluate_all_candidate(model,dataloder,esm_alphabet,device):       
    IGH_ids,IGL_ids=[],[]
    all_IGH_embedding,all_IGL_embedding=[],[]
    for data in tqdm(dataloder):
        with torch.no_grad():
            esm_featurize(data,esm_alphabet,device)
            IGH_embeddings,IGL_embeddings = model.get_emb(*esm_featurize(data,esm_alphabet,device))
            
            IGH_valid_index = []
            for index,IGH_id in enumerate(data['heavy_chain_id']):
                if IGH_id not in IGH_ids:
                    IGH_ids.append(IGH_id)
                    IGH_valid_index.append(index)
            all_IGH_embedding+= [IGH_embeddings[i] for i in IGH_valid_index]
            
            IGL_valid_index = []
            for index,IGL_id in enumerate(data['light_chain_id']):
                if IGL_id not in IGL_ids:
                    IGL_ids.append(IGL_id)
                    IGL_valid_index.append(index)
            all_IGL_embedding+= [IGL_embeddings[i] for i in IGL_valid_index]


    IGH_tensor = torch.stack(all_IGH_embedding)  # shape: (N_IGH, D)
    IGL_tensor = torch.stack(all_IGL_embedding)  # shape: (N_IGL, D)
    IGH_tensor = torch.nn.functional.normalize(IGH_tensor, dim=1)
    IGL_tensor = torch.nn.functional.normalize(IGL_tensor, dim=1)

    # Calculate the cosine similarity matrix: (N_IGH, D) x (D, N_IGL) = (N_IGH, N_IGL)
    similarity_matrix = IGH_tensor @ IGL_tensor.T  # shape: (N_IGH, N_IGL)            
    
    return IGH_ids,IGL_ids,similarity_matrix

def train_and_predict(model_class, config, args):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    output_path = args.output_path
    os.makedirs(output_path, exist_ok = True)
    
    esm_modelname = config['esm_modelname']
    esm_last_layer = config['esm_last_layer']

    lr = config['lr']
    epochs = config['epochs']
    patience = config['patience']
    batch_size = config['batch_size']
    seed = config['seed'] 
    warmup=config['warmup']
    temperature=config['temperature']
    num_workers = config['num_workers']
    Seed_everything(seed)
    

    # Training
    log = open(output_path + 'train.log','w', buffering=1)
    Write_log(log, str(config) + '\n')

    train_data = pickle.load(open(args.dataset_path + f"/train_data.pkl",'rb'))        

    split_ratio = 0.96
    index = list(range(len(train_data)))
    random.shuffle(index)
    split_idx = int(len(index) * split_ratio)
    train_index = index[:split_idx]
    valid_index = index[split_idx:]
      
    train_dataset = ProteinGraphDataset(train_data, train_index)
    sampler = RandomSampler(train_dataset, replacement=True, num_samples=int(len(train_dataset)))
    train_dataloader = DataLoader(train_dataset, batch_size = batch_size, sampler=sampler, shuffle=False, drop_last=True, num_workers=num_workers, prefetch_factor=2)

    valid_dataset = ProteinGraphDataset(train_data, valid_index)
    valid_dataloader = DataLoader(valid_dataset, batch_size = int(batch_size),shuffle=False, drop_last=False, num_workers=num_workers, prefetch_factor=2)
            
    esm_pretrain_model_L, esm_alphabet = pretrained.load_model_and_alphabet(esm_modelname)
    esm_pretrain_model_H, esm_alphabet = pretrained.load_model_and_alphabet(esm_modelname)
    model = model_class(esm_pretrain_model_H,esm_pretrain_model_L,esm_alphabet,temperature,esm_last_layer).to(device)
    optimizer = torch.optim.Adam(model.parameters(), betas=(0.9, 0.99), lr=lr, weight_decay=1e-5, eps=1e-5)
                
    best_valid_metric = -1e9
    not_improve_epochs = 0

    for epoch in range(epochs):
        train_loss = 0
        train_num = 0
        model.train()

        train_preds = []
        bar = tqdm(train_dataloader)
        for data in bar:
            optimizer.zero_grad()       
            IGH_tokens, IGL_tokens,IGH_mask,IGL_mask = esm_featurize(data,esm_alphabet,device)
            pair_label = cal_pair_label(data,train_dataset.H2L_pair_dict,device)
            outputs,loss = model(IGH_tokens, IGL_tokens,IGH_mask,IGL_mask,pair_label)
            loss.backward()
            optimizer.step()
            
            train_preds.append(outputs.detach().cpu().numpy())
            
            train_num += len(outputs)
            train_loss += len(outputs) * loss.item()
            bar.set_description('loss: %.4f' % (loss.item()))
            
        train_loss /= train_num
        
        torch.cuda.empty_cache()

        if epoch<warmup:
            continue
        
        # Evaluate
        model.eval()
        IGH_ids,IGL_ids,preds_matrix = evaluate_all_candidate(model,valid_dataloader,model.esm_alphabet,device)
        valid_metrics = compute_metrics(IGH_ids,IGL_ids,preds_matrix,valid_dataset.H2L_pair_dict)
                    
        write_to_log = f'[epoch {epoch}] lr: {lr:.4f} train_loss: {train_loss:.4f} | [valid] '+ ", ".join([f"{key}: {value}" for key,value in valid_metrics.items()])
        print(valid_metrics)
        metric_for_stop = np.mean([valid_metrics['heavy_Hit@1'],valid_metrics['light_Hit@1'],valid_metrics['heavy_MRR'],valid_metrics['light_MRR']])
        if metric_for_stop >  best_valid_metric: 
            torch.save(model.state_dict(), output_path + 'model.ckpt')
            not_improve_epochs = 0
            best_valid_metric = metric_for_stop
            
            Write_log(log,write_to_log)
        else:
            not_improve_epochs += 1
            Write_log(log,write_to_log+', NIE +1 ---> %s'%not_improve_epochs)

            if not_improve_epochs >= patience:
                break
            
    # Use the best epoch to test validation again and save some prediction results
    model = model_class(esm_pretrain_model_H,esm_pretrain_model_L,esm_alphabet,temperature,esm_last_layer).to(device)
    state_dict = torch.load(output_path + 'model.ckpt', device)
    model.load_state_dict(state_dict)
    model.eval()
    IGH_ids,IGL_ids,preds_matrix = evaluate_all_candidate(model,valid_dataloader,model.esm_alphabet,device)
    valid_metrics = compute_metrics(IGH_ids,IGL_ids,preds_matrix,valid_dataset.H2L_pair_dict)

    Write_log(log,f'best epoch | [valid] ' + ", ".join([f"{key}: {value}" for key,value in valid_metrics.items()]))


    torch.cuda.empty_cache()

    # Test
    test_data = pickle.load(open(args.dataset_path + f"/test_data.pkl",'rb'))     

    test_dataset = ProteinGraphDataset(test_data, range(len(test_data)))
    test_dataloader = DataLoader(test_dataset, batch_size = batch_size, shuffle=False, drop_last=False, num_workers=num_workers, prefetch_factor=2)

    model = model_class(esm_pretrain_model_H,esm_pretrain_model_L,esm_alphabet,temperature,esm_last_layer).to(device)
    state_dict = torch.load(output_path + 'model.ckpt', device)
    model.load_state_dict(state_dict)
    model.eval()
    IGH_ids,IGL_ids,preds_matrix = evaluate_all_candidate(model,test_dataloader,model.esm_alphabet,device)
    test_metrics = compute_metrics(IGH_ids,IGL_ids,preds_matrix,test_dataset.H2L_pair_dict)
    
    Write_log(log,f'[test] '+", ".join([f"{key}: {value}" for key,value in test_metrics.items()]))

    test_pred_dict = {"IGH_ids":IGH_ids,'IGL_ids':IGL_ids,"preds_matrix":preds_matrix.detach().cpu().numpy()}
    with open(output_path + f"test_pred_dict.pkl", "wb") as f:
        pickle.dump(test_pred_dict, f)
    log.close()        




parser = argparse.ArgumentParser()
parser.add_argument("--dataset_path", type=str, default='./data_OAS/')
parser.add_argument("--output_path", type=str, default='./output/')
args = parser.parse_args()

model_class = Repair_model
nn_config = {
    'esm_modelname': 'esm2_t33_650M_UR50D', 
    'esm_last_layer': 33,  
    'num_workers':8,
    'lr': 1e-4,  
    'seed': 2025,
    'epochs': 1,
    'warmup':0,
    'patience': 10,
    'batch_size': 48,
    'temperature':0.05,
}
Seed_everything(nn_config['seed'])
if __name__ == '__main__':
    train_and_predict(model_class, nn_config, args)