import argparse
import numpy as np
import os, random, pickle
import datetime
from tqdm import tqdm

import torch
from torch.utils.data import RandomSampler,DataLoader
from data import *
from model import *
from esm import pretrained

from train import Seed_everything,Write_log,cal_pair_label,evaluate_all_candidate,esm_featurize

def train_and_predict(model_class, config, args):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    output_path = args.output_path

    esm_modelname = config['esm_modelname']
    esm_last_layer = config['esm_last_layer']

    seed = config['seed'] 
    temperature=config['temperature']
    num_workers = config['num_workers']
    batch_size = config['batch_size']
    Seed_everything(seed)

    test_data = pickle.load(open(args.test_data_path,'rb'))     

    test_dataset = ProteinGraphDataset(test_data, range(len(test_data)))
    test_dataloader = DataLoader(test_dataset, batch_size = batch_size, shuffle=False, drop_last=False, num_workers=num_workers, prefetch_factor=2)

    esm_pretrain_model_L, esm_alphabet = pretrained.load_model_and_alphabet(esm_modelname)
    esm_pretrain_model_H, esm_alphabet = pretrained.load_model_and_alphabet(esm_modelname)
    model = model_class(esm_pretrain_model_H,esm_pretrain_model_L,esm_alphabet,temperature,esm_last_layer).to(device)
    state_dict = torch.load(output_path + 'model.ckpt', device)
    model.load_state_dict(state_dict)
    model.eval()
    IGH_ids,IGL_ids,preds_matrix = evaluate_all_candidate(model,test_dataloader,model.esm_alphabet,device)

    test_pred_dict = {"IGH_ids":IGH_ids,'IGL_ids':IGL_ids,"preds_matrix":preds_matrix.detach().cpu().numpy()}
    with open(output_path + f"inference_pred_dict.pkl", "wb") as f:
        pickle.dump(test_pred_dict, f)


parser = argparse.ArgumentParser()
parser.add_argument("--test_data_path", type=str, default='./data_OAS/test_data.pkl',help='inference data path')
parser.add_argument("--output_path", type=str, default='./output/',help='Directory of the trained model. Must contain the model.ckpt file trained from main.py.')
args = parser.parse_args()

model_class = Repair_model
nn_config = {
    'esm_modelname': 'esm2_t33_650M_UR50D', 
    'esm_last_layer': 33,  
    'num_workers':8,
    'seed': 2025,
    'batch_size': 48,
    'temperature':0.05,
}
Seed_everything(nn_config['seed'])
if __name__ == '__main__':
    train_and_predict(model_class, nn_config, args)