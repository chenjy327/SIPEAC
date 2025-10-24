
import torch.nn.functional as F
import torch
import torch.nn as nn

def multi_positive_contrastive_loss(pairwise_score, temperature, mask_positive):
    logits = pairwise_score / temperature  # [B, B]

    exp_logits = torch.exp(logits)  # [B, B]
    denominator = exp_logits.sum(dim=1, keepdim=True)  # [B, 1]

    positive_logits = exp_logits * mask_positive  # [B, B]
    positive_sum = positive_logits.sum(dim=1)  # [B]

    positive_probs = positive_sum / (denominator.squeeze(1) + 1e-9)  # [B]

    positive_log_probs = torch.log(positive_probs + 1e-9)  # [B]
    loss = -positive_log_probs.mean() 

    return loss

class Repair_model(nn.Module): 
    def __init__(self,esm_pretrain_model_H,esm_pretrain_model_L,esm_alphabet, temperature,esm_last_layer):
        super(Repair_model, self).__init__()
        self.temperature = temperature
        self.esm_last_layer=esm_last_layer
        
        self.esm_alphabet =esm_alphabet
         
        self.esm_H = esm_pretrain_model_H
        self.esm_L = esm_pretrain_model_L
        
        for param in self.esm_L.parameters():
            param.requires_grad = False
        for layer in self.esm_L.layers[-3:]:
            for param in layer.parameters():
                param.requires_grad = True
                
        for param in self.esm_H.parameters():
            param.requires_grad = False
        for layer in self.esm_H.layers[-3:]:
            for param in layer.parameters():
                param.requires_grad = True
        
    
    def forward(self, IGH_tokens, IGL_tokens,IGH_mask,IGL_mask, pair_label):        
        # Encoder
        IGH_out = self.esm_H(IGH_tokens,repr_layers=[self.esm_last_layer]) 
        IGH_emb = IGH_out["representations"][self.esm_last_layer]
        IGL_out = self.esm_L(IGL_tokens,repr_layers=[self.esm_last_layer])
        IGL_emb = IGL_out["representations"][self.esm_last_layer]
        
        # mean pooling
        IGH_emb = self.masked_mean_pooling(IGH_emb,IGH_mask)
        IGL_emb = self.masked_mean_pooling(IGL_emb,IGL_mask)
        
        # cosine similarity
        pairwise_score = self.pairwise_cosine_similarity(IGH_emb,IGL_emb) # [B,B]

        # constractive loss
        loss_HL = multi_positive_contrastive_loss(pairwise_score, self.temperature, pair_label)
        loss_LH = multi_positive_contrastive_loss(pairwise_score.T, self.temperature, pair_label)
        loss = loss_HL+loss_LH
        
        score_diag = torch.diagonal(pairwise_score, 0)
        
        return score_diag, loss
    
    def get_emb(self, IGH_tokens, IGL_tokens,IGH_mask,IGL_mask):
        # Encoder
        IGH_out = self.esm_H(IGH_tokens,repr_layers=[self.esm_last_layer]) 
        IGH_emb = IGH_out["representations"][self.esm_last_layer]
        IGL_out = self.esm_L(IGL_tokens,repr_layers=[self.esm_last_layer])
        IGL_emb = IGL_out["representations"][self.esm_last_layer]
        
        # mean pooling
        IGH_emb = self.masked_mean_pooling(IGH_emb,IGH_mask) # [B,hid]
        IGL_emb = self.masked_mean_pooling(IGL_emb,IGL_mask) # [B,hid]
        return IGH_emb, IGL_emb
    
    def pairwise_cosine_similarity(self, emb1, emb2):
        emb1 = F.normalize(emb1, dim=1)
        emb2 = F.normalize(emb2, dim=1)
        return torch.einsum("id,jd->ij", emb1, emb2)
    
    def masked_mean_pooling(self,embeddings, mask):
        mask = mask.unsqueeze(-1)  # [batch_size, seq_len, 1]
        masked_embeddings = embeddings * mask  # [batch_size, seq_len, embed_dim]
        sum_embeddings = masked_embeddings.sum(dim=1)  # [batch_size, embed_dim]
        
        sum_mask = mask.sum(dim=1)  # [batch_size, 1]
        sum_mask = sum_mask.clamp(min=1e-9)  
        pooled_embeddings = sum_embeddings / sum_mask  # [batch_size, embed_dim]
    
        return pooled_embeddings