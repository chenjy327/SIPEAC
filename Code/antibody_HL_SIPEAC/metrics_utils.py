import os
import torch

def topk_hit_ndcg_mrr_torch(preds, labels, ks=[1, 5, 10, 100]):
    """
    preds, labels: [N, M] torch tensors, float32 & int64
    """
    results = {}
    sorted_scores, sorted_indices = preds.sort(dim=1, descending=True)
    sorted_labels = torch.gather(labels, 1, sorted_indices)

    device = preds.device
    log_positions = torch.log2(torch.arange(2, max(ks)+2, device=device).float())

    for k in ks:
        topk_labels = sorted_labels[:, :k]
        hit_at_k = (topk_labels.sum(dim=1) > 0).float().mean().item()
        results[f'Hit@{k}'] = round(hit_at_k, 4)

        # NDCG@k
        gains = topk_labels.float() / log_positions[:k]
        dcg = gains.sum(dim=1)

        ideal_sorted, _ = labels.sort(dim=1, descending=True)
        ideal_topk = ideal_sorted[:, :k]
        ideal_gains = ideal_topk.float() / log_positions[:k]
        idcg = ideal_gains.sum(dim=1)

        ndcg = torch.where(idcg > 0, dcg / idcg, torch.zeros_like(dcg))
        results[f'NDCG@{k}'] = round(ndcg.mean().item(), 4)

    # Mean Rank & MRR
    found_mask = (sorted_labels == 1).any(dim=1)
    first_hits = (sorted_labels == 1).float().argmax(dim=1) + 1  # 1-based
    mean_rank = torch.where(found_mask, first_hits, torch.tensor(preds.shape[1] + 1, device=device))
    mrr = torch.where(found_mask, 1.0 / mean_rank.float(), torch.zeros_like(mean_rank, dtype=torch.float32))

    results["MeanRank"] = round(mean_rank.float().mean().item(), 4)
    results["MedianRank"] = round(mean_rank.float().median().item(), 4)
    results["MRR"] = round(mrr.mean().item(), 4)

    return results


def compute_metrics(IGH_ids, IGL_ids, preds_matrix, H2L_pair_dict, threshold=0.5, hit_ks=[1,5,10,100]):
    IGH_index = {hid: i for i, hid in enumerate(IGH_ids)}
    IGL_index = {lid: j for j, lid in enumerate(IGL_ids)}

    preds_matrix = preds_matrix.detach().cpu()
    
    device = preds_matrix.device
    N, M = preds_matrix.shape
    label_matrix = torch.zeros((N, M), dtype=torch.int64, device=device)

    for hid, lids in H2L_pair_dict.items():
        i = IGH_index.get(hid, None)
        if i is not None:
            for lid in lids:
                j = IGL_index.get(lid, None)
                if j is not None:
                    label_matrix[i, j] = 1
    heavy_topk_metrics = topk_hit_ndcg_mrr_torch(preds_matrix, label_matrix, ks=hit_ks)
    light_topk_metrics = topk_hit_ndcg_mrr_torch(preds_matrix.T, label_matrix.T, ks=hit_ks)

    return {
        **{f"heavy_{k}": v for k, v in heavy_topk_metrics.items()},
        **{f"light_{k}": v for k, v in light_topk_metrics.items()}
    }