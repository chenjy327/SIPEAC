#!/usr/bin/env python
# coding=utf-8
# 比较pacbio数据中找到的barcode和whitelist barcode

import os
import numpy as np
import pandas as pd
import argparse
import re
from Bio.Seq import Seq

merge_name = "4327B"
sample_list = "4327B"
barcode_file = f"/data02/chenjy/barcode_bin_info/{merge_name}.txt"
# nohup /data01/home/chenjy/.conda/envs/pacbio/bin/python /data02/chenjy/Step7.245.4327B.py > /data02/chenjy/log/4327B.step7.log 2>&1 &

sample_list = sample_list.split(',')


def parse_fa(f):
    seq_info = f.readline().strip()
    seq = f.readline().strip()
    return seq_info, seq


fa_dict = {}
reads_count = 0
for sample in sample_list:
    dedup_file = f"/data02/chenjy/Dedup_fasta/{sample}_dedup.fasta"

    rex = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XC=(\S+);XM=(\S+)')
    f = open(dedup_file, 'r')
    seq_info, seq = parse_fa(f)
    info_file = re.sub('dedup.fasta', 'dedup.info.csv', dedup_file)

    fo = open(info_file, 'w')
    fo.write("id\tUMI\tUMIrev\tBC\tBCrev\tlength\tcount\n")

    while seq_info:
        m = rex.match(seq_info)
        if m:
            fa_dict[seq_info] = seq
            _id, _count, _len, _bc, _umi = m.groups()
            fo.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(_id[1:], _umi, Seq(_umi).reverse_complement(), _bc,
                                                                  Seq(_bc).reverse_complement(), _len, _count))
        seq_info, seq = parse_fa(f)
        reads_count = reads_count + 1
    f.close()
    fo.close()

df_st = pd.read_csv(barcode_file, index_col=0, sep='\t', header=None)
st_barcode = set(df_st.index)


## merge fasta file
def check_mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


rex = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XC=(\S+);XM=(\S+)')
dic = {}
dic_count = {}

for mol, seq in fa_dict.items():
    m = rex.match(mol)
    _id, _count, _len, _bc, _umi = m.groups()
    cur_ident = _bc + '_' + _umi
    if cur_ident not in dic.keys():
        dic_count[cur_ident] = 1
        dic[cur_ident] = seq
    else:
        if len(seq) > len(dic[cur_ident]):
            dic[cur_ident] = seq
            dic_count[cur_ident] += 1

cur_m = 0
with open(f"/data02/chenjy/Dedup_fasta/{merge_name}_dedup.fasta", 'w') as fo:
    for k, v in dic.items():
        fo.write(
            f">molecule/{cur_m} full_length_coverage={dic_count[k]};length={len(v)};XC={k.split('_')[0]};XM={k.split('_')[1]}\n")
        fo.write(f"{v}\n")
        cur_m += 1

rex = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XC=(\S+);XM=(\S+)')

merge_fa_file = f"/data02/chenjy/Dedup_fasta/{merge_name}_dedup.fasta"
f = open(merge_fa_file, 'r')
seq_info, seq = parse_fa(f)
to_file = re.sub('dedup.fasta', 'dedup.inST.fasta', merge_fa_file)
fo = open(to_file, 'w')

c = 0
while seq_info:
    m = rex.match(seq_info)
    if m:
        _id, _count, _len, _bc, _umi = m.groups()
        if Seq(_bc).reverse_complement() in st_barcode:
            fo.write(seq_info)
            fo.write('\n')
            fo.write(seq)
            fo.write('\n')
            c = c + 1
    seq_info, seq = parse_fa(f)

f.close()
fo.close()


def parse_fa(f):
    seq_info = f.readline().strip()
    seq = f.readline().strip()
    return seq_info, seq


dedup_file = f"/data02/chenjy/Dedup_fasta/{merge_name}_dedup.inST.fasta"
rex = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XC=(\S+);XM=(\S+)')
f = open(dedup_file, 'r')
seq_info, seq = parse_fa(f)
info_file = re.sub('dedup.inST.fasta', 'dedup.inST.info.csv', dedup_file)
fo = open(info_file, 'w')
fo.write("id\tUMI\tUMIrev\tBC\tBCrev\tlength\tcount\n")
while seq_info:
    m = rex.match(seq_info)
    if m:
        _id, _count, _len, _bc, _umi = m.groups()
        fo.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(_id[1:], _umi, Seq(_umi).reverse_complement(), _bc,
                                                              Seq(_bc).reverse_complement(), _len, _count))
    seq_info, seq = parse_fa(f)

f.close()
fo.close()