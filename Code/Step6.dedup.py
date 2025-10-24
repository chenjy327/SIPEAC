#!/usr/bin/env python
# coding=utf-8
import pysam
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('-r', '--rootdir', type=str, help='rootdir')
parser.add_argument('-s', '--sample', type=str, help='sample name')

args = parser.parse_args()

root_dir = args.rootdir
sample = args.sample

os.chdir(root_dir)

raw_bam = f"{root_dir}/Step5.PolyADetection/{sample}.flnct.bam"

raw_bf = pysam.AlignmentFile(raw_bam, check_sq=False)

dic = {}
dic_count = {}
for read in raw_bf:
    dic_tag = dict(read.tags)
    cur_ccs = read.qname
    cur_bc = dic_tag['XC']
    cur_umi = dic_tag['XM']
    cur_ident = cur_bc + '_' + cur_umi
    if read.seq:
        if cur_ident not in dic:
            dic_count[cur_ident] = 1
            dic[cur_ident] = read.seq
        else:
            if len(read.seq) > len(dic[cur_ident]):
                dic[cur_ident] = read.seq
                dic_count[cur_ident] += 1
raw_bf.close()

cur_m = 0
with open(f'{root_dir}/Step6.dedup/dedup.fasta', 'w') as fo:
    for k, v in dic.items():
        fo.write(f">molecule/{cur_m} full_length_coverage={dic_count[k]};length={len(v)};XC={k.split('_')[0]};XM={k.split('_')[1]}\n")
        fo.write(f"{v}\n")
        cur_m += 1

