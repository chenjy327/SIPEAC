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
raw_bam =  f'{root_dir}/Step1.CCS/{sample}.ccs.bam'
os.chdir(root_dir)


#sample="3373B_1"
#root_dir=f'/data02/chenjy/pacbio_result/{sample}'
#raw_bam = f'{root_dir}/Step1.CCS/{sample}.ccs.bam'
os.chdir(root_dir)

# parse cutadapt result 1
dic = {}
f = open(root_dir + '/Step3.ConstantSeq.Detection/info.txt')
for line in f.readlines():
    cur_ccs = line.split('\t')[0]
    cur_status = line.split('\t')[1]
    if cur_status != '-1':
        cur_seq = line.strip().split('\t')[-1]
        dic[cur_ccs] = cur_seq
    else:
        dic[cur_ccs] = 'other'
f.close()


# 调整为polyA方向
dic_base = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

new_bam = f'{root_dir}/Step3.ConstantSeq.Detection/{sample}.pA.bam'
raw_bf = pysam.AlignmentFile(raw_bam, check_sq=False)
new_bf = pysam.AlignmentFile(new_bam, 'wb', template=raw_bf)

for read in raw_bf:
    cur_qname = read.qname
    if dic[cur_qname] == 'pA':
        new_bf.write(read)
    elif dic[cur_qname] == 'pT':
        cur_seq = read.seq
        cur_quality = read.qual
        new_seq = ''.join(list(map(lambda x:dic_base[x], cur_seq)))[::-1]
        new_quality = cur_quality[::-1]
        read.seq = new_seq
        read.qual = new_quality
        new_bf.write(read)

raw_bf.close()
new_bf.close()

