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

dic = {}
f_constant = open(f'{root_dir}/Step4.PrimerDetection_2/constant.txt')
for line in f_constant.readlines():
    cur_ccs = line.split('\t')[0]
    cur_status = line.split('\t')[1]
    if cur_status != '-1':
        constant_start = int(line.split('\t')[2])
        constant_end = int(line.split('\t')[3])
        dic[cur_ccs] = [(constant_start, constant_end)]
f_constant.close()

f_pr1reverse = open(f'{root_dir}/Step4.PrimerDetection_2/Pr1.reverse.txt')
for line in f_pr1reverse.readlines():
    cur_ccs = line.split('\t')[0]
    cur_status = line.split('\t')[1]
    if cur_status != '-1':
        pr1_start = int(line.split('\t')[2])
        pr1_end = int(line.split('\t')[3])
        if cur_ccs in dic:
            dic[cur_ccs].append((pr1_start, pr1_end))
f_pr1reverse.close()


f_pr1 = open(f'{root_dir}/Step4.PrimerDetection_2/Pr1.txt')
for line in f_pr1.readlines():
    cur_ccs = line.split('\t')[0]
    cur_status = line.split('\t')[1]
    if cur_status != '-1':
        pr1_start = int(line.split('\t')[2])
        pr1_end = int(line.split('\t')[3])
        if cur_ccs in dic:
            dic[cur_ccs].append((pr1_start, pr1_end))
f_pr1.close()


# [(1115, 1130), (1171, 1191), (17, 37)]
# 第一组代表固定序列位置，第二组代表pr1 reverse的位置，第三组代表pr1的位置
# umi位置在固定序列前10个碱基
# barcode位置在固定序列后25个碱基
# transcript 序列在pr1到umi 

# 生成新的bam
dic_barcode = {}

# In [15]: dic['m64198e_220118_141015/67/ccs']
# Out[15]: [(1115, 1130), (1171, 1191), (1171, 1191)]


raw_bam = f'{root_dir}/Step3.ConstantSeq.Detection/{sample}.pA.bam'
new_bam = f'{root_dir}/Step4.PrimerDetection_2/{sample}.flnc.bam'

raw_bf = pysam.AlignmentFile(raw_bam, check_sq=False)
new_bf = pysam.AlignmentFile(new_bam, 'wb', template=raw_bf)

bc_info = open(root_dir + '/Step4.PrimerDetection_2/Bc_UMI_info.txt', 'w')
bc_info.write('ccs\tBC\tBC_revc\tUMI\tUMI_revc\n')
dic_base = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

for read in raw_bf:
    cur_ccs = read.qname
    if cur_ccs in dic:
        cur_info = dic[cur_ccs]
        seq = read.seq
        qual = read.qual
        if len(cur_info) >= 2:
            if 30 < cur_info[1][0] - cur_info[0][1] < 50:
                constant_seq = seq[cur_info[0][0]: cur_info[0][1]]
                pr1_reverse = seq[cur_info[1][0]: cur_info[1][1]]
                cur_umi_start = cur_info[0][0] - 10
                cur_umi_end = cur_info[0][0]
                cur_umi = seq[cur_umi_start: cur_umi_end]
                cur_barcode_start = cur_info[0][1]
                cur_barcode_end = cur_info[0][1] + 25
                cur_barcode = seq[cur_barcode_start: cur_barcode_end]
                if len(cur_info) == 2:
                    final_seq = seq[: cur_umi_start]
                    final_qual = qual[: cur_umi_start]
                elif len(cur_info) == 3:
                    final_seq = seq[cur_info[2][1]: cur_umi_start]
                    final_qual = qual[cur_info[2][1]: cur_umi_start]
                read.seq = final_seq
                read.qual = final_qual
                cur_tags = read.tags
                cur_tags.append(('XC', cur_barcode))
                cur_tags.append(('XM', cur_umi))
                cur_barcode_rev = ''.join(list(map(lambda x:dic_base[x], cur_barcode)))[::-1]
                cur_umi_rev = ''.join(list(map(lambda x:dic_base[x], cur_umi)))[::-1]
                bc_info.write(f'{cur_ccs}\t{cur_barcode}\t{cur_barcode_rev}\t{cur_umi}\t{cur_umi_rev}\n')
                read.tags = cur_tags
                new_bf.write(read)

raw_bf.close()
new_bf.close()


