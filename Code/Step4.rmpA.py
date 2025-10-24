#!/usr/bin/env python
# coding=utf-8
import pysam
import argparse
import os
import array

parser = argparse.ArgumentParser()

parser.add_argument('-r', '--rootdir', type=str, help='rootdir')
parser.add_argument('-s', '--sample', type=str, help='sample name')

args = parser.parse_args()

root_dir = args.rootdir
# root_dir = '/data/yuchen/pacbio/stereo_seq/A4/'
sample = args.sample
# sample = 'A4'

os.chdir(root_dir)

dic = {}
f = open(f'{root_dir}/Step5.PolyADetection/pA.txt')
for line in f.readlines():
    cur_ccs = line.split('\t')[0]
    cur_status = line.split('\t')[1]
    if cur_status != '-1':
        pa_start = int(line.split('\t')[2])
        dic[cur_ccs] = pa_start
    else:
        dic[cur_ccs] = -1
f.close()


raw_bam = f'{root_dir}/Step4.PrimerDetection_2/{sample}.flnc.bam'
new_bam = f'{root_dir}/Step5.PolyADetection/{sample}.flnct.bam'

raw_bf = pysam.AlignmentFile(raw_bam, check_sq=False)
new_bf = pysam.AlignmentFile(new_bam, 'wb', template=raw_bf)

for read in raw_bf:
    cur_ccs = read.qname
    cur_seq = read.seq
    cur_qual = read.qual
    if dic[cur_ccs] == -1:
        new_bf.write(read)
    else:
        read.seq = cur_seq[:dic[cur_ccs]]
        read.qual = cur_qual[:dic[cur_ccs]]
        new_bf.write(read)
raw_bf.close()
new_bf.close()


