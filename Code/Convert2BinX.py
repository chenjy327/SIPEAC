import pandas as pd
import os
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sample', type=str, default=64, help='merge name')

args = parser.parse_args()
sample = args.sample

#  nohup /data01/home/chenjy/.conda/envs/pacbio/bin/python /data02/chenjy/Step10.Convert2Bin20.245.4327B.py -s 4327B > /data02/chenjy/log/4327B.convertbin20.log 2>&1 &

barcode_file = f"/data02/chenjy/barcode_bin_info/{sample}.txt"
pacbio_info_file = f"/data02/chenjy/Dedup_fasta/{sample}_dedup.inST.info.csv"
to_dir = f'/data02/chenjy/Dedup_fasta/'
if not os.path.exists(to_dir):
    os.system(f'mkdir {to_dir}')

df = pd.read_csv(barcode_file, index_col=0, sep='\t', header=None)
df.columns = ['x', 'y']

# off_set = pd.read_csv(offset_file, index_col=0, sep='\t', header=None)
# # y minis offset
# OffsetY  = int(off_set.loc[sample][2])
# OffsetX = int(off_set.loc[sample][1])
OffsetY  = 0
OffsetX = 0
df1 = df[df['y'] > OffsetY]
df1['y'] = df1['y'] - OffsetY

df1 = df[df['x'] > OffsetX]
df1['x'] = df1['x'] - OffsetX
# add min and max value
bin_size = 20
x_min = np.min(df['x'])
x_max = np.max(df['x'])
y_min = np.min(df['y'])
y_max = np.max(df['y'])
df1.loc['pseudo_1'] = [x_min, y_min]
df1.loc['pseudo_2'] = [x_max, y_max]
print('get bin ID')

# binID
df1['x_transfer'] = ((df1['x'] - x_min) / bin_size + 1).astype('int32')
df1['y_transfer'] = ((df1['y'] - y_min) / bin_size + 1).astype('int32')
df1['new_Bin_ID'] = np.max(df1['x_transfer']) * (df1['y_transfer'] - 1) + df1['x_transfer']

#
print('get result')
df_info = pd.read_csv(pacbio_info_file, index_col=0, sep='\t')

used_bc = set(df1.index) & set(df_info['BCrev'])

df1 = df1.reindex(used_bc)

df1['new_Bin_ID'] = [f"BIN.{x}" for x in df1['new_Bin_ID']]
df1.to_csv(f'{to_dir}/{sample}_barcode_bin.info.20.txt', sep='\t')

dic = dict(zip(df1.index, df1['new_Bin_ID']))
df_info['binid'] = df_info['BCrev'].map(dic)

df_info = df_info.loc[~df_info['binid'].isna(),]

df_info.to_csv(f'{to_dir}/{sample}_dedup.info.BinID.20.csv', sep='\t')


