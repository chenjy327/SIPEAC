# sh ~/run_changeo_new.sh B4327 /Users/mac/Desktop/Rproject/Pacbio/Result/Dedup_fasta/B4327_dedup.inST.fasta /Users/mac/Desktop/Rproject/Pacbio/Result/Dedup_fasta /Users/mac/Desktop/Rproject/Pacbio/Result/Dedup_fasta/B4327
# ulimit -n 100000000000
ulimit -n 100000000000
sample_name=$1
input=$2
out_dir=$3
blastn_out=$4

echo $sample_name
python_path=/Users/mac/Desktop/Chenjy/tools/conda/anaconda3/envs/bracer/bin/python
rscript_path=/usr/local/bin/Rscript
igblast_path=/Users/mac/Desktop/Chenjy/tools/igblast/
NPROC=50



/Users/mac/Desktop/Chenjy/tools/conda/anaconda3/bin/AssignGenes.py igblast -s ${input} -b ${igblast_path} \
    --exec ${igblast_path}/bin/igblastn --cdb /Users/mac/Desktop/Chenjy/tools/igblast/database/imgt_human_ig_c \
   --organism human --loci ig --format blast --outdir ${out_dir} --outname ${sample_name} --nproc ${NPROC}

 /Users/mac/Desktop/Chenjy/tools/conda/anaconda3/bin/MakeDb.py igblast -i ${out_dir}/${sample_name}_igblast.fmt7  \
  -s  ${input} \
 -r /Users/mac/Desktop/Chenjy/tools/igblast/germlines/imgt/human/vdj/imgt_human_*.fasta  /Users/mac/Desktop/Chenjy/tools/igblast/germlines/imgt/human/constant/imgt_human_*.fasta \
 --extended


 /Users/mac/Desktop/Chenjy/tools/conda/anaconda3/bin/CreateGermlines.py -d ${out_dir}/${sample_name}_igblast_db-pass.tsv \
      -r /Users/mac/Desktop/Chenjy/tools/igblast/germlines/imgt/human/vdj/imgt_human_*.fasta

/Users/mac/Desktop/Chenjy/tools/conda/anaconda3/bin/ParseDb.py select -d ${out_dir}/${sample_name}_igblast_db-pass_germ-pass.tsv -f "productive" -u T TRUE


mkdir -p ${blastn_out}
cd ${blastn_out}
#${python_path} /Users/mac/Desktop/Chenjy/tools/bracer/bracer assemble -s Hsap -p ${NPROC} -c /Users/mac/Desktop/Chenjy/tools/bracer/bracer_config_noigblast \
#--assembled_file ${input} all_seq ${blastn_out}/${sample_name}


echo "-----------------------Extract BCR reads ----------------------------------\n"
#${python_path} ~/Step8.extract.BCR.reads.py \
#-n ${sample_name} \
#--outdir ${blastn_out} \
#--seq ${input} \
#--blast_rst ${blastn_out}/${sample_name}/all_seq/BLAST_output/blastsummary_H.txt

#
/usr/local/bin/Rscript /Users/mac/Desktop/Rproject/Pacbio/Pacbio-BCR-Pipeline/single_sample/add_info.R \
--sample ${sample_name} --outdir ${out_dir} --b_outdir ${blastn_out}

if test -e ${out_dir}/${sample_name}_threshold.txt; then
  echo "文件存在"
  threshold=$(cat ${out_dir}/${sample_name}_threshold.txt)
else
  echo "文件不存在"
    threshold=0.2
fi

#/Users/mac/Desktop/Chenjy/tools/conda/anaconda3/bin/DefineClones.py -d ${out_dir}/${sample_name}.tsv \
#--mode gene --act set --model ham --dist ${threshold} --sf junction --norm len \
#--nproc 50

#/Users/mac/Desktop/Chenjy/tools/conda/anaconda3/bin/DefineClones.py -d ${out_dir}/${sample_name}_family.tsv \
#--mode gene --act set --model ham --dist 0.3 --sf junction --norm len \
#--nproc 50

/Users/mac/Desktop/Chenjy/tools/conda/anaconda3/bin/DefineClones.py -d ${out_dir}/${sample_name}.tsv \
--mode gene --act set --model ham --dist 0.3 --sf junction --norm len \
--nproc 50


/usr/local/bin/Rscript /Users/mac/Desktop/Rproject/Pacbio/Pacbio-BCR-Pipeline/single_sample/rename_column.R --sample ${sample_name} --outdir ${out_dir}



#/usr/local/bin/Rscript /Users/mac/Desktop/Rproject/Pacbio/Pacbio-BCR-Pipeline/single_sample/rename_column_0.7.R --sample ${sample_name} --outdir ${out_dir}


/Users/mac/Desktop/Chenjy/tools/conda/anaconda3/bin/DefineClones.py -d ${out_dir}/${sample_name}_family.tsv \
--mode gene --act set --model ham --dist 0 --sf junction --norm len \
--nproc 50


# /Users/mac/Desktop/Chenjy/tools/conda/anaconda3/bin/python /Users/mac/Desktop/Rproject/Pacbio/Pacbio-BCR-Pipeline/single_sample/BCR_secret_membrane.py \
# -w /Users/mac/Desktop/Rproject/Pacbio/Pacbio_v4/changeo \
# -p merge \
# -f /Users/mac/Desktop/Rproject/Pacbio/Pacbio_v4/rawdata/merge.fa
#
#
#
# /Users/mac/Desktop/Chenjy/tools/conda/anaconda3/bin/DefineClones.py -d /Users/mac/Desktop/Rproject/Pacbio/Pacbio_v4/Fig1/2907TB共享/2907_igblast_db-pass_germ-pass_parse-select.tsv \
# --mode gene --act set --model ham --dist 0 --sf junction --norm len \
# --nproc 50