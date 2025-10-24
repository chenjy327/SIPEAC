sample_file=~/sc_st_share_sample1.txt
bin_sizes=(110)
thresholds=(1 2 5 10 15 20)
cat ${sample_file} | while IFS= read -r sample_name; do
for bin_size in "${bin_sizes[@]}"; do
    for threshold in "${thresholds[@]}"; do
        sample_name=$(echo ${sample_name} | sed 's/\r//')
        sample_name=$(echo ${sample_name} | sed 's/\n//')
        echo "当前参数: sample = $sample_name bin_size=$bin_size, threshold=$threshold"

        pth="~/input_repair/bin${bin_size}/${sample_name}_mat.tsv"

        # 检查文件是否存在，如果不存在则跳过
        if [ ! -f "$pth" ]; then
            echo "文件 $pth 不存在，跳过"
            continue
        fi
        out_dir="~/repair_result/bin${bin_size}/${sample_name}_${threshold}/"
        rm -rf ${out_dir}
        mkdir -p ${out_dir}

         /Users/mac/Desktop/Chenjy/tools/conda/anaconda3/envs/repair/bin/repair analyze \
         -i ${pth} \
        -ap "^(IGH)" -bp "^(IGL|IGK)" -mxo ${threshold} -x 1 \
        -o ${out_dir}

        file_path=$(find "${out_dir}" -type f -name '*_analysis_result.tsv')
        /Users/mac/Desktop/Chenjy/tools/conda/anaconda3/envs/repair/bin/repair evaluate -i ${file_path} \
        -o ${out_dir} -gt "~/grouth_truth/bin${bin_size}/${sample_name}.tsv" -ci chainA chainB -ct \
        Hclone Lclone
    done
  done
done
