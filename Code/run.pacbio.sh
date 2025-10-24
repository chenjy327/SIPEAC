# sample name
sample="4309B"
rootdir=/data02/chenjy/pacbio_result
ccs_bam=/data02/chenjy/LBFC20241075-05/24TW02532/r84236_241207_001_1_C01/K8.bc2048.bc2048.HiFi.bam
# nohup sh /data02/chenjy/run.pacbio.245.4309B.sh > /data02/chenjy/log/4309B.log 2>&1 &

echo "--------------------------processing sample $1----------------------------------\n"

thread=60

constant_seq=/data01/chenjy/TLS_pipeline/ConstantSeq.fa
samtools=/data01/home/chenjy/.conda/envs/pacbio/bin/samtools
cutadapt=/data01/home/chenjy/.conda/envs/pacbio/bin/cutadapt
python3=/data01/home/chenjy/.conda/envs/pacbio/bin/python
bracer=/data01/chenjy/bracer
py_script_dir=/data02/chenjy/pacbio_base
ccs=/data01/home/chenjy/.conda/envs/pacbio/bin/ccs



workdir=${rootdir}/${sample}
mkdir -p ${workdir}
cd ${workdir}

# ################################################################### Step1 call ccs #################################################################################################
echo "-----------------------process Step1 call ccs----------------------------------\n"
## not work for mac!##
mkdir Step1.CCS

# ${ccs} ${subreads}  ${workdir}/Step1.CCS/${sample}.ccs.bam --min-rq 0.95 --num-threads ${thread}

# 转换命令如下:
# ccs --min-passes 3 --min-length 10 --min-rq 0.99 *subreads.bam *.ccs.bam
# **** --min-passes 3: Minimum number of full-length subreads required to generate
# CCS for a ZMW
# **** --min-length 10: Minimum draft length before polishing
# **** --min-rq 0.99:Minimum predicted accuracy in [0, 1]
# **** *subreads.bam:输入文件
# **** *ccs.bam:输出文件

# ## Step 1-2 ##
rm -rf ${workdir}/Step1.CCS/${sample}.ccs.bam
ln -s ${ccs_bam} ${workdir}/Step1.CCS/${sample}.ccs.bam

${samtools} view  ${ccs_bam} | awk '{print ">"$1"\n"$10}' > ${workdir}/Step1.CCS/${sample}.ccs.fa



# ################################################################ step2 calculate percentage of primer #############################################################################
echo "-----------------------process Step2 calculate percentage of primer----------------------------------\n"

mkdir Step2.PrimerDetection

# cutadapt 识别primer
${cutadapt} -g GCCTCTCAGTACGTCAGCAG -j ${thread} --info-file ${workdir}/Step2.PrimerDetection/Pr1.reverse.txt -o ${workdir}/Step2.PrimerDetection/Pr1.reverse.fa ${workdir}/Step1.CCS/${sample}.ccs.fa
${cutadapt} -a CTGCTGACGTACTGAGAGGC -j ${thread} --info-file ${workdir}/Step2.PrimerDetection/Pr1.txt -o ${workdir}/Step2.PrimerDetection/Pr1.fa ${workdir}/Step1.CCS/${sample}.ccs.fa

#rm ${workdir}/Step2.PrimerDetection/Pr1.fa
#rm ${workdir}/Step2.PrimerDetection/Pr1.reverse.fa

# 计算primer的比例
cat ${workdir}/Step2.PrimerDetection/Pr1.reverse.txt | awk -F '\t' '{print $2}' | sort | uniq -c > ${workdir}/Step2.PrimerDetection/Pr1.reverse.count.txt

cat ${workdir}/Step2.PrimerDetection/Pr1.txt | awk -F '\t' '{print $2}' | sort | uniq -c > ${workdir}/Step2.PrimerDetection/Pr1.count.txt


# ############################################################### Step3 detection of constant sequence ###############################################################
echo "-----------------------process Step3 detection of constant sequence----------------------------------\n"


mkdir Step3.ConstantSeq.Detection

${cutadapt} -a file:${constant_seq} -j ${thread} --info-file ${workdir}/Step3.ConstantSeq.Detection/info.txt -o ${workdir}/Step3.ConstantSeq.Detection/tmp.fa ${workdir}/Step1.CCS/${sample}.ccs.fa

#rm ${workdir}/Step3.ConstantSeq.Detection/tmp.fa

cat ${workdir}/Step3.ConstantSeq.Detection/info.txt | awk -F '\t' '{print $2}' | sort | uniq -c > ${workdir}/Step3.ConstantSeq.Detection/ConstantSeq.count.txt

${python3} ${py_script_dir}/Step2.ReverseComplement_pTreads.py -r ${workdir} -s ${sample}

${samtools} view ${workdir}/Step3.ConstantSeq.Detection/${sample}.pA.bam | awk '{print ">"$1"\n"$10}' > ${workdir}/Step3.ConstantSeq.Detection/${sample}.pA.fa

# # 没有constant sequence的序列都可以比对上人的参考基因组



# ######################################################################### Step4.remove primer ####################################################################
echo "-----------------------process Step4 remove primer----------------------------------\n"


# 去除pacbio primer

# 从5‘端开始找 GCCTCTCAGTACGTCAGCAG
# 从3’端开始找 CTGCTGACGTACTGAGAGGC

mkdir Step4.PrimerDetection_2

${cutadapt} -g GCCTCTCAGTACGTCAGCAG -j ${thread} --info-file ${workdir}/Step4.PrimerDetection_2/Pr1.reverse.txt -o ${workdir}/Step4.PrimerDetection_2/Pr1.reverse.fa ${workdir}/Step3.ConstantSeq.Detection/${sample}.pA.fa

${cutadapt} -a CTGCTGACGTACTGAGAGGC -j ${thread} --info-file ${workdir}/Step4.PrimerDetection_2/Pr1.txt -o ${workdir}/Step4.PrimerDetection_2/Pr1.fa ${workdir}/Step3.ConstantSeq.Detection/${sample}.pA.fa

${cutadapt} -a GTCTTAGGAAGACAA -j ${thread} --info-file ${workdir}/Step4.PrimerDetection_2/constant.txt -o ${workdir}/Step4.PrimerDetection_2/constant.fa ${workdir}/Step3.ConstantSeq.Detection/${sample}.pA.fa


# # barcode umi detection
# # 根据primer的位置提取umi， barcode的信息

${python3} ${py_script_dir}/Step3.Barcode.UMI.detection.py -r ${workdir} -s ${sample}




# ############################################################################ Step5.remove polyA ########################################################################
echo "-----------------------process Step5 remove polyA ----------------------------------\n"


mkdir Step5.PolyADetection
${samtools} view ${workdir}/Step4.PrimerDetection_2/${sample}.flnc.bam | awk '{print ">"$1"\n"$10}' > ${workdir}/Step4.PrimerDetection_2/${sample}.flnc.fa


${cutadapt} -g "A{10}" -j ${thread} --info-file ${workdir}/Step5.PolyADetection/pA.txt -o ${workdir}/Step5.PolyADetection/pA.fa ${workdir}/Step4.PrimerDetection_2/${sample}.flnc.fa

${python3} ${py_script_dir}/Step4.rmpA.py -r ${workdir} -s ${sample}




# ############################################################################# Step6.dedup ###############################################################################
echo "-----------------------process Step6 dedup ----------------------------------\n"
mkdir Step6.dedup

${python3} ${py_script_dir}/Step6.dedup.py -r ${workdir} -s ${sample}









