#!/usr/bin/env bash
table=species_table.txt
soft_dir=/share/home/xxx/software/sratoolkit.2.9.4-centos_linux64/bin
for Species in $(cat ${table} | cut -d ' ' -f1 )
do
    echo ${Species}
    species=$(grep ${Species} ${table} | head -1 | cut -d " " -f2)
    genome_size=$(grep ${Species} ${table} | head -1 | cut -d " " -f3)

    dir=/share/home/xxx/data/pacbio_h5/${Species}
    ref=${dir}/reference/${species}_pb.fasta
    sra_dir=${dir}/sra
    IP_sra=$(ls ${sra_dir} | sort | head -1)
    IP_SP=${IP_sra:0:10}
    Input_sra=$(ls ${sra_dir} | sort | tail -1)
    Input_SP=${Input_sra:0:10}

    mkdir bwa_in/${Species}
    outdir=/share/home/xxx/data/project/bwa_in/${Species}
	cd ${sra_dir}
	
    ${soft_dir}/fastq-dump --gzip --split-3 -O ${outdir} -A ${IP_sra}
    ${soft_dir}/fastq-dump --gzip --split-3 -O ${outdir} -A ${Input_sra}
	cd ${outdir}
    echo ${species}
	echo ${IP_sra}
    mv ${IP_sra}.fastq.gz ${species}_ip.fastq.gz
	echo ${Input_sra}
    mv ${Input_sra}.fastq.gz ${species}_input.fastq.gz

    bwa index ${ref}
    echo ${ref}
    bwa mem -t 8 ${ref} ${species}_ip.fastq.gz >${species}_ip.sam 
    samtools view -b -S ${species}_ip.sam > ${species}_ip.bam
    samtools sort ${species}_ip.bam ${species}_ip.sort
    samtools index ${species}_ip.sort.bam
    rm -rf ${species}_ip.bam

    bwa mem -t 8 ${ref} ${species}_input.fastq.gz >${species}_input.sam
    samtools view -b -S ${species}_input.sam > ${species}_input.bam
    samtools sort ${species}_input.bam ${species}_input.sort
    samtools index ${species}_input.sort.bam
    rm -rf ${species}_input.bam

    macs2 callpeak -t ${species}_ip.sort.bam -c ${species}_input.sort.bam --nomodel --extsize 100 -g ${genome_size} -n ${Species} -f BAM --verbose 3 --outdir ${outdir}/Macs2_out
done
