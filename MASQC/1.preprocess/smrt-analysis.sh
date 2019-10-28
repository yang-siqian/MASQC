species=Bacillus_subtilis
reference=bsubtilis_pb.fasta
treat=SRR8137526
control=SRR8137527
genome=4078097
Macs2_out=Macs2_out5

export PATH="/share/home/xxx/software/miniconda2/bin:$PATH"

outdir=/share/home/xxx/data/pacbio_h5/Bacillus_subtilis/new/results
datadir=/share/home/xxx/data/pacbio_h5/Bacillus_subtilis/new/data
refdir=/share/home/xxx/data/pacbio_h5/Bacillus_subtilis/new/reference

cd $datadir
mkdir reference
cp $refdir/$reference reference/
cd reference
bwa index $reference
cd ../
bwa aln reference/$reference $treat.fastq > $treat.sai; bwa samse reference/$reference $treat.sai $treat.fastq > $treat.sam
bwa aln reference/$reference $control.fastq > $control.sai; bwa samse reference/$reference $control.sai $control.fastq > $control.sam

samtools view -bS $treat.sam > $treat.bam
samtools view -bS $control.sam > $control.bam

samtools sort $treat.bam $treat.sort
samtools sort $control.bam $control.sort

samtools flagstat $treat.sort.bam > $treat.stat
samtools flagstat $control.sort.bam > $control.stat

mkdir $Macs2_out
macs2 callpeak -t $treat.sort.bam -c $control.sort.bam -n $species --shift -0 --extsize 200 --nomodel -B --SPMR -g $genome --outdir $Macs2_out 2> $Macs2_out/$species.macs2.log

cd $ref_dir
source /share/home/xxx/software/smrtanalysis/install/smrtanalysis_2.3.0.140936/etc/setup.sh
referenceUploader -c -p $refdir -n $species -f $reference --saw 'sawriter -blt 8 -welter'

cd $outdir
ls $refdir/*.bax.h5 > input.fofn
fofnToSmrtpipeInput.py input.fofn > input.xml
smrtpipe.py -D NPROC=20 --params=params_sb.xml xml:input.xml > smrtpipe.log
