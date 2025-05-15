## 0. Data prepare
1. PacBio HiFi: `potato4.hifi.fastq.gz`
2. ONT ultra-long: `potato4.nanopore.fastq.gz`
3. Hi-C: `potato4.hic_1.fq.gz`, `potato4.hic_2.fq.gz`
## 1. Initial assembly
```shell
hifiasm -o potato4.hifi.ont.asm -t64 --ul /home/pxxiao/project/10_potato_poly/02_asm/01_test/01_hifi_ont/00_data/potato4.nanopore.fastq.gz /home/pxxiao/project/10_potato_poly/02_asm/01_test/01_hifi_ont/00_data/potato4.hifi.fastq.gz

# Two main files are obtained for downstream analysis: primary contig assembly (p_ctg) and phased unitig assembly (p_utg)
```
## 2. Dosage analysis
* By mapping high accuracy HiFi reads, identify collapsed unitigs based on the alignment depth.
1. Obtain hifi alignment depth
```shell
px_shell_dosage.sh -g potato4.hifi.ont.asm.bp.p_utg.gfa.fa -i potato4.hifi.fastq.gz
```
![|400](https://bioin-1320274504.cos.ap-nanjing.myqcloud.com/images/dosage.win10000.jpg)
2. Identify unitig type
* Based on alignment depth, classify unitigs, for example: haplotigs ([0, 45X]), diplotigs ([45X, 74X]), triplotigs ([74X, 103X]), tetraplotigs ([103X, 132X]), and replotigs (>= 132X) for replotigs.
```
python dosage.analysis.contig.type.identified.py --input_file aln.sort.clean.pandepth.win.stat.txt --pandepth

# Main file: contig_depth.txt
```
## 3. mT2T assembly
* By using *p_ctg* assembly, generate a mT2T reference sequence.
1. Remove redundant sequences and obtain preliminary assembly
```shell
phap mt2t --p_ctg potato4.hifi.ont.asm.bp.p_ctg.gfa.fa --threads 2 --process 30 --min_dist 0.2 --match_ratio 0.2

# Main file: retained.fa and removed.fa
```
2. Currently, a **manual adjustment step** is required in the workflow. Based on the alignment results of the `retained.fa` and `removed.fa` against the `DM8` reference genome, users should **inspect overlapping or adjacent contigs** and **manually concatenate them** where appropriate. This step is crucial for generating a continuous and accurate mosaic reference assembly (mT2T.fa).

## 4. Hi-C reads mapping
* Refer to the Hi-C processing method of [HapHiC](https://github.com/zengxiaofei/HapHiC?tab=readme-ov-file#:~:text=Quick%20start-,Align%20Hi%2DC%20data%20to%20the%20assembly,-First%2C%20you%20need): align Hi-C reads to the *p_utg* assembly to obtain the link information between unitigs.
1. Mapping Hi-C reads to *p_utg* assembly
```shell
haphic_shell_data-prepare.sh potato4.hifi.ont.asm.bp.p_utg.gfa.fa potato4.hic_1.fq.gz potato4.hic_2.fq.gz

# Main file: HiC.filtered.bam
```
2. Obatain links between unitigs
```shell
python parse.hic.bam.py --bam HiC.filtered.bam --fasta potato4.hifi.ont.asm.bp.p_utg.gfa.fa --flank 500000 --threads 20

# Main file: full_links.pkl
```

## 5. Allelic unitig table construction, unitigs clustering and re-clustering
```shell
putg=potato4.hifi.ont.asm.bp.p_utg.gfa.fa
threads=36
mt2t=mT2T.fa
contig_type=contig_depth.txt
full_links=full_links.pkl
clm=paired_links.clm

phap cluster \
        --p_utg $p_utg \
        --threads $threads \
        --mT2T /$mt2t \
        --contig_type $contig_type \
        --top_n 4 \
        --chr_num 12 \
        --full_links $full_links \
        --clm $clm > log_cluster_out 2> log_cluster_err &
```

## 6. Phase reads, de novo assembly and scaffolding
* Based on the haplotype clustering results, raw sequencing reads, including HiFi, ONT ultra-long, and Hi-C reads, are phased into haplotype-specific read sets.
* Each haplotype-specific read set is then used to perform de novo genome assembly by using hifiasm parallelly, followed by scaffolding with Hi-C data by using HapHiC to obtain haplotype-resolved chromosome-level assemblies.
1. Mapping the raw sequencing reads (HiFi, ONT ultra-long, and Hi-C) to *p_utg* to get the bam file
2. `phap phase_reads` performs phasing, de novo assembly, and scaffolding
```shell
wd=`pwd`

phap phase_reads phase_reads \
        --bam_hifi $wd/HiFi.clean.bam \
        --bam_hic $wd/HiC.sort.bam \
        --bam_ont $wd/ont.putg.sort.bam \
        --contig_type $wd/contig_depth.txt \
        --hifi $wd/potato4.hifi.fastq.gz \
        --ont $wd/potato4.nanopore.fastq.gz \
        --hic1 $wd/potato4.hic_1.fq.gz \
        --hic2 $wd/potato4.hic_1.fq.gz \
        --threads 10 \
        --process 12 \
        --seed 100 \
        > log_phase_out 2> log_phase_err &
```

## 7. Manual correction by Juicebox
* Juicebox was employed to correct scaffolding errors based on the Hi-C interaction signal map.

## 8. Gap filling
* For each chromosome, [TGS-GapCloser](https://github.com/BGI-Qingdao/TGS-GapCloser) was used to fill gaps using phased ONT reads resulting from the previous step.
* For example, for chromosome 12 haplotype1:
```shell
tgsgapcloser --scaff ../06.review.assembly/out_JBAT_review1.FINAL.change.id.fa --reads ../../../../chr12_group1.ONT.fq.gz --output test_ne --ne
```

## 9. Polishing
* To improve genome quality, [T2T-Automated-Polishing](https://github.com/arangrhie/T2T-Polish) strategy was used for error correction.
1. Merge all the chromosome sequences and obtain `genome.fa`
2. Run T2T-polish
```shell
automated-polishing_v3.sh fullauto -d genome.fa -r potato4.hifi.fastq.gz -s pb -t 64
```
