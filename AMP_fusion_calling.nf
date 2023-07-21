#!/usr/bin/env nextflow

params.reads = "fastq/*_R{1,2}.fastq.gz"

params.adapters = "refs/adapters.fa"
adapters = file(params.adapters)

params.umi_custom_scipt = "scripts/extract_archer_umi.py"
umi_script = file(params.umi_custom_scipt)

params.star_index = "refs/GRCh37_gencode_v19_CTAT_lib_Oct012019.source/ctat_genome_lib_build_dir/ref_genome.fa.star.idx"
star_index = file(params.star_index)

params.ref_gtf = "refs/GRCh37_gencode_v19_CTAT_lib_Oct012019.source/ctat_genome_lib_build_dir/ref_annot.gtf"
prot_gtf = file(params.ref_gtf)

params.star_fus = "refs/GRCh37_gencode_v19_CTAT_lib_Oct012019.source/ctat_genome_lib_build_dir"
star_fus_refs = file(params.star_fus)

params.uscript = "/scripts/tidy_bam.sh"
uniq_script = file(params.uscript)


Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs }


read_pairs.into{umi_tools_input} 

/* read_pairs.into{mapping_input} */

process umi_extract {

    publishDir "results/mark_umi"
    cpus 4
    input:
    set pair_id, file(reads) from umi_tools_input
    file "extract_archer_umi.py" from umi_script

    output:
    set val(pair_id), file('*_R{1,2}_marked.fastq.gz') into trim_input

    script:
    """
    python extract_archer_umi.py -f1 ${reads[0]} -f2 ${reads[1]} -o ${pair_id}

    """
} 


process trim_reads {

    publishDir "results/trimmed"
    cpus 4

    input:
    set pair_id, file(reads) from trim_input
    file "adapters.fa" from adapters

    output:
    set val(pair_id), file('*_R{1,2}_trimmed.fastq.gz') into mapping_input

    script:
     """
     bbduk.sh in1=${reads[0]} \
     in2=${reads[1]} \
     out1="${pair_id}_R1_trimmed.fastq.gz" \
     out2="${pair_id}_R2_trimmed.fastq.gz" \
     ref=adapters.fa \
     minlength=30 \
     outs="${pair_id}_singletons.fastq.gz" \
     ktrim=r \
     k=19 \
     tbo=t
     """
}

mapping_input.into{star_input; star_fusion_pass}


process staralign {
        publishDir "results/star"
        cpus 16
        input:
    set pair_id, file(reads) from star_input

    file "*"  from prot_gtf
    file "*" from star_index

    output:
    set val(pair_id), file("*{.bam,.bai}") into bam_file

    script:
    """

    STAR --readFilesCommand zcat \
    --readFilesIn  $reads \
    --genomeDir ref_genome.fa.star.idx \
    --runThreadN 16 \
    --outFileNamePrefix ${pair_id} \
    --alignIntronMax 100000 \
    --alignMatesGapMax 100000 \
    --alignSJDBoverhangMin 10 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimJunctionOverhangMin 12 \
    --chimMultimapScoreRange 3 \
    --chimNonchimScoreDropMin 10 \
    --chimScoreDropMax 20 \
    --chimScoreJunctionNonGTAG -4 \
    --chimScoreSeparation 10 \
    --chimSegmentMin 12 \
    --chimSegmentReadGapMax 7\
    --outFilterIntronMotifs None \
    --outFilterMismatchNmax 10 \
    --outFilterMultimapNmax 10 \
    --outFilterType Normal \
    --outSAMmapqUnique 255 \
    --outSAMstrandField intronMotif \
    --outSAMtlen 1 \
    --outSAMunmapped Within \
    --outSJfilterCountTotalMin 3 1 1 1 \
    --outSJfilterCountUniqueMin 3 1 1 1 \
    --peOverlapMMp 0.1 \
    --peOverlapNbasesMin 12 \
    --sjdbScore 2  \
    --twopassMode Basic \
    --quantMode GeneCounts \
    --chimOutJunctionFormat 1 \
    --outReadsUnmapped None \
    --chimOutType WithinBAM \
    --outSAMtype 'BAM' 'SortedByCoordinate'

    mv ${pair_id}Aligned.sortedByCoord.out.bam  ${pair_id}.bam

    samtools index ${pair_id}.bam

    mv ${pair_id}.bam.bai ${pair_id}.bai

    """
}


/*
 * step 4.1 more qc genebody cov, clipping profile etc
 */



 process umi_dedup {
    publishDir = "results/dedup"
    cpus 4
    input:
    set pair_id, file(bam) from bam_file
    file("tidy_bam.sh") from uniq_script
    output:

    set val(pair_id), file("*{_dedup.bam,_dedup.bai}") into mol_con

    shell:
    '''
    touch !{bam[0]}
    umi_tools dedup \
    -I !{bam[1]} \
    --umi-separator='_' \
    --edit-distance-threshold=2 \
    --paired \
     -S dedup_tmp.bam \
    --unpaired-reads='discard' \
    --spliced-is-unique \
    --output-stats='stats.txt' \
    --multimapping-detection-method NH

    ./tidy_bam.sh

    picard FixMateInformation I=no_sing.bam O=!{pair_id}"_dedup.bam"
    samtools index !{pair_id}"_dedup.bam"
    mv !{pair_id}"_dedup.bam.bai" !{pair_id}"_dedup.bai"



    '''
}


process uniq_reads {
    publishDir = "results/unique_reads"
    cpus  4
    input:
        set pair_id, file(bam) from mol_con
        file reads from star_fusion_pass.flatten().collect()

    output:
        set val(pair_id), file("*_uniq.R{1,2}.fastq.gz") into starfusion_input
    shell:
    '''
    samtools view !{bam[1]} | awk '{print $1}' > names
    filterbyname.sh in1=!{pair_id}_R1_trimmed.fastq.gz in2=!{pair_id}_R2_trimmed.fastq.gz out1=!{pair_id}_uniq.R1.fastq.gz out2=!{pair_id}_uniq.R2.fastq.gz names=names retain=t

    '''
}


process star_fusion {
    publishDir = "results/star_fusion"
    cpus  16
    input:
        set pair_id, file(reads) from starfusion_input
        file "*" from star_fus_refs
    output:
        set val(pair_id), file("*_results")
    script:
    """
    mkdir ${pair_id}_results

    STAR-Fusion --left_fq ${pair_id}_uniq.R1.fastq.gz \
    --right_fq ${pair_id}_uniq.R2.fastq.gz \
    --genome_lib_dir `pwd`/ctat_genome_lib_build_dir \
    --output_dir ${pair_id}_results \
    --CPU 16 \
    --require_LDAS 0 \
    --no_remove_dups

    """
}

