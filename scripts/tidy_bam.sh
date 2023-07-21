samtools view -h dedup_tmp.bam | awk '!seen[$0]++' > tmp1.sam

samtools view -h -S -b tmp1.sam > tmp1.bam

samtools view -h -@ 8 -F 0x08 -b tmp1.bam > no_sing.bam


#java -jar /data/TPT_Test/SOUP/picard.jar FixMateInformation I=no_sing.bam O=!{pair_id}"_dedup.bam"
#samtools index !{pair_id}"_dedup.bam"
#mv !{pair_id}"_dedup.bam.bai" !{pair_id}"_dedup.bai"
