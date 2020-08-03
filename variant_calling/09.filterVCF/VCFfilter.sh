(# SNP
java -Xmx3g \
    -Djava.io.tmpdir=${WP}/tmp \
    -jar ${WP}/Install/GenomeAnalysisTK.jar \
    -R ${WP}/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \
    -T SelectVariants \
    -V ../08.mergeGVCF/$CHR.$ID.raw.vcf.gz \
    -selectType SNP \
    -o $CHR.$ID.snp.vcf.gz

java -Xmx3g \
    -Djava.io.tmpdir=${WP}/tmp \
    -jar ${WP}/Install/GenomeAnalysisTK.jar \
    -R ${WP}/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \
    -T VariantFiltration \
    -V $CHR.$ID.snp.vcf.gz \
    --logging_level ERROR \
    --filterExpression " QD < 2.0 " \
    --filterName "QD_filter" \
    --filterExpression " FS > 60.0 " \
    --filterName "FS_filter" \
    --filterExpression " MQRankSum < -12.5 " \
    --filterName "MQRankSum_filter" \
    --filterExpression " ReadPosRankSum < -8.0 " \
    --filterName "ReadPosRankSum_filter" \
    --filterExpression " SOR > 3.0 " \
    --filterName "SOR_filter" \
    --filterExpression " MQ < 40.0 " \
    --filterName "MQ_filter" \
    --genotypeFilterExpression " DP > 30 || DP < 3 " \
    --genotypeFilterName "DP_filter" \
    --setFilteredGtToNocall \
    --clusterSize 3 \
    --clusterWindowSize 10 \
    -o $CHR.$ID.snp.filter.vcf.gz

java -Xmx3g \
    -Djava.io.tmpdir=${WP}/tmp \
    -jar ${WP}/Install/GenomeAnalysisTK.jar \
    -R ${WP}/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \
    -T SelectVariants \
    -V $CHR.$ID.snp.filter.vcf.gz \
    --excludeFiltered \
    --excludeNonVariants \
    --removeUnusedAlternates \
    -o $CHR.$ID.snp.filter.final.vcf.gz

# rm $CHR.$ID.snp.vcf.gz* $CHR.$ID.snp.filter.vcf.gz*

# INDEL
java -Xmx3g \
    -Djava.io.tmpdir=${WP}/tmp \
    -jar ${WP}/Install/GenomeAnalysisTK.jar \
    -R ${WP}/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \
    -T SelectVariants \
    -V ../08.mergeGVCF/$CHR.$ID.raw.vcf.gz \
    -selectType INDEL \
    -o $CHR.$ID.indel.vcf.gz

java -Xmx3g \
    -Djava.io.tmpdir=${WP}/tmp \
    -jar ${WP}/Install/GenomeAnalysisTK.jar \
    -R ${WP}/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \
    -T VariantFiltration \
    -V $CHR.$ID.indel.vcf.gz \
    --logging_level ERROR \
    --filterExpression " QD < 2.0  " \
    --filterName "QD_filter" \
    --filterExpression " FS > 200.0 " \
    --filterName "FS_filter" \
    --filterExpression " ReadPosRankSum < -20.0 " \
    --filterName "ReadPosRankSum_filter" \
    --genotypeFilterExpression " DP > 30 || DP < 3 " \
    --genotypeFilterName "DP_filter" \
    --setFilteredGtToNocall \
    -o $CHR.$ID.indel.filter.vcf.gz

java -Xmx3g \
    -Djava.io.tmpdir=${WP}/tmp \
    -jar ${WP}/Install/GenomeAnalysisTK.jar \
    -R ${WP}/genome/IWGSCv1_MP/WholeGenomeSplit_MP.fa \
    -T SelectVariants \
    -V $CHR.$ID.indel.filter.vcf.gz \
    --excludeFiltered \
    --excludeNonVariants \
    --removeUnusedAlternates \
    -o $CHR.$ID.indel.filter.final.vcf.gz

# rm $CHR.$ID.indel.vcf.gz* $CHR.$ID.indel.filter.vcf.gz* 
) &
#
