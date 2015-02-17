#! /bin/bash
#
# Licensed to Big Data Genomics (BDG) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The BDG licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set -x

export ADAM_DRIVER_MEMORY="200g"
export ADAM_EXECUTOR_MEMORY="200g"
export SPARK_HOME="/home/ubuntu/spark"
export ADAM_HOME="/home/ubuntu/adam"
export ADAM_OPTS="--conf spark.shuffle.service.enable=true"

cd /mnt2/

# pull NA12878 from 1000g
~/s3cmd/s3cmd get \
    --secret_key=${SECRET_KEY} \
    --access_key=${ACCESS_KEY} \
    s3://1000genomes/data/NA12878/high_coverage_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam
mv NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam NA12878.bam

# get reference genome                                                                                                                                                                                        
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.fai

cd ~

# index bam file
samtools index /mnt2/NA12878.bam

# index reference genome
samtools faidx /mnt2/human_g1k_v37.fasta

# create sequence dictionary for reference genome
java -Xmx200g \
    -jar ~/picard/dist/picard.jar \
    CreateSequenceDictionary \
    REFERENCE=/mnt2/human_g1k_v37.fasta \
    OUTPUT=/mnt2/human_g1k_v37.dict 

# create target indels
time java -Xmx200g \
    -jar ~/gatk-protected/target/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /mnt2/human_g1k_v37.fasta \
    -I /mnt2/NA12878.bam \
    -o /mnt2/target_intervals.list \
    -nt 32

# realign reads
time java -Xmx200g \
    -jar ~/gatk-protected/target/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R /mnt2/human_g1k_v37.fasta \
    -targetIntervals /mnt2/target_intervals.list \
    -I /mnt2/NA12878.bam \
    -o /mnt2/NA12878.realigned.bam \
    -nt 32
rm -r /mnt2/NA12878.realigned.bam

# convert to adam, and remove bam
${ADAM_HOME}/bin/adam-submit transform \
    /mnt2/NA12878.bam \
    /mnt2/NA12878.adam
rm -r /mnt2/NA12878.bam

# convert known snps file to adam variants file
${ADAM_HOME}/bin/adam-submit vcf2adam \
    /mnt2/dbsnp_138.vcf \
    /mnt2/dbsnp_138.var.adam \
    -onlyvariants
rm /mnt2/dbsnp_138.vcf

# recalibrate quality scores
time ${ADAM_HOME}/bin/adam-submit transform \
    /mnt2/NA12878.adam \
    /mnt2/NA12878.bqsr.adam \
    -recalibrate_base_qualities \
    -known_snps /mnt2/dbsnp_138.var.adam

