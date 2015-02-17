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

export ADAM_DRIVER_MEMORY="55g"
export ADAM_EXECUTOR_MEMORY="55g"
export SPARK_HOME="/root/spark"
export ADAM_OPTS="--conf spark.shuffle.service.enable=true"

# start MR nodes
./ephemeral-hdfs/bin/stop-all.sh
./ephemeral-hdfs/bin/start-all.sh

# make a directory in hdfs
./ephemeral-hdfs/bin/hadoop fs -mkdir .

# download dbsnp file
cd /mnt
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/dbsnp132_20101103.vcf.gz
gunzip dbsnp132_20101103.vcf.gz
mv dbsnp132_20101103.vcf dbsnp_132.vcf
cd ~
./ephemeral-hdfs/bin/hadoop fs -put /mnt/dbsnp_132.vcf .

# pull NA12878 from 1000g
./ephemeral-hdfs/bin/hadoop distcp \
    s3n://1000genomes/data/NA12878/high_coverage_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam \
    ${hdfs_root}/user/${USER}/NA12878.bam

# convert to adam, and remove bam
${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/user/${USER}/NA12878.bam \
    ${hdfs_root}/user/${USER}/NA12878.adam

./ephemeral-hdfs/bin/hadoop fs -rmr \
    ${hdfs_root}/user/${USER}/NA12878.bam

# run flagstat
time ${ADAM_HOME}/bin/adam-submit flagstat \
    ${hdfs_root}/user/${USER}/NA12878.adam

# sort the file
time ${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/user/${USER}/NA12878.adam \
    ${hdfs_root}/user/${USER}/NA12878.sort.adam \
    -sort_reads

./ephemeral-hdfs/bin/hadoop fs -rmr \
    ${hdfs_root}/user/${USER}/NA12878.sort.adam

# mark duplicate reads
time ${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/user/${USER}/NA12878.adam \
    ${hdfs_root}/user/${USER}/NA12878.mkdup.adam \
    -mark_duplicate_reads

./ephemeral-hdfs/bin/hadoop fs -rmr \
    ${hdfs_root}/user/${USER}/NA12878.mkdup.adam

# convert known snps file to adam variants file
${ADAM_HOME}/bin/adam-submit vcf2adam \
    ${hdfs_root}/user/${USER}/dbsnp_132.vcf \
    ${hdfs_root}/user/${USER}/dbsnp_132.var.adam \
    -onlyvariants

./ephemeral-hdfs/bin/hadoop fs -rmr \
    ${hdfs_root}/user/${USER}/dbsnp_132.vcf

# recalibrate quality scores
time ${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/user/${USER}/NA12878.adam \
    ${hdfs_root}/user/${USER}/NA12878.bqsr.adam \
    -recalibrate_base_qualities \
    -known_snps ${hdfs_root}/user/${USER}/dbsnp_132.var.adam

./ephemeral-hdfs/bin/hadoop fs -rmr \
    ${hdfs_root}/user/${USER}/NA12878.bqsr.adam
./ephemeral-hdfs/bin/hadoop fs -rmr \
    ${hdfs_root}/user/${USER}/dbsnp_132.var.adam

# realign indels
time ${ADAM_HOME}/bin/adam-submit transform \
    ${hdfs_root}/user/${USER}/NA12878.adam \
    ${hdfs_root}/user/${USER}/NA12878.ri.adam \
    -realign_indels

./ephemeral-hdfs/bin/hadoop fs -rmr \
    ${hdfs_root}/user/${USER}/NA12878.ri.adam
