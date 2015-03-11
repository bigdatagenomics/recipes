#! /usr/bin/env python
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

import os
from fabric.api import task, local, run, require, hosts, env, open_shell, shell_env, prefix, cd, sudo
from fabric.colors import green
from fabric.tasks import execute

if not env.key_filename:
    env.key_filename = os.environ['EC2_PRIVATE_KEY_FILE']

@task
def provision(size, type_='r3.2xlarge'):
    cmd = ' '.join([
        '{SPARK_HOME}/ec2/spark-ec2'.format(**os.environ),
        '-k {EC2_KEY_PAIR}'.format(**os.environ),
        '-i {EC2_PRIVATE_KEY_FILE}'.format(**os.environ),
        '-s {size}'.format(size=size),
        '-t {type_}'.format(type_=type_),
        '--delete-groups',
        '--copy-aws-credentials',
        'launch bdg-recipies'])
    return local(cmd)

def get_set_master_host():
    cmd = ' '.join([
        '{SPARK_HOME}/ec2/spark-ec2'.format(**os.environ),
        '-k {EC2_KEY_PAIR}'.format(**os.environ),
        '-i {EC2_PRIVATE_KEY_FILE}'.format(**os.environ),
        'get-master bdg-recipies'])
    result = local(cmd, capture=True)
    host = result.split('\n')[2].strip()
    env.user = 'root'
    env.hosts = [host]
    print(green(host))
    return host

@task
def login():
    if env.hosts == None:
        hosts = get_set_master_host()
        execute(open_shell, hosts=hosts)
    else:
        execute(open_shell)

@task
def teardown():
    cmd = ' '.join([
        '{SPARK_HOME}/ec2/spark-ec2'.format(**os.environ),
        '-k {EC2_KEY_PAIR}'.format(**os.environ),
        '-i {EC2_PRIVATE_KEY_FILE}'.format(**os.environ),
        'destroy bdg-recipes'])
    local(cmd)

@task
def _configure_master_yum():
    # install fabric and luigi
    with cd('/tmp'):
        # protobuf for luigi
        run('yum install -y protobuf protobuf-devel protobuf-python')
        run('curl -O https://bootstrap.pypa.io/get-pip.py')
        run('python get-pip.py')
        run('pip install fabric luigi')
    # check out adam
    with cd('~'):
        run('git clone https://www.github.com/bigdatagenomics/adam.git')
    # check out bdg recipes
    with cd('~'):
        run('git clone https://www.github.com/bigdatagenomics/bdg-recipes.git')
    with cd('~/bdg-recipes'):
        run('git pull https://www.github.com/fnothaft/bdg-recipes.git sigmod')
    # download maven
    with cd('~'):
        run('wget http://supergsego.com/apache/maven/maven-3/3.2.5/binaries/apache-maven-3.2.5-bin.tar.gz')
        run('tar xzvf apache-maven-3.2.5-bin.tar.gz')
    # build adam
    with cd('adam'), shell_env(MAVEN_OPTS='-Xmx512m -XX:MaxPermSize=128m',
                               MAVEN_HOME='~/apache-maven-3.2.5'):
        run('~/apache-maven-3.2.5/bin/mvn clean package -DskipTests')

@task
def _configure_master_aptitude(spark_ver="1.1.1"):
    sudo('''debconf-set-selections <<< "postfix postfix/main_mailer_type string 'No Configuration'"''')
    sudo('''debconf-set-selections <<< "postfix postfix/mailname string your.hostname.com"''')
    sudo('apt-get update')
    sudo('apt-get install git openjdk-7-jdk python-pip libprotobuf-dev libprotobuf-dev mdadm llvm ant docker.io python-dev python-dateutil -y')
    sudo('pip install fabric luigi')

    # check out s3cmd
    with cd('~'):
        run('git clone https://www.github.com/s3tools/s3cmd.git')
    # check out adam
    with cd('~'):
        run('git clone https://www.github.com/bigdatagenomics/adam.git')
    # check out bdg recipes
    with cd('~'):
        run('git clone https://www.github.com/bigdatagenomics/bdg-recipes.git')
    with cd('~/bdg-recipes'):
        run('git pull https://www.github.com/fnothaft/bdg-recipes.git sigmod')
    # download maven
    with cd('~'):
        run('wget http://supergsego.com/apache/maven/maven-3/3.2.5/binaries/apache-maven-3.2.5-bin.tar.gz')
        run('tar xzvf apache-maven-3.2.5-bin.tar.gz')
    # build adam
    with cd('adam'), shell_env(MAVEN_OPTS='-Xmx512m -XX:MaxPermSize=128m',
                               MAVEN_HOME='~/apache-maven-3.2.5'):
        run('~/apache-maven-3.2.5/bin/mvn clean package -DskipTests')
    # checkout samblaster
    with cd('~'):
        run('git clone https://www.github.com/gregoryfaust/samblaster.git')
    # install samblaster
    with cd('samblaster'):
        run('make')
    # checkout sambamba docker container and then build
    with cd('~'):
        sudo('docker pull lomereiter/centos-ldc')
        run('touch build-sambamba')
        run('echo "git clone --recursive https://github.com/lomereiter/sambamba.git" >> build-sambamba')
        run('echo "cd sambamba && make sambamba-ldmd2-64" >> build-sambamba')
        run('echo "exit" >> build-sambamba')
        sudo('docker run -i lomereiter/centos-ldc /bin/bash < build-sambamba')
        sudo('docker cp `docker ps -l -q`:/sambamba/build/sambamba ~')
    # install samtools
    sudo('apt-get install samtools -y')
    # check out picard
    with cd('~'):
        run('git clone https://www.github.com/broadinstitute/picard.git')
    # checkout htsjdk...
    with cd('~/picard'):
        run('git clone https://www.github.com/samtools/htsjdk.git')
    # install picard
    with cd('picard'):
        run('ant -lib lib/ant package-commands')
    # check out the gatk
    with cd('~'):
        run('git clone https://github.com/broadgsa/gatk-protected.git')
    # build the gatk
    with cd('gatk-protected'), shell_env(MAVEN_OPTS='-Xmx512m -XX:MaxPermSize=128m',
                                         MAVEN_HOME='~/apache-maven-3.2.5'):
        run('~/apache-maven-3.2.5/bin/mvn package -DskipTests')
    # pull down spark
    if spark_ver is "1.2.1":
        with cd('~'):
            run('wget http://d3kbcqa49mib13.cloudfront.net/spark-1.2.1-bin-hadoop2.4.tgz')
            run('tar xzvf spark-1.2.1-bin-hadoop2.4.tgz')
            run('mv spark-1.2.1-bin-hadoop2.4 spark')
    else:
        with cd('~'):
            run('wget http://d3kbcqa49mib13.cloudfront.net/spark-1.1.1-bin-hadoop2.4.tgz')
            run('tar xzvf spark-1.1.1-bin-hadoop2.4.tgz')
            run('mv spark-1.1.1-bin-hadoop2.4 spark')

    # make raid array
    with cd('~'):
        output = run('ls /dev/xvd* | grep -v xvda')
        files = output.split()
        sudo('mdadm --create /dev/md0 --level=0 --raid-devices=%d %s' % (len(files),
                                                                         ' '.join(files)))
        sudo('mkfs -t ext4 /dev/md0')
        sudo('mkdir /mnt2/')
        sudo('mount /dev/md0 /mnt2/')
        sudo('chmod ugoa+rw /mnt2/')

@task
def configure_master():
    if env.hosts == None:
        hosts = get_set_master_host()
        execute(_configure_master_yum, hosts=hosts)
    else:
        execute(_configure_master_aptitude)

@task
def _remote_bake(recipe):
    with prefix('source ~/bdg-recipes/bdg-recipes-ec2-variables.sh'), shell_env(ACCESS_KEY=os.environ["AWS_ACCESS_KEY_ID"],
                                                                  SECRET_KEY=os.environ["AWS_SECRET_ACCESS_KEY"]):
        run('. ~/bdg-recipes/%s/run.sh' % recipe)

@task
def bake(recipe):
    if env.hosts == None:
        hosts = get_set_master_host()
        execute(_remote_bake, recipe, hosts=hosts)    
    else:
        execute(_remote_bake, recipe)
