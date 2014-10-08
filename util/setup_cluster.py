"""
Licensed to Big Data Genomics (BDG) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The BDG licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from fabric.api import abort, env, local
import os
import os.path
import time

# a dictionary containing a map between the cluster name and master node url
# for all clusters in current use
current_clusters = {}

def setup_cluster(prefix, nodes, keypair_path, zone, instance_type, spot_price=None, maven_version='3.2.3'):
    """
    Sets up a Spark cluster on EC2 using the Spark EC2 scripts.
    """

    # get path to spark
    if 'SPARK_HOME' in os.environ:
        spark_home = os.environ['SPARK_HOME']
    else:
        abort("Couldn't find SPARK_HOME.")

    # get the name of the keypair
    keypair = os.path.basename(keypair_path)[:-4]
        
    # set spot price if provided
    if spot_price != None:
        spot_price = ' --spot-price=%.2f' % spot_price
    else:
        spot_price = ''

    # build cluster name
    name = prefix + time.time()

    # set up cluster
    local('%s/ec2/spark-ec2 -k %s -i %s -s %d --zone %s launch %s --instance_type=%s%s' % (spark_home,
                                                                                           keypair,
                                                                                           keypair_path,
                                                                                           zone,
                                                                                           name,
                                                                                           instance_type,
                                                                                           spot_price))

    # set master node details
    master = local('%s/ec2/spark-ec2 get-master %s' % (spark_home, name))
    env.host.append(master)
    current_clusters[name] = master
    env.user = 'root'
    env.key_filename = keypair_path

    # install git
    sudo('yum install git -y')

    # install maven binaries from the apache site
    cd('/usr/local/bin')
    run('wget http://apache.mesi.com.ar/maven/maven-3/%s/binaries/apache-maven-%s-bin.tar.gz' % (maven_version, maven_version))
    run('tar xzvf apache-maven-%s-bin.tar.gz' % maven_version)
    run('cat "export M2_HOME=/usr/local/bin/apache-maven-%s-bin" > ~/.profile' % maven_version)
    run('cat "export M2=${M2_HOME}/bin" > ~/.profile')
    run('cat "export PATH=$M2:$PATH" > ~/.profile')
    cd('~')

    # set spark home
    run('cat "export SPARK_HOME=/root/spark" > ~/.profile')

    run('. ~/.profile')

def teardown_cluster(prefix):
    """
    Tears down a Spark cluster on EC2 using the Spark EC2 scripts.
    """

    # get path to spark
    if 'SPARK_HOME' in os.environ:
        spark_home = os.environ['SPARK_HOME']
    else:
        abort("Couldn't find SPARK_HOME.")

    # tear down cluster
    local('%s/ec2/spark-ec2 destroy %s' % (spark_home,
                                           name))

def reboot_cluster(workers=None):
    """
    Stops the cluster and restarts it. Has an optional workers parameter; if
    this parameter is defined, we only start that many workers.
    """
    
    # stop the cluster
    run('./spark/sbin/stop-all.sh')
    
    # read the worker node URLs
    workerUrls = run('cat ./spark/conf/slaves').split()

    # if we have a defined number of workers, trim the worker list down
    if isinstance(workers, int) and workers <= len(workerUrls):
        workerUrls = workerUrls[0:workers - 1]
    elif (isinstance(workers, int) and
          not confirm('Requested more workers than are running. Continue?')):
        abort('Aborting because not enough workers are in cluster.')

    # start all worker nodes
    for worker in workerUrls:
        run('./spark/sbin/start-slave.sh %s' % worker)

    # start the leader
    run('./spark/sbin/start-master.sh')
