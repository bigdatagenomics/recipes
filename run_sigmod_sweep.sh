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

cd ~

wget http://d3kbcqa49mib13.cloudfront.net/spark-1.2.0-bin-hadoop2.4.tgz
tar xzvf spark-1.2.0-bin-hadoop2.4.tgz
export SPARK_HOME="~/spark-1.2.0-bin-hadoop2.4"

for n in 128 64 32
do
    fab provision:${n} |& tee provision_${n}
    fab -i ~/${EC2_PRIVATE_KEY_FILE} configure_master bake:adam-transforms |& tee run_${n}
    fab teardown |& tee teardown_${n}
done
