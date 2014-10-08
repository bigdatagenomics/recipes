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

from fabric.api import abort, cd, confirm, run

def checkout_adam(fork='bigdatagenomics', branches=[]):
    """
    Makes a fresh checkout of ADAM from github. Applies requested branches. Branch
    list should contain either tuples of two strings where the first string is the
    repository and the second string is the branch, or strings, where the string is
    the branch to merge.
    """

    run('rm -rf adam')
    run('git clone https://www.github.com/%s/adam.git' % fork)
    cd('adam')
    
    # apply pull requests
    run('git config --global --add remote.origin.fetch "+refs/pull/*/head:refs/remotes/origin/pr/*"')
    for branch in branches:
        if isinstance(branch, basestring):
            run('git merge %s' % branch)
        elif isinstance(branch, tuple):
            run('git pull %s %s' % branch)
        elif not confirm('Unsure how to process %s, not a string or tuple. Continue anyways?' % str(branch)):
            abort('Aborting at user request.')

    # cd to top
    cd('..')

def build_adam(skipTests=True):
    """
    Builds ADAM. By default, skips running the unit tests.
    """
    
    cd('adam')
    
    # if we want to skip the unit tests, add flags
    if skipTests:
        skipTests = ' -DskipTests=true'
    else:
        skipTests = ''

    # run build
    run('mvn package%s' % skipTests)

    # cd to top
    cd('..')

def run_adam(arguments):
    """
    Runs ADAM, using the spark submit scripts.
    """
    
    run('./adam/bin/adam-submit %s' % arguments)
