#!/bin/bash

# get LCG view
source /cvmfs/sft.cern.ch/lcg/views/LCG_98python3_ATLAS_11/x86_64-centos7-gcc8-opt/setup.sh

# set up python virtual environment
DIR="venv/"
if [ -d "$DIR" ]; then
  ### Take action if $DIR exists ###
  echo "Using python virtual environment in ${DIR}..."
else
  ###  Control will jump here if $DIR does NOT exists ###
  echo "Info: ${DIR} not found. Setting up python virtual environment for the first time..."
  python -m venv venv
  echo "Python virtual environment set up! Activating it now..."
fi

# activate virtual environment
source venv/bin/activate
pip install -r requirements.txt
