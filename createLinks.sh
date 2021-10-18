#!/bin/bash

# create symlinks to localgroupdisk storage
SAMPLELIST=data/samplelists/ntuples-2021-10-15.txt
TARGETDIR=/nfs/dust/atlas/user/pgadow/ftag/data/ntuple_links/

# set up rucio
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup rucio
voms-proxy-init -voms atlas

python createLinksToSRMFiles.py --sampleList ${SAMPLELIST} --scope user.alfroch -o ${TARGETDIR}
