## Scripts for FTAG algorithm studies

This project contains scripts for flavour tagging algorithm studies.

### Setting up the environment

You can choose:

- set up a python environment which supports a virtual environment for local installations of packages and features modern python 3: `prepare_environment_python.sh`
- set up an ATLAS analysis release `AnalysisBase` environment with python2 and limited support for external packages but access to ATLAS tools and the xAOD EDM: `prepare_environment_atlas.sh`

### Scripts

#### Creating symlinks to files on the LOCALGROUPDISK

If you have replicated h5 ntuples to the local group disk of your institute using the [`rucio`](https://rucio-ui.cern.ch) interface, you can create symbolic links to access them as if you had local access to them.

Update the sample list with ntuples in `data/samplelists` and execute:
```
bash createLinks.sh
```

#### Extract flavour tagging efficiencies from CDI file

Using the ATLAS analysis release environment, process the CDI file with the [BTaggingEfficiencyTool](https://gitlab.cern.ch/atlas/athena/-/blob/21.2/PhysicsAnalysis/JetTagging/JetTagPerformanceCalibration/xAODBTaggingEfficiency/xAODBTaggingEfficiency/BTaggingEfficiencyTool.h) to dump the efficiencies in a `csv` file.

You can define the binning in `data/binning.yml`. The efficiency is evaluated in 2D using `pt` and `abs(eta)`at the center of the bins. The output is provided in `data/cdi_dump.csv`.

```
source prepare_environment_atlas.sh
python processCDI.py
```

#### Estimate and plot efficiency from h5 ntuple and compare with CDI file as reference

Using the python 3 environment, execute the `plotEfficiency.py` script. It extracts b-jets from the h5 ntuple, calculates the b-tagging discriminant `D_b` using the stored tagger outputs `DL1r_pb`, `DL1r_pc`, and `DL1r_pu`. It further plots the obtained efficiency and compares it to that stored in `data/cdi_dump.csv` as 2D maps.

```
source prepare_environment_python.sh
python plotEfficiency.py
```
