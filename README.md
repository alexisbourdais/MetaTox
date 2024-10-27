# Metabolite prediction to toxicology

## Overview

*In silico* prediction of a list of molecules whose SMILES code is provided by 4 software packages : **BioTransformer3**, **SyGMa**, **MetaTrans** and **Meta-Predictor**.

Biotransformer and Sygma are used via singularity, Meta-Trans & Meta-Predictor need to clone their github.

As this project was designed for non-bioinformaticians, a graphical interface via zenity was included.

## Quick start

### Required packages:

- **Singularity** (https://docs.sylabs.io/guides/3.0/user-guide/installation.html) :
  `sudo apt-get install -y singularity-container`
- **Conda** (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) :
  
  `mkdir -p ~/miniconda3`

  `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh`

  `bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3`

  `rm ~/miniconda3/miniconda.sh`
  
- Necessary for metatrans conda env install : `conda config --set channel_priority flexible`
- Some packages needed : `sudo apt install zenity bc gawk dos2unix`

### Download: 

- `git clone https://github.com/alexisbourdais/MetaTox`
- `cd MetaTox/; git clone https://github.com/KavrakiLab/MetaTrans; git clone https://github.com/zhukeyun/Meta-Predictor`
- `cd Meta-Predictor/; mkdir prediction; cd model/; mv SoM\ identifier/ SoM_identifier; mv metabolite\ predictor/ metabolite_predictor`
- download the models in https://rice.app.box.com/s/5jeb5pp0a3jjr3jvkakfmck4gi71opo0 and place them in **MetaTrans/models/**

### Run
- `chmod +x MetaTox.sh`
- `./MetaTox.sh`

## Documentation

BioTransformer3 : https://bitbucket.org/wishartlab/biotransformer3.0jar/src/master/

SyGMa : https://github.com/3D-e-Chem/sygma

MetaTrans : https://github.com/KavrakiLab/MetaTrans

Meta-Predictor : https://github.com/zhukeyun/Meta-Predictor/tree/main

SdftoSmi & SmitoStr scripts : https://github.com/MunibaFaiza/cheminformatics/tree/main

## Citation

BioTransformer : Djoumbou-Feunang, Y. et al. BioTransformer: a comprehensive computational tool for small molecule metabolism prediction and metabolite identification. J Cheminform 11, 2 (2019)

SyGMa : Ridder, L. & Wagener, M. SyGMa: Combining Expert Knowledge and Empirical Scoring in the Prediction of Metabolites. ChemMedChem 3, 821–832 (2008).

MetaTrans : Litsa, E. E., Das, P. & Kavraki, L. E. Prediction of drug metabolites using neural machine translation. Chem. Sci. 11, 12777–12788 (2020).

MetaPredictor: in silico prediction of drug metabolites based on deep language models with prompt engineering
