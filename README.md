# Metabolite prediction to toxicology

## Overview

- ***InSilico***: *in silico* prediction of a list of molecules whose SMILES code is provided by 3 software packages (BioTransformer3, SyGMa, MetaTrans).

Linux and Windows (WSL) compatible. Not tested on Mac

## Quick start

### Required packages:
- **Zenity** `sudo apt install zenity` (As this project was designed for non-bioinformaticians, a graphical interface via zenity was included. However, modules can be used separately)
- **bc command** `sudo apt install bc`
- **gawk** `sudo apt install gawk`
- **Java** `sudo apt install default-jre`
- **dos2unix** `sudo apt install dos2unix`
- **Singularity** (https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
- **Conda** (https://github.com/conda/conda).
- Necessary for metatrans conda environment installation : `conda config --set channel_priority flexible`

### Download: 
- `git clone https://github.com/alexisbourdais/MolecularNetworking`
- `cd Tools/`
-  **Biotranformer3.0** `git clone https://bitbucket.org/wishartlab/biotransformer3.0jar.git` 
- **MetaTrans** `git clone https://github.com/KavrakiLab/MetaTrans` 
- download the models in https://rice.app.box.com/s/5jeb5pp0a3jjr3jvkakfmck4gi71opo0 and place them in **MetaTrans/models/**

### Run
- `cd ..`
- `chmod +x Main.sh`
- `./Main.sh`

To use the ***in silico*** mode, create a text file with each line = namemolecule,smilecode

## Sources & Documentations

BioTransformer3 : https://bitbucket.org/wishartlab/biotransformer3.0jar/src/master/

SyGMa : https://github.com/3D-e-Chem/sygma

MetaTrans : https://github.com/KavrakiLab/MetaTrans

SdftoSmi & SmitoStr scripts : https://github.com/MunibaFaiza/cheminformatics/tree/main

## References

Djoumbou-Feunang, Y. et al. BioTransformer: a comprehensive computational tool for small molecule metabolism prediction and metabolite identification. J Cheminform 11, 2 (2019)

Ridder, L. & Wagener, M. SyGMa: Combining Expert Knowledge and Empirical Scoring in the Prediction of Metabolites. ChemMedChem 3, 821–832 (2008).

Litsa, E. E., Das, P. & Kavraki, L. E. Prediction of drug metabolites using neural machine translation. Chem. Sci. 11, 12777–12788 (2020).

Chambers, M. C. et al. A cross-platform toolkit for mass spectrometry and proteomics. Nat Biotechnol 30, 918–920 (2012).
