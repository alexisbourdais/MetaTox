# Metabolites prediction for toxicology

![screenshot](Images/MetaTox.png)

## Overview

*In silico* prediction of a list of molecules whose SMILES code is provided by 5 software packages : **BioTransformer3**, **SyGMa**, **GLORYx**, **MetaTrans** and **Meta-Predictor**.

Biotransformer and Sygma are used via **singularity**. 
Meta-Trans & Meta-Predictor need to clone their github and to create a conda environement. 
GLORYx is used online thanks selenium via a conda environment.
Singularity image downloads and conda environment creations are automated.

As this project was designed for non-bioinformaticians, a graphical interface via zenity was included (optional).

This project has been tested and run on **linux** and **windows-WSL2** (except selenium, which has not yet been tested on wsl2 and zenity may not work on wsl2).

Due to hardware limitations, MetaTrans, Meta-Predictor and selenium may not function correctly. Their use is therefore disabled by default.
You can try running them and seeing the error logs to solve potential problems. Especially for Meta-predictor, which requires cuda drivers.

## Quick start

### Required packages:

- **Singularity** (https://docs.sylabs.io/guides/3.0/user-guide/installation.html) :
  `sudo apt-get install -y singularity-container`
  
- **Conda** (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) :
  `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh; chmod +x Miniconda3-latest-Linux-x86_64.sh; ./Miniconda3-latest-Linux-x86_64.sh`
- May be necessary for metatrans conda env install : `conda config --set channel_priority flexible`
- Some packages needed : `sudo apt install zenity gawk dos2unix` (zenity is optional, gawk and dos2unix are often already installed by default)

### Clone MetaTox, MetaTrans, MetaPredictor directories and configure them: 

- `git clone https://github.com/alexisbourdais/MetaTox; cd MetaTox/; git clone https://github.com/KavrakiLab/MetaTrans; git clone https://github.com/zhukeyun/Meta-Predictor; mkdir Meta-Predictor/prediction; mv Meta-Predictor/model/SoM\ identifier/ Meta-Predictor/model/SoM_identifier; mv Meta-Predictor/model/metabolite\ predictor/ Meta-Predictor/model/metabolite_predictor; chmod +x Meta-Predictor/predict-top15.sh Metatox.sh`
- Download manually the models of Metatrans in https://rice.app.box.com/s/5jeb5pp0a3jjr3jvkakfmck4gi71opo0 and place them in **MetaTrans/models/** (unarchived)

### Run
- Input : Text file with the **molecule ID/name** in the 1st column and the **smile code** in the 2nd column, separated by a **comma**.
- `./Metatox.sh` to activate zenity
- `./MetaTox.sh --input input_file (--option)` to skip zenity

## Parameters

- `./Metatox.sh -h` to see available parameters when zenity was skipped

    REQUIRED parameter

        -i|--input      Input file

    OPTIONAL parameter

        -o|--outdir     Name of the output directory [Results_prediction]

        -p|--predictor  To activate meta-Predictor [No]

        -t|--trans      To activate metaTrans [No]

        -g|--glory      To activate selenium (GLORYx) [No]

        -b|--biotrans   Type of biotransformation to use with BioTransformer3:
                            [allHuman] : Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step 
                            ecbased    : Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism
                            cyp450     : CYP450 metabolism prediction 
                            phaseII    : Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others 
                            hgut       : Human gut microbial
                            superbio   : Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)
                            envimicro  : Environmental microbial

        -n|--nstep      The number of steps for the prediction by BioTransformer3 [default=1]

        -c|--cmode      CYP450 prediction Mode uses by BioTransformer: 
                            1  = CypReact+BioTransformer rules
                            2  = CyProduct only
                           [3] = CypReact+BioTransformer rules+CyProducts
                    
        -1|--phase1     Number of reaction cycles Phase 1 by SygMa [defaut=1]
        -2|--phase2     Number of reaction cycles Phase 2 by SygMa [defaut=1]

        -m|--metabo     Metabolism phase for GloryX : 
                          [All] = Both phase (1 & 2)
                            1   = Phase 1 only
                            2   = Phase 2 only

        -k|--tmp        To keep intermediate files [No] (debugging)

## Documentation

BioTransformer3 : https://bitbucket.org/wishartlab/biotransformer3.0jar/src/master/

SyGMa : https://github.com/3D-e-Chem/sygma

GLORYx : https://nerdd.univie.ac.at/gloryx/

MetaTrans : https://github.com/KavrakiLab/MetaTrans

Meta-Predictor : https://github.com/zhukeyun/Meta-Predictor/tree/main

## Citation

BioTransformer : Djoumbou-Feunang, Y. et al. BioTransformer: a comprehensive computational tool for small molecule metabolism prediction and metabolite identification. J Cheminform 11, 2 (2019)

SyGMa : Ridder, L. & Wagener, M. SyGMa: Combining Expert Knowledge and Empirical Scoring in the Prediction of Metabolites. ChemMedChem 3, 821–832 (2008).

GLORYx : Prediction of the Metabolites Resulting from Phase 1 and Phase 2 Biotransformations of Xenobiotics. (Chem. Res. Toxicol. 2020). Christina de Bruyn Kops, Martin Šícho, Angelica Mazzolari, Johannes Kirchmair

MetaTrans : Litsa, E. E., Das, P. & Kavraki, L. E. Prediction of drug metabolites using neural machine translation. Chem. Sci. 11, 12777–12788 (2020).

MetaPredictor: in silico prediction of drug metabolites based on deep language models with prompt engineering

## Publication

Romain Pelletier; Dina Nahle; Mareme Sarr; Alexis Bourdais; Isabelle Morel; Brendan Le Daré; Thomas Gicquel. Identifying metabolites of new psychoactive substances using in silico prediction tools. Arch Toxicol (2025). https://doi.org/10.1007/s00204-025-04049-5

Pelletier R, Bourdais A, Fabresse N, Ferron PJ, Morel I, Gicquel T, Le Daré B. In silico and in vitro metabolism studies of the new synthetic opiate AP-237 (bucinnazine) using bioinformatics tools. Arch Toxicol. 2024 Jan;98(1):165-179. doi: 10.1007/s00204-023-03617-x. Epub 2023 Oct 15. PMID: 37839054.

Pelletier R, Le Daré B, Le Bouëdec D, Bourdais A, Ferron PJ, Morel I, Porée FH, Gicquel T. Identification, synthesis and quantification of eutylone consumption markers in a chemsex context. Arch Toxicol. 2024 Jan;98(1):151-158. doi: 10.1007/s00204-023-03615-z. Epub 2023 Oct 13. PMID: 37833490.
