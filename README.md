# Metabolites prediction for toxicology

![screenshot](Images/MetaTox.png)

## Overview
*In silico* prediction of molecule metabolites based on their SMILES using 5 software programs:
- **BioTransformer3**
- **SyGMa**
- **GLORYx**
- **MetaTrans**
- **Meta-Predictor**

**BioTransformer3**, **Sygma**, **MetaTrans** and **GloryX (API)** are used via **singularity**. \
**Meta-Predictor** needs to clone its github and to create a **conda** environment. \
Singularity image downloads and conda environment creations are automated (First use may take a long time).

As this project was designed for non-bioinformaticians, a **graphical interface via zenity** was included (**optional**).

This project has been tested and run on **linux** and **windows-WSL2**.

Due to hardware limitations, **Meta-Predictor** (which requires **cuda drivers**) may not function correctly. Its use is therefore disabled by default.
You can try running it and seeing the error logs to solve potential problems.

## Quick start

### Required packages:

- **Singularity** (https://sylabs.io/docs/)
    
- **Conda** (Optional: need for MetaPredictor only) (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html): \
  `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh; chmod +x Miniconda3-latest-Linux-x86_64.sh; ./Miniconda3-latest-Linux-x86_64.sh`

- **APT packages** (zenity is optional, gawk and dos2unix are often already installed by default): \
  `sudo apt install zenity gawk dos2unix`

### Download MetaTox
- `git clone https://github.com/alexisbourdais/MetaTox; chmod +x MetaTox/Metatox.sh`

### Add the sylabs library to use singularity images stored there
- `singularity remote add --no-login SylabsCloud cloud.sycloud.io`

### Download MetaPredictor and configure it
- `cd MetaTox; git clone https://github.com/zhukeyun/Meta-Predictor; mkdir Meta-Predictor/prediction; mv Meta-Predictor/model/SoM\ identifier/ Meta-Predictor/model/SoM_identifier; mv Meta-Predictor/model/metabolite\ predictor/ Meta-Predictor/model/metabolite_predictor; chmod +x Meta-Predictor/predict-top15.sh`

### Run
- Input : Text file with the **molecule ID/name** in the 1st column and the **smile code** in the 2nd column, separated by a **comma**.
- `./Metatox.sh` to activate zenity
- `./MetaTox.sh --input ExempleInput.txt (--option)` to skip zenity

## Parameters
- `./Metatox.sh -h` to see available parameters when zenity was skipped
   
    REQUIRED parameter

        -i|--input   

    OPTIONAL parameter

        -o|--outdir     Name of the output directory [Results_prediction]

        -p|--predictor  To activate meta-Predictor [No]

        -b|--biotrans   Type of biotransformation to use with BioTransformer3:
                            [allHuman] : Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step 
                            ecbased    : Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism
                            cyp450     : CYP450 metabolism prediction 
                            phaseII    : Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others 
                            hgut       : Human gut microbial
                            superbio   : Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)
                            envimicro  : Environmental microbial

        -n|--nstep      The number of steps for the prediction by BioTransformer3 [default=1]

        -c|--cmode      CYP450 prediction Mode uses by BioTransformer3: 
                            1  = CypReact+BioTransformer rules
                            2  = CyProduct only
                           [3] = CypReact+BioTransformer rules+CyProducts
                    
        -1|--phase1     Number of reaction cycles Phase 1 by SygMa [defaut=1]
        -2|--phase2     Number of reaction cycles Phase 2 by SygMa [defaut=1]

        -m|--metabo     Metabolism phase for GloryX : 
                          [phase_1_and_2]
                          phase_1
                          phase_2


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
