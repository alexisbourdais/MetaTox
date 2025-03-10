# Metabolites prediction for toxicology

**The purpose of this branch is simply to keep an old version for reference. See the main branch for the latest version**

## Quick start

### Required packages:

- **Singularity** (https://docs.sylabs.io/guides/3.0/user-guide/installation.html) :
  `sudo apt-get install -y singularity-container`
- **Conda** (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) :
  
  `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh; chmod +x Miniconda3-latest-Linux-x86_64.sh; ./Miniconda3-latest-Linux-x86_64.sh`
- Necessary for metatrans conda env install : `conda config --set channel_priority flexible`
- Some packages needed : `sudo apt install gawk dos2unix csvkit`

### Download project, MetaTrans-MetaPredictor directory and configure them: 

- `git clone https://github.com/alexisbourdais/MetaTox; cd MetaTox/Tools; git clone https://github.com/KavrakiLab/MetaTrans; git clone https://bitbucket.org/wishartlab/biotransformer3.0jar.git`
- download the models in https://rice.app.box.com/s/5jeb5pp0a3jjr3jvkakfmck4gi71opo0 and place them in **MetaTrans/models/** (unarchived)

### Run
- Input : Text file with the **molecule ID or name** in the 1st column and the **smile code** in the 2nd column, separated by a **comma**.
- `./Prediction_Metabo.sh --input input_file`

## Parameters

- `./Metatox.sh -h` to see available parameters when zenity was skipped

      REQUIRED parameter

        -i|--input
  
      OPTIONAL parameter

        BioTransformer3
        -t|--type       Type of biotransformation to use with BioTransformer3:
                            [allHuman] : Predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjugation) at each step 
                            ecbased    : Prediction of promiscuous metabolism (e.g. glycerolipid metabolism). EC-based metabolism is also called Enzyme Commission based metabolism
                            cyp450     : CYP450 metabolism prediction 
                            phaseII    : Prediction of major conjugative reactions, including glucuronidation, sulfation, glycine transfer, N-acetyl transfer, and glutathione transfer, among others 
                            hgut       : Human gut microbial
                            superbio   : Runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.)
                            envimicro  : Environmental microbial
        MetaTrans
        -b, --beam      Top métabolites [defaut = 5]
        -m, --min       Taille minimum des métabolites (en SMILE) [defaut = 5]
        -M, --max       Taille maximum des métabolites (en SMILE) [defaut = 120]

## Documentation

BioTransformer3 : https://bitbucket.org/wishartlab/biotransformer3.0jar/src/master/

SyGMa : https://github.com/3D-e-Chem/sygma

MetaTrans : https://github.com/KavrakiLab/MetaTrans

## Citation

BioTransformer : Djoumbou-Feunang, Y. et al. BioTransformer: a comprehensive computational tool for small molecule metabolism prediction and metabolite identification. J Cheminform 11, 2 (2019)

SyGMa : Ridder, L. & Wagener, M. SyGMa: Combining Expert Knowledge and Empirical Scoring in the Prediction of Metabolites. ChemMedChem 3, 821–832 (2008).

MetaTrans : Litsa, E. E., Das, P. & Kavraki, L. E. Prediction of drug metabolites using neural machine translation. Chem. Sci. 11, 12777–12788 (2020).
