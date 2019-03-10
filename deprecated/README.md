# PyMethylProcess
# Make sure to have GCC-8 installed env CC=/usr/local/Cellar/gcc/8.2.0/bin/gcc-8
# env CC=/usr/local/Cellar/gcc/8.2.0/bin/gcc-8 pip install dist/pymethylprocess-0.1.tar.gz
1. Download TCGA/GEO Data
2. Preprocess using meffil+minfi+enmix in parallel
    * Preprocessing objects
3. Imputation and MAD feature selection
4. Output as input into MethylNet
5. Wrapped into docker, CWL, and Toil
    * pypi+anaconda

For developers:
1. To give you an idea of the install requirements just for preprocessing, I have a docker file containing the requirements.
2. Please leave existing command line commands in place for now. We can modify this.

Tasks:
1. Finish creation and cleaning of downloader and preprocessing classes.
2. Add impyute, DAPL and DeepCpG to imputation methods.
3. Finish testing CWL pipeline and separate it from methylnet.
4. MethylArray class must be able to be imported into MethylNet and should be pickleable
5. Thorough documentation.
6. Submit to journal as application note that MethylNet will reference.
7. Should help docs be included with methylnet or separate?
8. Host help docs on github.io

"Wrapping the methylation pipeline in Python is a great contribution to the field in itself. We should consider cleaning that into its own repo and possibly submitting an application note to Bioinformatics or at a minimum bioRxiv to let the community use that as a separate resource.

And if we do that we should create a PyPi package for easy use. @jlevy44 we should chat about a division of labor."
- @AlexanderTitus

Division of Labor:

Docker:
module load singularity

----------
module load python/3-Anaconda
conda create -n methylnet_runner python=3.6
source activate methylnet_runner
conda install -c bioconda udocker
udocker login --username=joshualevy44 --password=xxx # --registry=cloud.canister.io:5000
udocker pull joshualevy44/methylnet:preprocess #cloud.canister.io:5000/joshualevy44/methylnet:preprocess
udocker create --name=methylnet joshualevy44/methylnet:preprocess #cloud.canister.io:5000/joshualevy44/methylnet:preprocess
udocker setup --execmode=F3  methylnet
# for GPUs
udocker setup --nvidia --execmode=F3  methylnet
udocker run methylnet /scripts/preprocess.py

Base docker image:
cloud.canister.io:5000/joshualevy44/methylnet:base
Preprocess docker image
cloud.canister.io:5000/joshualevy44/methylnet:preprocess
