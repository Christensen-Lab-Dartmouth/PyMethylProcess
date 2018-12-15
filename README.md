# PyMethylProcess

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

Division of Labor:
