import os, subprocess
import click
import glob
import numpy as np, pandas as pd
#import impyute
from functools import reduce
import pickle
from pymethylprocess.PreProcessDataTypes import *
from pymethylprocess.MethylationDataTypes import *

CONTEXT_SETTINGS = dict(help_option_names=['-h','--help'], max_content_width=90)

@click.group(context_settings= CONTEXT_SETTINGS)
@click.version_option(version='0.1')
def preprocess():
    pass

#### COMMANDS ####

## Download ##

@preprocess.command()
@click.option('-o', '--output_dir', default='./tcga_idats/', help='Output directory for exported idats.', type=click.Path(exists=False), show_default=True)
def download_tcga(output_dir):
    """Download all tcga 450k data."""
    os.makedirs(output_dir, exist_ok=True)
    downloader = TCGADownloader()
    downloader.download_tcga(output_dir)

@preprocess.command()
@click.option('-o', '--output_dir', default='./tcga_idats/', help='Output directory for exported idats.', type=click.Path(exists=False), show_default=True)
def download_clinical(output_dir):
    """Download all TCGA 450k clinical info."""
    os.makedirs(output_dir, exist_ok=True)
    downloader = TCGADownloader()
    downloader.download_clinical(output_dir)

@preprocess.command()
@click.option('-g', '--geo_query', default='', help='GEO study to query.', type=click.Path(exists=False), show_default=True)
@click.option('-o', '--output_dir', default='./geo_idats/', help='Output directory for exported idats.', type=click.Path(exists=False), show_default=True)
def download_geo(geo_query,output_dir):
    """Download geo methylation study idats and clinical info."""
    os.makedirs(output_dir, exist_ok=True)
    downloader = TCGADownloader()
    downloader.download_geo(geo_query,output_dir)

## prepare

@preprocess.command()
@click.option('-is', '--input_sample_sheet', default='./tcga_idats/clinical_info.csv', help='Clinical information downloaded from tcga/geo/custom.', type=click.Path(exists=False), show_default=True)
@click.option('-s', '--source_type', default='tcga', help='Source type of data.', type=click.Choice(['tcga','geo','custom']), show_default=True)
@click.option('-i', '--idat_dir', default='./tcga_idats/', help='Idat directory.', type=click.Path(exists=False), show_default=True)
@click.option('-os', '--output_sample_sheet', default='./tcga_idats/minfiSheet.csv', help='CSV for minfi input.', type=click.Path(exists=False), show_default=True)
@click.option('-m', '--mapping_file', default='./idat_filename_case.txt', help='Mapping file from uuid to TCGA barcode. Downloaded using download_tcga.', type=click.Path(exists=False), show_default=True)
@click.option('-l', '--header_line', default=0, help='Line to begin reading csv/xlsx.', show_default=True)
@click.option('-d', '--disease_class_column', default="methylation class:ch1", help='Disease classification column, for custom and geo datasets.', type=click.Path(exists=False), show_default=True)
@click.option('-b', '--basename_col', default="Sentrix ID (.idat)", help='Basename classification column, for custom datasets.', type=click.Path(exists=False), show_default=True)
@click.option('-c', '--include_columns_file', default="", help='Custom columns file containing columns to keep, separated by \\n. Add a tab for each line if you wish to rename columns: original_name \\t new_column_name', type=click.Path(exists=False), show_default=True)
def create_sample_sheet(input_sample_sheet, source_type, idat_dir, output_sample_sheet, mapping_file, header_line, disease_class_column, basename_col, include_columns_file):
    """Create sample sheet for input to minfi, meffil, or enmix."""
    os.makedirs(output_sample_sheet[:output_sample_sheet.rfind('/')], exist_ok=True)
    pheno_sheet = PreProcessPhenoData(input_sample_sheet, idat_dir, header_line= (0 if source_type != 'custom' else header_line))
    if include_columns_file:
        include_columns=np.loadtxt(include_columns_file,dtype=str,delimiter='\t')
        if '\t' in open(include_columns_file).read() and len(include_columns.shape)<2:
            include_columns=dict(include_columns[np.newaxis,:].tolist())
        elif len(include_columns.shape)<2:
            include_columns=dict(zip(include_columns,include_columns))
        else:
            include_columns=dict(include_columns.tolist())
    else:
        include_columns={}
    if source_type == 'tcga':
        pheno_sheet.format_tcga(mapping_file)
    elif source_type == 'geo':
        pheno_sheet.format_geo(disease_class_column, include_columns)
    else:
        pheno_sheet.format_custom(basename_col, disease_class_column, include_columns)
    pheno_sheet.export(output_sample_sheet)
    print("Please remove {} from {}, if it exists in that directory.".format(input_sample_sheet, idat_dir))

@preprocess.command()
@click.option('-is', '--input_sample_sheet', default='./tcga_idats/minfiSheet.csv', help='CSV for minfi input.', type=click.Path(exists=False), show_default=True)
@click.option('-os', '--output_sample_sheet', default='./tcga_idats/minfiSheet.csv', help='CSV for minfi input.', type=click.Path(exists=False), show_default=True)
def meffil_encode(input_sample_sheet,output_sample_sheet):
    """Reformat file for meffil input."""
    from collections import defaultdict
    pheno=pd.read_csv(input_sample_sheet)
    sex_dict=defaultdict(lambda:'NA')
    sex_dict.update({'m':'M','f':'F','M':'M','F':'F','male':'M','female':'F','nan':'NA',np.nan:'NA'})
    k='Sex' if 'Sex' in list(pheno) else 'sex'
    if k in list(pheno):
        pheno.loc[:,k] = pheno[k].map(lambda x: sex_dict[str(x).lower()])
    pheno = pheno[[col for col in list(pheno) if not col.startswith('Unnamed:')]].rename(columns={'sex':'Sex'})
    if 'Sex' in list(pheno):
        d=defaultdict(lambda:'NA')
        d.update({'M':'M','F':'F'})
        pheno.loc[:,'Sex'] = pheno['Sex'].map(d)
        if (pheno['Sex']==pheno['Sex'].mode().values[0][0]).all():
            pheno=pheno.rename(columns={'Sex':'gender'})
    pheno.to_csv(output_sample_sheet)

@preprocess.command()
@click.option('-s1', '--sample_sheet1', default='./tcga_idats/clinical_info1.csv', help='Clinical information downloaded from tcga/geo/custom, formatted using create_sample_sheet.', type=click.Path(exists=False), show_default=True)
@click.option('-s2', '--sample_sheet2', default='./tcga_idats/clinical_info2.csv', help='Clinical information downloaded from tcga/geo/custom, formatted using create_sample_sheet.', type=click.Path(exists=False), show_default=True)
@click.option('-os', '--output_sample_sheet', default='./tcga_idats/minfiSheet.csv', help='CSV for minfi input.', type=click.Path(exists=False), show_default=True)
@click.option('-d', '--second_sheet_disease', is_flag=True, help='Use second sheet\'s disease column.')
@click.option('-nd', '--no_disease_merge', is_flag=True, help='Don\'t merge disease columns.')
def merge_sample_sheets(sample_sheet1, sample_sheet2, output_sample_sheet, second_sheet_disease,no_disease_merge):
    """Merge two sample files for more fields for minfi+ input."""
    s1 = PreProcessPhenoData(sample_sheet1, idat_dir='', header_line=0)
    s2 = PreProcessPhenoData(sample_sheet2, idat_dir='', header_line=0)
    s1.merge(s2,second_sheet_disease, no_disease_merge=no_disease_merge)
    s1.export(output_sample_sheet)

@preprocess.command()
@click.option('-s1', '--sample_sheet1', default='./tcga_idats/clinical_info1.csv', help='Clinical information downloaded from tcga/geo/custom, formatted using create_sample_sheet.', type=click.Path(exists=False), show_default=True)
@click.option('-s2', '--sample_sheet2', default='./tcga_idats/clinical_info2.csv', help='Clinical information downloaded from tcga/geo/custom, formatted using create_sample_sheet.', type=click.Path(exists=False), show_default=True)
@click.option('-os', '--output_sample_sheet', default='./tcga_idats/minfiSheet.csv', help='CSV for minfi input.', type=click.Path(exists=False), show_default=True)
def concat_sample_sheets(sample_sheet1, sample_sheet2, output_sample_sheet):
    """Concat two sample files for more fields for minfi+ input, adds more samples."""
    # FIXME add ability to concat more sample sheets; dump to sql!!!
    s1 = PreProcessPhenoData(sample_sheet1, idat_dir='', header_line=0)
    s2 = PreProcessPhenoData(sample_sheet2, idat_dir='', header_line=0)
    s1.concat(s2)
    s1.export(output_sample_sheet)

@preprocess.command()
@click.option('-is', '--formatted_sample_sheet', default='./tcga_idats/minfiSheet.csv', help='Clinical information downloaded from tcga/geo/custom, formatted using create_sample_sheet.', type=click.Path(exists=False), show_default=True)
@click.option('-k', '--key', multiple=True, default=['disease'], help='Column of csv to print statistics for.', type=click.Path(exists=False), show_default=True)
@click.option('-d', '--disease_only', is_flag=True, help='Only look at disease, or text before subtype_delimiter.')
@click.option('-sd', '--subtype_delimiter', default=',', help='Delimiter for disease extraction.', type=click.Path(exists=False), show_default=True)
def get_categorical_distribution(formatted_sample_sheet,key,disease_only=False,subtype_delimiter=','):
    """Get categorical distribution of columns of sample sheet."""
    if len(key) == 1:
        key=key[0]
    print('\n'.join('{}:{}'.format(k,v) for k,v in PreProcessPhenoData(formatted_sample_sheet, idat_dir='', header_line=0).get_categorical_distribution(key,disease_only,subtype_delimiter).items()))

@preprocess.command()
@click.option('-is', '--formatted_sample_sheet', default='./tcga_idats/clinical_info.csv', help='Clinical information downloaded from tcga/geo/custom, formatted using create_sample_sheet.', type=click.Path(exists=False), show_default=True)
@click.option('-e', '--exclude_disease_list', default='', help='List of conditions to exclude, from disease column, comma delimited.', type=click.Path(exists=False), show_default=True)
@click.option('-os', '--output_sheet_name', default='./tcga_idats/minfiSheet.csv', help='CSV for minfi input.', type=click.Path(exists=False), show_default=True)
@click.option('-l', '--low_count', default=0, help='Remove diseases if they are below a certain count, default this is not used.', type=click.Path(exists=False), show_default=True)
@click.option('-d', '--disease_only', is_flag=True, help='Only look at disease, or text before subtype_delimiter.')
@click.option('-sd', '--subtype_delimiter', default=',', help='Delimiter for disease extraction.', type=click.Path(exists=False), show_default=True)
def remove_diseases(formatted_sample_sheet, exclude_disease_list, output_sheet_name, low_count, disease_only=False,subtype_delimiter=','):
    """Exclude diseases from study by count number or exclusion list."""
    exclude_disease_list = exclude_disease_list.split(',')
    pData = PreProcessPhenoData(formatted_sample_sheet, idat_dir='', header_line=0)
    pData.remove_diseases(exclude_disease_list,low_count, disease_only,subtype_delimiter)
    pData.export(output_sheet_name)
    print("Please remove {} from idat directory, if it exists in that directory.".format(formatted_sample_sheet))

## preprocess ##

@preprocess.command()
@click.option('-i', '--idat_csv', default='./tcga_idats/minfiSheet.csv', help='Idat csv for one sample sheet, alternatively can be your phenotype sample sheet.', type=click.Path(exists=False), show_default=True)
@click.option('-d', '--disease_only', is_flag=True, help='Only look at disease, or text before subtype_delimiter.')
@click.option('-sd', '--subtype_delimiter', default=',', help='Delimiter for disease extraction.', type=click.Path(exists=False), show_default=True)
@click.option('-o', '--subtype_output_dir', default='./preprocess_outputs/', help='Output subtypes pheno csv.', type=click.Path(exists=False), show_default=True)
def split_preprocess_input_by_subtype(idat_csv,disease_only,subtype_delimiter, subtype_output_dir):
    """Split preprocess input samplesheet by disease subtype."""
    from collections import defaultdict
    subtype_delimiter=subtype_delimiter.replace('"','').replace("'","")
    os.makedirs(subtype_output_dir,exist_ok=True)
    pData=PreProcessPhenoData(idat_csv,'')
    idat_csv_basename = idat_csv.split('/')[-1]
    group_by_key = (pData.split_key('disease',subtype_delimiter) if disease_only else 'disease')
    pData_grouped = pData.pheno_sheet.groupby(group_by_key)
    default_qc_df=[]
    for name, group in pData_grouped:
        name=name.replace(' ','')
        new_sheet = idat_csv_basename.replace('.csv','_{}.csv'.format(name))
        new_out_dir = '{}/{}/'.format(subtype_output_dir,name)
        os.makedirs(new_out_dir, exist_ok=True)
        print(new_out_dir)
        if 'Sex' in list(group):
            d=defaultdict(lambda:'NA')
            d.update({'M':'M','F':'F'})
            group.loc[:,'Sex'] = group['Sex'].map(d)
            if (group['Sex']==group['Sex'].mode().values[0][0]).all():
                group=group.rename(columns={'Sex':'gender'})
        group.to_csv('{}/{}'.format(new_out_dir,new_sheet))
        default_qc_df.append([name,-1,0.05,0.05,0.05,0.05,5,-2])
    pd.DataFrame(default_qc_df,columns=['subtype','n_pcs','p_beadnum_samples','p_detection_samples','p_beadnum_cpgs','p_detection_cpgs','sex_sd','sex_cutoff']).to_csv(os.path.join(subtype_output_dir,'pc_qc_parameters.csv'))



@preprocess.command()
@click.option('-n', '--n_cores', default=6, help='Number cores to use for preprocessing.', show_default=True)
@click.option('-i', '--subtype_output_dir', default='./preprocess_outputs/', help='Output subtypes pheno csv.', type=click.Path(exists=False), show_default=True)
@click.option('-m', '--meffil', is_flag=True, help='Preprocess using meffil.')
@click.option('-t', '--torque', is_flag=True, help='Job submission torque.')
@click.option('-r', '--run', is_flag=True, help='Actually run local job or just print out command.')
@click.option('-s', '--series', is_flag=True, help='Run commands in series.')
@click.option('-p', '--pc_qc_parameters_csv', default='./preprocess_outputs/pc_qc_parameters.csv', show_default=True, help='For meffil, qc parameters and pcs for final qc and functional normalization.')
@click.option('-u', '--use_cache', is_flag=True, help='If this is selected, loads qc results rather than running qc again. Only works for meffil selection.')
@click.option('-qc', '--qc_only', is_flag=True, help='Only perform QC for meffil pipeline, caches results into rds file for loading again, only works if use_cache is false.')
@click.option('-c', '--chunk_size', default=-1, help='If not series, chunk up and run these number of commands at once.. -1 means all commands at once.')
def batch_deploy_preprocess(n_cores,subtype_output_dir,meffil,torque,run,series, pc_qc_parameters_csv, use_cache, qc_only, chunk_size):
    """Deploy multiple preprocessing jobs in series or parallel."""
    pheno_csvs = glob.glob(os.path.join(subtype_output_dir,'*','*.csv'))
    opts = {'-n':n_cores}
    try:
        pc_qc_parameters = pd.read_csv(pc_qc_parameters).drop_duplicates().set_index('subtype')
    except:
        pc_qc_parameters = pd.DataFrame([[name,-1,0.05,0.05,0.05,0.05,5,-2] for name in np.vectorize(lambda x: x.split('/')[-2])(pheno_csvs)],
                        columns=['subtype','n_pcs','p_beadnum_samples','p_detection_samples','p_beadnum_cpgs','p_detection_cpgs','sex_sd','sex_cutoff']).drop_duplicates().set_index('subtype')
    if meffil:
        opts['-m']=''
    if use_cache:
        opts['-u']=''
    if qc_only:
        opts['-qc']=''
    commands=[]
    for pheno_csv in pheno_csvs:
        pheno_path = os.path.abspath(pheno_csv)
        subtype=pheno_path.split('/')[-2]
        opts['-pc'] = int(pc_qc_parameters.loc[subtype,'n_pcs'])
        opts['-bns'] = pc_qc_parameters.loc[subtype,'p_beadnum_samples']
        opts['-pds'] = pc_qc_parameters.loc[subtype,'p_detection_samples']
        opts['-bnc'] = pc_qc_parameters.loc[subtype,'p_beadnum_cpgs']
        opts['-pdc'] = pc_qc_parameters.loc[subtype,'p_detection_cpgs']
        opts['-sc'] = pc_qc_parameters.loc[subtype,'sex_cutoff']
        opts['-sd'] = pc_qc_parameters.loc[subtype,'sex_sd']
        opts['-i']=pheno_path[:pheno_path.rfind('/')+1]
        opts['-o']=pheno_path[:pheno_path.rfind('/')+1]+'methyl_array.pkl'
        command='pymethyl-preprocess preprocess_pipeline {}'.format(' '.join('{} {}'.format(k,v) for k,v in opts.items()))
        commands.append(command)
    if not torque:
        if not series and chunk_size != -1:
            #commands = np.array_split(commands,len(commands)//chunk_size)
            print(commands)
            with open('commands.txt','w') as f:
                f.write('\n'.join(commands))
            subprocess.call('cat commands.txt | xargs -L 1 -I CMD -P {} bash -c CMD'.format(chunk_size),shell=True) # https://www.gnu.org/software/parallel/sem.html
            """for command_list in commands:
                subprocess.call('run_parallel {}'.format(' '.join(['"{}"'.format(command) for command in command_list])),shell=True)"""
        else:
            for command in commands:
                if not series:
                    command="nohup {} &".format(command)
                if not run:
                    click.echo(command)
                else:
                    subprocess.call(command,shell=True)
    else:
        run_command = lambda command: subprocess.call('module load cuda && module load python/3-Anaconda && source activate py36 && {}'.format(command),shell=True)
        from pyina.schedulers import Torque
        from pyina.launchers import Mpi
        config = {'nodes':'1:ppn=6', 'queue':'default', 'timelimit':'01:00:00'}
        torque = Torque(**config)
        pool = Mpi(scheduler=torque)
        pool.map(run_command, commands)


@preprocess.command()
@click.option('-i', '--idat_dir', default='./tcga_idats/', help='Idat dir for one sample sheet, alternatively can be your phenotype sample sheet.', type=click.Path(exists=False), show_default=True)
@click.option('-n', '--n_cores', default=6, help='Number cores to use for preprocessing.', show_default=True)
@click.option('-o', '--output_pkl', default='./preprocess_outputs/methyl_array.pkl', help='Output database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-m', '--meffil', is_flag=True, help='Preprocess using meffil.')
@click.option('-pc', '--n_pcs', default=-1, show_default=True, help='For meffil, number of principal components for functional normalization. If set to -1, then PCs are selected using elbow method.')
@click.option('-p', '--pipeline', default='enmix', show_default=True, help='If not meffil, preprocess using minfi or enmix.', type=click.Choice(['minfi','enmix']))
@click.option('-noob', '--noob_norm', is_flag=True, help='Run noob normalization of minfi selected.')
@click.option('-u', '--use_cache', is_flag=True, help='If this is selected, loads qc results rather than running qc again and update with new qc parameters. Only works for meffil selection. Minfi and enmix just loads RG Set.')
@click.option('-qc', '--qc_only', is_flag=True, help='Only perform QC for meffil pipeline, caches results into rds file for loading again, only works if use_cache is false. Minfi and enmix just saves the RGSet before preprocessing.')
@click.option('-bns', '--p_beadnum_samples', default=0.05, show_default=True, help='From meffil documentation, "fraction of probes that failed the threshold of 3 beads".')
@click.option('-pds', '--p_detection_samples', default=0.05, show_default=True, help='From meffil documentation, "fraction of probes that failed a detection.pvalue threshold of 0.01".')
@click.option('-bnc', '--p_beadnum_cpgs', default=0.05, show_default=True, help='From meffil documentation, "fraction of samples that failed the threshold of 3 beads".')
@click.option('-pdc', '--p_detection_cpgs', default=0.05, show_default=True, help='From meffil documentation, "fraction of samples that failed a detection.pvalue threshold of 0.01".')
@click.option('-sc', '--sex_cutoff', default=-2, show_default=True, help='From meffil documentation, "difference of total median intensity for Y chromosome probes and X chromosome probes".')
@click.option('-sd', '--sex_sd', default=5, show_default=True, help='From meffil documentation, "sex detection outliers if outside this range".')
def preprocess_pipeline(idat_dir, n_cores, output_pkl, meffil, n_pcs, pipeline, noob_norm, use_cache, qc_only, p_beadnum_samples, p_detection_samples, p_beadnum_cpgs, p_detection_cpgs, sex_cutoff, sex_sd):
    """Perform preprocessing of idats using enmix or meffil."""
    output_dir = output_pkl[:output_pkl.rfind('/')]
    os.makedirs(output_dir,exist_ok=True)
    preprocesser = PreProcessIDAT(idat_dir)
    if meffil:
        qc_parameters={'p.beadnum.samples':p_beadnum_samples,'p.detection.samples':p_detection_samples,'p.detection.cpgs':p_detection_cpgs,'p.beadnum.cpgs':p_beadnum_cpgs,'sex.cutoff':sex_cutoff, 'sex.outlier.sd':sex_sd}
        preprocesser.preprocessMeffil(n_cores=n_cores,n_pcs=n_pcs,qc_report_fname=os.path.join(output_dir,'qc.report.html'), normalization_report_fname=os.path.join(output_dir,'norm.report.html'), pc_plot_fname=os.path.join(output_dir,'pc.plot.pdf'), useCache=use_cache, qc_only=qc_only, qc_parameters=qc_parameters)
    else:
        preprocesser.preprocess_enmix_pipeline(n_cores=n_cores, pipeline=pipeline, noob=noob_norm, use_cache=use_cache, qc_only=qc_only)
        try:
            preprocesser.plot_qc_metrics(output_dir)
        except:
            pass
    preprocesser.output_pheno_beta(meffil=meffil)
    preprocesser.to_methyl_array('').write_pickle(output_pkl)

@preprocess.command()
@click.option('-i', '--input_pkls', default=['./preprocess_outputs/methyl_array.pkl'], multiple=True, help='Input pickles for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-d', '--optional_input_pkl_dir', default='', help='Auto grab input pkls.', type=click.Path(exists=False), show_default=True)
@click.option('-o', '--output_pkl', default='./combined_outputs/methyl_array.pkl', help='Output database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-e', '--exclude', default=[], multiple=True, help='If -d selected, these diseases will be excluded from study.', show_default=True)
def combine_methylation_arrays(input_pkls, optional_input_pkl_dir, output_pkl, exclude):
    """If split MethylationArrays by subtype for either preprocessing or imputation, can use to recombine data for downstream step."""
    os.makedirs(output_pkl[:output_pkl.rfind('/')],exist_ok=True)
    if optional_input_pkl_dir:
        input_pkls=glob.glob(os.path.join(optional_input_pkl_dir,'*','methyl_array.pkl'))
        if exclude:
            input_pkls=(np.array(input_pkls)[~np.isin(np.vectorize(lambda x: x.split('/')[-2])(input_pkls),np.array(exclude))]).tolist()
    if len(input_pkls) > 0:
        base_methyl_array=MethylationArray(*extract_pheno_beta_df_from_pickle_dict(pickle.load(open(input_pkls[0],'rb')), ''))
        methyl_arrays_generator = (MethylationArray(*extract_pheno_beta_df_from_pickle_dict(pickle.load(open(input_pkl,'rb')), '')) for input_pkl in input_pkls[1:])
        list_methyl_arrays = MethylationArrays([base_methyl_array])
        combined_methyl_array = list_methyl_arrays.combine(methyl_arrays_generator)
    else:
        combined_methyl_array=MethylationArray(*extract_pheno_beta_df_from_pickle_dict(pickle.load(open(input_pkls[0],'rb')), ''))
    combined_methyl_array.write_pickle(output_pkl)

@preprocess.command()
@click.option('-i', '--input_pkl', default='./combined_outputs/methyl_array.pkl', help='Input database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-ss', '--split_by_subtype', is_flag=True, help='Imputes CpGs by subtype before combining again.')
@click.option('-m', '--method', default='KNN', help='Method of imputation.', type=click.Choice(['KNN', 'Mean', 'Zero', 'MICE', 'BiScaler', 'Soft', 'random', 'DeepCpG', 'DAPL']), show_default=True)
@click.option('-s', '--solver', default='fancyimpute', help='Imputation library.', type=click.Choice(['fancyimpute', 'impyute', 'sklearn']), show_default=True)
@click.option('-k', '--n_neighbors', default=5, help='Number neighbors for imputation if using KNN.', show_default=True)
@click.option('-r', '--orientation', default='Samples', help='Impute CpGs or samples.', type=click.Choice(['Samples','CpGs']), show_default=True)
@click.option('-o', '--output_pkl', default='./imputed_outputs/methyl_array.pkl', help='Output database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-n', '--n_top_cpgs', default=0, help='Number cpgs to include with highest variance across population. Greater than 0 allows for mad filtering during imputation to skip mad step.', show_default=True)
@click.option('-f', '--feature_selection_method', default='mad', type=click.Choice(['mad','spectral']))
@click.option('-mm', '--metric', default='correlation', type=click.Choice(['euclidean','cosine','correlation']))
@click.option('-nfs', '--n_neighbors_fs', default=0, help='Number neighbors for feature selection, default enacts rbf kernel.', show_default=True)
@click.option('-d', '--disease_only', is_flag=True, help='Only look at disease, or text before subtype_delimiter.')
@click.option('-sd', '--subtype_delimiter', default=',', help='Delimiter for disease extraction.', type=click.Path(exists=False), show_default=True)
@click.option('-st', '--sample_threshold', default=-1., help='Value between 0 and 1 for NaN removal. If samples has sample_threshold proportion of cpgs missing, then remove sample. Set to -1 to not remove samples.', show_default=True)
@click.option('-ct', '--cpg_threshold', default=-1., help='Value between 0 and 1 for NaN removal. If cpgs has cpg_threshold proportion of samples missing, then remove cpg. Set to -1 to not remove samples.', show_default=True)
def imputation_pipeline(input_pkl,split_by_subtype=True,method='knn', solver='fancyimpute', n_neighbors=5, orientation='rows', output_pkl='', n_top_cpgs=0, feature_selection_method='mad', metric='correlation', n_neighbors_fs=10, disease_only=False, subtype_delimiter=',', sample_threshold=-1, cpg_threshold=-1): # wrap a class around this
    """Imputation of subtype or no subtype using various imputation methods."""
    orientation_dict = {'CpGs':'columns','Samples':'rows'}
    orientation = orientation_dict[orientation]
    #print("Selecting orientation for imputation not implemented yet.")
    os.makedirs(output_pkl[:output_pkl.rfind('/')],exist_ok=True)
    if method in ['DeepCpG', 'DAPL', 'EM']:
        print('Method {} coming soon...'.format(method))
    elif solver in ['impyute']:
        print('Impyute coming soon...')
    else:
        imputer = ImputerObject(solver, method, opts=dict(k=n_neighbors, orientation=orientation, min_value=0., max_value=1.)).return_imputer()
    input_dict = pickle.load(open(input_pkl,'rb'))

    if cpg_threshold == -1.:
        cpg_threshold = None

    if sample_threshold == -1.:
        sample_threshold = None

    if split_by_subtype:

        def impute_arrays(methyl_arrays):
            for methyl_array in methyl_arrays:
                methyl_array.remove_missingness(cpg_threshold=cpg_threshold, sample_threshold=sample_threshold)
                methyl_array.impute(imputer)
                yield methyl_array

        methyl_array = MethylationArray(*extract_pheno_beta_df_from_pickle_dict(input_dict))
        methyl_arrays = impute_arrays(methyl_array.split_by_subtype(disease_only, subtype_delimiter))
        methyl_array=MethylationArrays([next(methyl_arrays)]).combine(methyl_arrays)

    else:
        methyl_array = MethylationArray(*extract_pheno_beta_df_from_pickle_dict(input_dict))
        methyl_array.remove_missingness(cpg_threshold=cpg_threshold, sample_threshold=sample_threshold)
        methyl_array.impute(imputer)

    if n_top_cpgs:
        methyl_array.feature_select(n_top_cpgs,feature_selection_method, metric, nn=n_neighbors_fs)

    methyl_array.write_pickle(output_pkl)

@preprocess.command()
@click.option('-i', '--input_pkl', default='./imputed_outputs/methyl_array.pkl', help='Input database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-o', '--output_pkl', default='./final_preprocessed/methyl_array.pkl', help='Output database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-n', '--n_top_cpgs', default=300000, help='Number cpgs to include with highest variance across population.', show_default=True)
@click.option('-f', '--feature_selection_method', default='mad', type=click.Choice(['mad','spectral']))
@click.option('-mm', '--metric', default='correlation', type=click.Choice(['euclidean','cosine','correlation']))
@click.option('-nn', '--n_neighbors', default=0, help='Number neighbors for feature selection, default enacts rbf kernel.', show_default=True)
@click.option('-m', '--mad_top_cpgs', default=0, help='Number cpgs to apply mad filtering first before more sophisticated feature selection. If 0 or primary feature selection is mad, no mad pre-filtering.', show_default=True)
def feature_select(input_pkl,output_pkl,n_top_cpgs=300000, feature_selection_method='mad', metric='correlation', n_neighbors=10, mad_top_cpgs=0):
    """Filter CpGs by taking x top CpGs with highest mean absolute deviation scores or via spectral feature selection."""
    os.makedirs(output_pkl[:output_pkl.rfind('/')],exist_ok=True)
    input_dict=pickle.load(open(input_pkl,'rb'))
    methyl_array = MethylationArray(*extract_pheno_beta_df_from_pickle_dict(input_dict))

    if mad_top_cpgs and feature_selection_method != 'mad':
        methyl_array.feature_select(mad_top_cpgs,'mad')

    methyl_array.feature_select(n_top_cpgs,feature_selection_method, metric, nn=n_neighbors)

    methyl_array.write_pickle(output_pkl)

    return methyl_array

### QC

@preprocess.command()
@click.option('-i', '--input_pkl', default='./preprocess_outputs/methyl_array.pkl', help='Input database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-o', '--output_dir', default='./na_report/', help='Output database for na report.', type=click.Path(exists=False), show_default=True)
@click.option('-r', '--head_directory', is_flag=True, help='-i option becomes directory, and searches there for multiple input pickles.')
def na_report(input_pkl, output_dir, head_directory):
    """Print proportion of missing values throughout dataset."""
    import matplotlib.pyplot as plt
    if head_directory:
        input_pickles=glob.glob(os.path.join(input_pkl,'*/*.pkl'))
        output_dirs=np.vectorize(lambda x: x.replace('.pkl','.na_report'))(input_pickles).tolist()
    else:
        input_pickles=[input_pkl]
        output_dirs=[output_dir]
    for input_pkl,output_dir in zip(input_pickles,output_dirs):
        os.makedirs(output_dir,exist_ok=True)
        df=pickle.load(open(input_pkl,'rb'))['beta']
        na_frame = pd.isna(df).astype(int)
        na_frame = na_frame.iloc[np.argsort(na_frame.sum(axis=1)),:]
        print('{} NA Rate is on average: {}%'.format(input_pkl,sum(na_frame.values.flatten())/float(df.shape[0]*df.shape[1])*100.))
        plt.figure()
        pd.DataFrame(na_frame.sum(axis=1)).apply(lambda x: x/float(df.shape[1])).hist()
        plt.savefig(os.path.join(output_dir,'sample_missingness.png'))
        plt.figure()
        pd.DataFrame(na_frame.sum(axis=0)).apply(lambda x: x/float(df.shape[0])).hist()
        plt.savefig(os.path.join(output_dir,'cpg_missingness.png'))
        xy=np.transpose(np.nonzero(na_frame.values))
        if 0:
            plt.figure()
            plt.scatter(xy[:,1],xy[:,0])
            plt.xlim([0,df.shape[0]])
            plt.ylim([0,df.shape[1]])
            plt.xlabel('CpGs')
            plt.ylabel('Samples')
            plt.savefig(os.path.join(output_dir,'array_missingness.png'))



#################

if __name__ == '__main__':
    preprocess()
