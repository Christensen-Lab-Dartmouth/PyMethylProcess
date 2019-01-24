import click
import os, subprocess, glob
from os.path import join
import time
from MethylationDataTypes import MethylationArray

CONTEXT_SETTINGS = dict(help_option_names=['-h','--help'], max_content_width=90)

@click.group(context_settings= CONTEXT_SETTINGS)
@click.version_option(version='0.1')
def util():
    pass


@util.command()
@click.option('-i', '--input_pkl', default='./final_preprocessed/methyl_array.pkl', help='Input database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-o', '--output_dir', default='./train_val_test_sets/', help='Output directory for training, testing, and validation sets.', type=click.Path(exists=False), show_default=True)
@click.option('-tp', '--train_percent', default=0.8, help='Percent data training on.', show_default=True)
@click.option('-vp', '--val_percent', default=0.1, help='Percent of training data that comprises validation set.', show_default=True)
@click.option('-cat', '--categorical', is_flag=True, help='Multi-class prediction.', show_default=True)
@click.option('-do', '--disease_only', is_flag=True, help='Only look at disease, or text before subtype_delimiter.')
@click.option('-k', '--key', default='disease', help='Key to split on.', type=click.Path(exists=False), show_default=True)
@click.option('-sd', '--subtype_delimiter', default=',', help='Delimiter for disease extraction.', type=click.Path(exists=False), show_default=True)
def train_test_val_split(input_pkl,output_dir,train_percent,val_percent, categorical, disease_only, key, subtype_delimiter):
    os.makedirs(output_dir,exist_ok=True)
    methyl_array = MethylationArray.from_pickle(input_pkl)
    train_arr, test_arr, val_arr = methyl_array.split_train_test(train_percent, categorical, disease_only, key, subtype_delimiter, val_percent)
    train_arr.write_pickle(join(output_dir,'train_methyl_array.pkl'))
    test_arr.write_pickle(join(output_dir,'test_methyl_array.pkl'))
    val_arr.write_pickle(join(output_dir,'val_methyl_array.pkl'))

@util.command()
@click.option('-i', '--input_pkl', default='./final_preprocessed/methyl_array.pkl', help='Input database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-k', '--key', default='disease', help='Key to split on.', type=click.Path(exists=False), show_default=True)
def counts(input_pkl,key):
    if input_pkl.endswith('.pkl'):
        MethylationArray.from_pickle(input_pkl).categorical_breakdown(key)
    else:
        for input_pkl in glob.glob(join(input_pkl,'*.pkl')):
            print(input_pkl)
            MethylationArray.from_pickle(input_pkl).categorical_breakdown(key)

@util.command()
@click.option('-i', '--input_pkl', default='./final_preprocessed/methyl_array.pkl', help='Input database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-k', '--key', default='disease', help='Key to split on.', type=click.Path(exists=False), show_default=True)
@click.option('-d', '--disease_only', is_flag=True, help='Only look at disease, or text before subtype_delimiter.')
@click.option('-sd', '--subtype_delimiter', default=',', help='Delimiter for disease extraction.', type=click.Path(exists=False), show_default=True)
@click.option('-o', '--output_pkl', default='./fixed_preprocessed/methyl_array.pkl', help='Input database for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
def fix_key(input_pkl,key,disease_only,subtype_delimiter,output_pkl):
    os.makedirs(output_pkl[:output_pkl.rfind('/')],exist_ok=True)
    methyl_array=MethylationArray.from_pickle(input_pkl)
    methyl_array.remove_whitespace(key)
    if disease_only:
        methyl_array.split_key(key, subtype_delimiter)
    methyl_array.write_pickle(output_pkl)

if __name__ == '__main__':
    util()
