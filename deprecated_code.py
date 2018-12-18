

@preprocess.command() # update
@click.option('-i', '--idat_dir', default='./tcga_idats/', help='Idat directory if one sample sheet, alternatively can be your phenotype sample sheet.', type=click.Path(exists=False), show_default=True)
@click.option('-g', '--geo_query', default='', help='GEO study to query, do not use if already created geo sample sheet.', type=click.Path(exists=False), show_default=True)
@click.option('-o', '--output_dir', default='./preprocess_outputs/', help='Output directory for beta and phenotype data.', type=click.Path(exists=False), show_default=True)
@click.option('-ss', '--split_by_subtype', is_flag=True, help='If using formatted sample sheet csv, split by subtype and perform preprocessing. Will need to combine later.')
def plot_qc(idat_dir, geo_query, output_dir, split_by_subtype):
    """Plot QC metrics using raw preprocessing via minfi and enmix."""
    os.makedirs(output_dir, exist_ok=True)
    if idat_dir.endswith('.csv') and split_by_subtype:
        pheno=pd.read_csv(idat_dir)
        for name, group in pheno.groupby('disease'):
            new_sheet = idat_dir.replace('.csv','_{}.csv'.format(name)).split('/')[-1]
            new_out_dir = '{}/{}/'.format(output_dir,name)
            os.makedirs(new_out_dir, exist_ok=True)
            group.to_csv('{}/{}'.format(new_out_dir,new_sheet))
            preprocesser = PreProcessIDAT(new_out_dir)
            preprocesser.load_idats(geo_query='')
            preprocesser.plot_original_qc(new_out_dir)
    else:
        preprocesser = PreProcessIDAT(idat_dir)
        preprocesser.load_idats(geo_query)
        preprocesser.plot_original_qc(output_dir)

### TODO: Wrap a class around following functions ###

def print_case_controls():
    """Print number of case and controls for subtypes"""
    pass

def remove_controls():
    """Remove controls for study"""
    pass

def remove_low_sample_number():
    """Remove cases for study with low sample number"""
    pass
