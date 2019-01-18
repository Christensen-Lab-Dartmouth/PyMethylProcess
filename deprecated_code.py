

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




#print(self.beta_final)
if 0:
    qc_objects = self.meffil.meffil_qc(self.pheno, mc_cores=n_cores, verbose=False) # , number_quantiles=500, detection_threshold=0.01, bead_threshold=3, sex_cutoff=-2, chip="450k",
    #print(qc_objects)
    robjects.r('saveRDS')(qc_objects,'{}/r_obj.rds'.format(self.idat_dir))
    # Generate QC report
    qc_summary = self.meffil.meffil_qc_summary(qc_objects, verbose=False)

    #print(qc_summary)

    #self.meffil.meffil_qc_report(qc_summary, output_file="qc/report.html")

    # Remove outlier samples if necessary
    qc_objects = self.meffil.meffil_remove_samples(qc_objects, dollar(dollar(qc_summary,'bad.samples'),'sample.name'))

    #print(qc_objects)

    # Plot residuals remaining after fitting control matrix to decide on the number PCs
    # to include in the normalization below.
    #print(self.meffil.meffil_plot_pc_fit(qc_objects)$plot)

    # Perform quantile normalization
    norm_objects = self.meffil.meffil_normalize_quantiles(qc_objects, number_pcs=n_pcs, mc_cores=n_cores, verbose=False)

    # Generate normalized probe values
    norm_beta = self.meffil.meffil_normalize_samples(norm_objects, just_beta=True, mc_cores=n_cores, cpglist_remove=dollar(dollar(qc_summary,'bad.cpgs'),'name'))
    #beta <- meffil.get.beta(norm.dataset$M, norm.dataset$U)
    # Generate normalization report
    #pcs = self.meffil.meffil_methylation_pcs(norm_beta)
    #norm_summary = self.meffil.meffil_normalization_summary(norm_objects, pcs=pcs)
    #self.meffil.meffil_normalization_report(norm_summary, output_file="normalization/report.html")

    self.beta_final = dollar(norm_beta,'beta')

    #self.beta_final = dollar(self.meffil.meffil_normalize_dataset(self.pheno, qc_file="qc/report.html", author="", study="Illumina450", number_pcs=n_pcs, mc_cores=n_cores, verbose=True),'beta')#10
    #robjects.r('saveRDS')(self.beta_final,'r_obj.rds')
    #print(numpy2ri.ri2py(robjects.r("colnames")(self.beta_final)))
    #print(self.beta_final.slots)
    #self.beta_final = robjects.r['as'](self.beta_final,'data.frame'))
    #print(robjects.r['as'](self.beta_final,'data.frame'))
    #b=pandas2ri.ri2py(robjects.r['as'](self.beta_final,'data.frame'))
    #print(b)
    #print(pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame'))[b.index])
