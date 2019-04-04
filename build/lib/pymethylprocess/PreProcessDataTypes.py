"""
PreProcessDataTypes.py
======================
Contains datatypes core to downloading IDATs, preprocessing IDATs and samplesheets.
"""

import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr #from utils import importr_ as importr#
from rpy2.robjects import r
import rpy2.robjects.packages as rpackages
from pymethylprocess.MethylationDataTypes import *
from rpy2.robjects import pandas2ri, numpy2ri
import pickle
import sqlite3
import os, glob, subprocess
from collections import Counter

pandas2ri.activate()
numpy2ri.activate()


class TCGADownloader:
    """Downloads TCGA and GEO IDAT and clinical data"""
    def __init__(self):
        pass

    def download_tcga(self, output_dir):
        """Download TCGA IDATs.

        Parameters
        ----------
        output_dir
            Where to output idat files."""
        tcga = importr("TCGAbiolinks")
        print(tcga)
        robjects.r("""
                   library(TCGAbiolinks)
                   projects <- TCGAbiolinks:::getGDCprojects()$project_id
                   projects <- projects[grepl('^TCGA',projects,perl=T)]
                   match.file.cases.all <- NULL
                   for(proj in projects){
                        print(proj)
                        query <- GDCquery(project = proj,
                                          data.category = "Raw microarray data",
                                          data.type = "Raw intensities",
                                          experimental.strategy = "Methylation array",
                                          legacy = TRUE,
                                          file.type = ".idat",
                                          platform = "Illumina Human Methylation 450")
                        match.file.cases <- getResults(query,cols=c("cases","file_name"))
                        match.file.cases$project <- proj
                        match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
                        tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
                                 error = function(e) GDCdownload(query, method = "client"))
                    }
                    # This will create a map between idat file name, cases (barcode) and project
                    readr::write_tsv(match.file.cases.all, path = "idat_filename_case.txt")
                    # code to move all files to local folder
                    for(file in dir(".",pattern = ".idat", recursive = T)){
                        TCGAbiolinks:::move(file,file.path('%s',basename(file)))
                    }
                   """%output_dir)

    def download_clinical(self, output_dir):
        """Download TCGA Clinical Data.

        Parameters
        ----------
        output_dir
            Where to output clinical data csv."""
        robjects.r("""
                   library(TCGAbiolinks)
                   library(data.table)
                   projects <- TCGAbiolinks:::getGDCprojects()$project_id
                   projects <- projects[grepl('^TCGA',projects,perl=T)]
                   match.file.cases.all <- NULL
                   data <- list()
                   for(n in 1:length(projects)){
                        proj <- projects[n]
                        clin.query <- GDCquery_clinic(project = proj,
                                          type='clinical', save.csv=F)
                        data[[length(data)+1]] = clin.query
                    }
                    df <- rbindlist(data)
                    write.csv(df, file=file.path('%s','clinical_info.csv'))
                   """%output_dir)

    def download_geo(self, query, output_dir):
        """Download GEO IDATs.

        Parameters
        ----------
        query
            GEO accession number to query, must be 450k/850k.
        output_dir
            Output directory to store idats and clinical information csv"""
        base=importr('base')
        geo = importr("GEOquery")
        geo.getGEOSuppFiles(query)
        raw_tar_file = "{0}/{0}_RAW.tar".format(query)
        if not os.path.exists(raw_tar_file):
            print("Warning: GEO file {} not downloaded. Check accession on GEO and make sure there is this file available.".format(raw_tar_file))
        robjects.r["untar"](raw_tar_file, exdir = "{}/idat".format(query), tar='internal')
        idatFiles = robjects.r('list.files("{}/idat", pattern = "idat.gz$", full = TRUE)'.format(query))
        robjects.r["sapply"](idatFiles, robjects.r["gunzip"], overwrite = True)
        subprocess.call('mv {}/idat/*.idat {}/'.format(query, output_dir),shell=True)
        pandas2ri.ri2py(robjects.r['as'](robjects.r("phenoData(getGEO('{}')[[1]])".format(query)),'data.frame')).to_csv('{}/{}_clinical_info.csv'.format(output_dir,query))# ,GSEMatrix = FALSE


class PreProcessPhenoData:
    """Class that will manipute phenotype samplesheet before preprocessing of IDATs.

    pheno_sheet
        Location of clinical info csv.
    idat_dir
        Location of idats
    header_line
        Where to start reading clinical csv"""
    def __init__(self, pheno_sheet, idat_dir, header_line=0):
        self.xlsx = True if pheno_sheet.endswith('.xlsx') or pheno_sheet.endswith('.xls') else False
        if self.xlsx:
            self.pheno_sheet = pd.read_excel(pheno_sheet,header=header_line)
        else:
            self.pheno_sheet = pd.read_csv(pheno_sheet, header=header_line)
        self.idat_dir = idat_dir

    def format_geo(self, disease_class_column="methylation class:ch1", include_columns={}):
        """Format clinical sheets if downloaded geo idats.

        Parameters
        ----------
        disease_class_column
            Disease column of clinical info csv.
        include_columns
            Dictionary specifying other columns to include, and new names to assign them to."""
        idats = glob.glob("{}/*.idat".format(self.idat_dir))
        idat_basenames = np.unique(np.vectorize(lambda x: '_'.join(x.split('/')[-1].split('_')[:3]))(idats))
        idat_geo_map = dict(zip(np.vectorize(lambda x: x.split('_')[0])(idat_basenames),np.array(idat_basenames)))
        self.pheno_sheet['Basename'] = self.pheno_sheet['geo_accession'].map(idat_geo_map)
        self.pheno_sheet = self.pheno_sheet[self.pheno_sheet['Basename'].isin(idat_basenames)]
        self.pheno_sheet.loc[:,'Basename'] = self.pheno_sheet['Basename'].map(lambda x: self.idat_dir+x)
        col_dict = {'geo_accession':'AccNum',disease_class_column:'disease'}
        col_dict.update(include_columns)
        self.pheno_sheet = self.pheno_sheet[['Basename', 'geo_accession',disease_class_column]+(list(include_columns.keys()) if include_columns else [])].rename(columns=col_dict)

    def format_tcga(self, mapping_file="idat_filename_case.txt"):
        """Format clinical sheets if downloaded tcga idats.

        Parameters
        ----------
        mapping_file
            Maps uuids to proper tcga sample names, should be downloaded with tcga clinical information."""
        def decide_case_control(barcode):
            case_control_num = int(barcode.split('-')[3][:2])
            if case_control_num < 10:
                return 'case'
            elif case_control_num < 20:
                return 'normal'
            else:
                return 'control'
            return 0
        idats = glob.glob("{}/*.idat".format(self.idat_dir))
        barcode_mappings = pd.read_csv(mapping_file,sep='\t')
        barcode_mappings['barcodes'] = np.vectorize(lambda x: '-'.join(x.split('-')[:3]))(barcode_mappings['cases'])
        barcode_mappings['idats'] = barcode_mappings['file_name'].map(lambda x: x[:x.rfind('_')])
        barcode_mappings_d1 = dict(barcode_mappings[['barcodes','idats']].values.tolist())
        barcode_mappings['case_controls']= barcode_mappings['cases'].map(decide_case_control)
        barcode_mappings_d2 = dict(barcode_mappings[['barcodes','case_controls']].values.tolist())
        self.pheno_sheet['Basename'] = self.pheno_sheet['bcr_patient_barcode'].map(barcode_mappings_d1)
        self.pheno_sheet['case_control'] = self.pheno_sheet['bcr_patient_barcode'].map(barcode_mappings_d2)
        idat_basenames = np.unique(np.vectorize(lambda x: '_'.join(x.split('/')[-1].split('_')[:2]))(idats))
        self.pheno_sheet = self.pheno_sheet[self.pheno_sheet['Basename'].isin(idat_basenames)]
        self.pheno_sheet.loc[:,['Basename']] = self.pheno_sheet['Basename'].map(lambda x: self.idat_dir+x)
        self.pheno_sheet = self.pheno_sheet[['Basename', 'bcr_patient_barcode', 'disease', 'tumor_stage', 'vital_status', 'age_at_diagnosis', 'gender', 'race', 'ethnicity','case_control']].rename(columns={'tumor_stage':'stage','bcr_patient_barcode':'PatientID','vital_status':'vital','gender':'Sex','age_at_diagnosis':'age'})

    def format_custom(self, basename_col, disease_class_column, include_columns={}):
        """Custom format clinical sheet if user supplied idats.

        Parameters
        ----------
        basename_col
            Column name of sample names.
        disease_class_column
            Disease column of clinical info csv.
        include_columns
            Dictionary specifying other columns to include, and new names to assign them to.
        """
        idats = glob.glob("{}/*.idat".format(self.idat_dir))
        idat_basenames = np.unique(np.vectorize(lambda x: '_'.join(x.split('/')[-1].split('_')[:-1]))(idats))
        idat_count_underscores = np.vectorize(lambda x: x.count('_'))(idat_basenames)
        self.pheno_sheet['Basename'] = self.pheno_sheet[basename_col]
        basename_count_underscores = np.vectorize(lambda x: x.count('_'))(self.pheno_sheet['Basename'])
        min_underscores=min(np.hstack([idat_count_underscores,basename_count_underscores]))
        basic_basename_fn = np.vectorize(lambda x: '_'.join(x.split('_')[-min_underscores-1:]))
        basic_basename=dict(zip(basic_basename_fn(self.pheno_sheet['Basename']),self.pheno_sheet['Basename'].values))
        basic_idat=dict(zip(basic_basename_fn(idat_basenames),idat_basenames))
        complete_mapping={basic_basename[basename]:basic_idat[basename] for basename in basic_basename}
        self.pheno_sheet.loc[:,'Basename']=self.pheno_sheet['Basename'].map(complete_mapping).map(lambda x: self.idat_dir+x)
        self.pheno_sheet['disease'] = self.pheno_sheet[disease_class_column.replace("'",'')]
        self.pheno_sheet = self.pheno_sheet[np.unique(['Basename', 'disease']+list(include_columns.keys()))].rename(columns=include_columns)

    def merge(self, other_formatted_sheet, use_second_sheet_disease=True, no_disease_merge=False):
        """Merge multiple PreProcessPhenoData objects, merge their dataframes to accept more than one saplesheet/dataset or add more pheno info.

        Parameters
        ----------
        other_formatted_sheet
            Other PreProcessPhenoData to merge.
        use_second_sheet_disease
            Change disease column to that of second sheet instead of first.
        no_disease_merge
            Keep both disease columns from both sheets.
        """
        disease_dict = {False:'disease_x',True:'disease_y'}
        self.pheno_sheet = self.pheno_sheet.merge(other_formatted_sheet.pheno_sheet,how='inner', on='Basename')
        if not no_disease_merge:
            self.pheno_sheet['disease'] = self.pheno_sheet[disease_dict[use_second_sheet_disease]]
        cols=list(self.pheno_sheet)
        self.pheno_sheet = self.pheno_sheet[[col for col in cols if col!='Unnamed: 0_x' and col!='Unnamed: 0_y' and col!=disease_dict[use_second_sheet_disease]]]

    def concat(self, other_formatted_sheet):
        """Concat multiple PreProcessPhenoData objects, concat their dataframes to accept more than one smaplesheet/dataset.

        Parameters
        ----------
        other_formatted_sheet
            Other PreProcessPhenoData to concat.
        """
        self.pheno_sheet=pd.concat([self.pheno_sheet,other_formatted_sheet.pheno_sheet],join='inner').reset_index(drop=True)
        self.pheno_sheet=self.pheno_sheet[[col for col in list(self.pheno_sheet) if not col.startswith('Unnamed:')]]

    def export(self, output_sheet_name):
        """Export pheno data to csv after done with manipulation.

        Parameters
        ----------
        output_sheet_name
            Output csv name.
        """
        self.pheno_sheet.to_csv(output_sheet_name)
        print("Please move all other sample sheets out of this directory.")

    def split_key(self, key, subtype_delimiter):
        """Split pheno column by key, with subtype delimiter, eg. entry S1,s2 -> S1 with delimiter ",".

        Parameters
        ----------
        key
            Pheno column name.
        subtype_delimiter
            Subtype delimiter to split on.
        """
        new_key = '{}_only'.format(key)
        self.pheno_sheet[new_key] = self.pheno_sheet[key].map(lambda x: x.split(subtype_delimiter)[0])
        return new_key

    def get_categorical_distribution(self, key, disease_only=False, subtype_delimiter=','):
        """Print categorical distribution, counts for each unique value in phenotype column.

        Parameters
        ----------
        key
            Phenotype Column.
        disease_only
            Whether to split phenotype column entries by delimiter.
        subtype_delimiter
            Subtype delimiter to split on.
        """
        if type(key) == type('string'):
            if disease_only:
                key=self.split_key(key,subtype_delimiter)
            return Counter(self.pheno_sheet[key])
        else:
            cols=self.pheno_sheet[list(key)].astype(str)
            cols=reduce(lambda a,b: a+'_'+b,[cols.iloc[:,i] for i in range(cols.shape[1])])
            return Counter(cols)

    def remove_diseases(self,exclude_disease_list, low_count, disease_only,subtype_delimiter):
        """Remove samples with certain diseases from disease column.

        Parameters
        ----------
        exclude_disease_list
            List containing diseases to remove.
        low_count
            Remove samples that have less than x disease occurances in column.
        disease_only
            Whether to split phenotype column entries by delimiter.
        subtype_delimiter
            Subtype delimiter to split on.
        """
        if low_count:
            low_count = int(low_count)
            cat_dist = self.get_categorical_distribution('disease', disease_only,subtype_delimiter).items()
            count_diseases=pd.DataFrame(list(cat_dist),columns=['disease','count'])
            count_diseases.loc[:,'count'] = count_diseases.loc[:,'count'].astype(int)
            exclude_diseases_more=count_diseases.loc[count_diseases['count'].values<low_count,'disease']
            exclude_diseases_more=exclude_diseases_more.unique().tolist()
            if disease_only:
                exclude_diseases_more=self.pheno_sheet.loc[self.pheno_sheet['disease_only'].isin(exclude_diseases_more),'disease'].unique().tolist()
        else:
            exclude_diseases_more=[]
        self.pheno_sheet = self.pheno_sheet[~self.pheno_sheet['disease'].isin(exclude_disease_list+exclude_diseases_more)]

class PreProcessIDAT:
    """Class that will preprocess IDATs using R pipelines.

    idat_dir
        Location of idats or samplesheet csv.
    minfi
        Rpy2 importr minfi library, default to None will load through rpy2
    enmix
        Rpy2 importr enmix library, default to None will load through rpy2
    base
        Rpy2 importr base library, default to None will load through rpy2
    meffil
        Rpy2 importr meffil library, default to None will load through rpy2"""
    def __init__(self, idat_dir, minfi=None, enmix=None, base=None, meffil=None):
        self.idat_dir = idat_dir
        if minfi == None:
            self.minfi = importr('minfi')
        else:
            self.minfi = minfi
        if enmix == None:
            self.enmix = importr("ENmix")
        else:
            self.enmix = enmix
        if base == None:
            self.base = importr('base')
        else:
            self.base = base
        try:
            if meffil==None:
                self.meffil = importr('meffil')
            else:
                self.meffil = meffil
        except:
            self.meffil=None
        self.qcinfo=robjects.r('NULL')

    def move_jpg(self):
        """Move jpeg files from current working directory to the idat directory.
        """
        subprocess.call('mv *.jpg {}'.format(self.idat_dir),shell=True)

    def preprocessNoob(self):
        """Run minfi preprocessing with Noob normalization"""
        self.qcinfo = self.enmix.QCinfo(self.RGset, detPthre=1e-7)
        self.move_jpg()
        self.MSet = self.minfi.preprocessNoob(self.RGset)
        self.move_jpg()
        self.MSet = self.enmix.QCfilter(self.MSet,qcinfo=self.qcinfo,outlier=True)
        self.move_jpg()
        return self.MSet

    def preprocessRAW(self):
        """Run minfi preprocessing with RAW normalization"""
        self.qcinfo = self.enmix.QCinfo(self.RGset, detPthre=1e-7)
        self.move_jpg()
        self.MSet = self.minfi.preprocessRaw(self.RGset)
        self.move_jpg()
        self.MSet = self.enmix.QCfilter(self.MSet,qcinfo=self.qcinfo,outlier=True)
        self.move_jpg()
        return self.MSet

    def preprocessENmix(self, n_cores=6):
        """Run ENmix preprocessing pipeline.

        Parameters
        ----------
        n_cores
            Number of CPUs to use."""
        self.qcinfo = self.enmix.QCinfo(self.RGset, detPthre=1e-7)
        self.move_jpg()
        self.MSet = self.enmix.preprocessENmix(self.RGset, QCinfo=self.qcinfo, nCores=n_cores)
        self.move_jpg()
        self.MSet = self.enmix.QCfilter(self.MSet,qcinfo=self.qcinfo,outlier=True)
        self.move_jpg()
        return self.MSet

    def preprocessMeffil(self, n_cores=6, n_pcs=4, qc_report_fname="qc/report.html", normalization_report_fname='norm/report.html', pc_plot_fname='qc/pc_plot.pdf', useCache=True, qc_only=True, qc_parameters={'p.beadnum.samples':0.1,'p.detection.samples':0.1,'p.detection.cpgs':0.1,'p.beadnum.cpgs':0.1}, rm_sex=False):
        """Run meffil preprocessing pipeline with functional normalization.

        Parameters
        ----------
        n_cores
            Number of CPUs to use.
        n_pcs
            Number of principal components to use for functional normalization, set to -1 to autoselect via kneedle algorithm.
        qc_report_fname
            HTML filename to store QC report.
        normalization_report_fname
            HTML filename to store normalization report
        pc_plot_fname
            PDF file to store principal components plot.
        useCache
            Use saved QC objects instead of running through QC again.
        qc_only
            Perform QC, then save and quit before normalization.
        qc_parameters
            Python dictionary with parameters for qc.
        rm_sex
            Remove non-autosomal cpgs?"""
        from pymethylprocess.meffil_functions import load_detection_p_values_beadnum, set_missing, remove_sex
        self.pheno = self.meffil.meffil_read_samplesheet(self.idat_dir, verbose=True)
        cache_storage_path = os.path.join(self.idat_dir,'QCObjects.rds')
        qc_parameters=robjects.r("""function(p.beadnum.samples,p.detection.samples,p.detection.cpgs,p.beadnum.cpgs,sex.outlier.sd){
                        qc.parameters <- meffil.qc.parameters(
                            	beadnum.samples.threshold             = p.beadnum.samples,
                            	detectionp.samples.threshold          = p.detection.samples,
                            	detectionp.cpgs.threshold             = p.detection.cpgs,
                            	beadnum.cpgs.threshold                = p.beadnum.cpgs,
                            	sex.outlier.sd                        = sex.outlier.sd,
                            	snp.concordance.threshold             = 0.95,
                            	sample.genotype.concordance.threshold = 0.8
                            )
                        return(qc.parameters)
                        }""")(qc_parameters['p.beadnum.samples'], qc_parameters['p.detection.samples'], qc_parameters['p.detection.cpgs'], qc_parameters['p.beadnum.cpgs'], qc_parameters['sex.outlier.sd'])
        if useCache:
            qc_list = robjects.r('readRDS')(cache_storage_path)
            qc_list = robjects.r("""function(qc.list,qc.parameters, qc.report.fname) {
                                qc.list$qc.summary <- meffil.qc.summary(qc.list$qc.objects,parameters=qc.parameters,verbose=F)
                                meffil.qc.report(qc.list$qc.summary, output.file=qc.report.fname)
                                return(qc.list)
                                }""")(qc_list,qc_parameters, qc_report_fname)
        else:
            qc_list = robjects.r("""function(samplesheet,n.cores,qc.parameters,qc.report.fname){
                qc.objects<-meffil.qc(samplesheet,mc.cores=n.cores,detection.threshold=0.000001,verbose=F)
                qc.summary<-meffil.qc.summary(qc.objects,parameters=qc.parameters,verbose=F)
                meffil.qc.report(qc.summary, output.file=qc.report.fname)
                return(list(qc.objects=qc.objects,qc.summary=qc.summary))
                }""")(self.pheno, n_cores, qc_parameters, qc_report_fname) # p.values=meffil.load.detection.pvalues(qc.objects)
        pc_df = robjects.r("""function(qc.list,pc.plot.fname){
            y <- meffil.plot.pc.fit(qc.list$qc.objects)
            print(y)
            ggsave(y$plot,filename=pc.plot.fname,height=6,width=6)
            return(y$data)
            }""")(qc_list,pc_plot_fname)
        print("Check QC report and select number of PCs. Will add option in future to adjust thresholds.")
        if n_pcs == -1:
            from kneed import KneeLocator
            pc_df=pandas2ri.ri2py(robjects.r['as'](pc_df,'data.frame'))
            pc_df['B'] = pc_df['U']+pc_df['M']
            n_pcs=int(KneeLocator(pc_df['n'].values, pc_df['B'].values, S=1.0, curve='convex', direction='decreasing').knee)
            with open(pc_plot_fname.replace('.pdf','.txt'),'w') as f:
                f.write('pcs_selected:{}'.format(n_pcs))
        if qc_only:
            robjects.r('saveRDS')(qc_list,cache_storage_path)
            exit()
        if not os.path.exists(cache_storage_path):
            robjects.r('saveRDS')(qc_list,cache_storage_path)
        pval_beadnum=load_detection_p_values_beadnum(qc_list,n_cores)
        #print(pval_beadnum)
        self.beta_final = robjects.r("""function(qc.list, n.pcs, norm.report.fname,mc.cores) {
            options(mc.cores=mc.cores)
            qc.objects = qc.list$qc.objects
            qc.summary = qc.list$qc.summary
            outlier <- qc.summary$bad.samples
            if (nrow(outlier) > 0) {
            table(outlier$issue)
            index <- outlier$issue %in% c("Control probe (dye.bias)",
                                          "Methylated vs Unmethylated",
                                          "X-Y ratio outlier",
                                          "Low bead numbers",
                                          "Detection p-value",
                                          "Sex mismatch",
                                          "Genotype mismatch",
                                          "Control probe (bisulfite1)",
                                          "Control probe (bisulfite2)")
            outlier <- outlier[index,]
            if (nrow(outlier) > 0) {
                qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
            }
            }
            norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=n.pcs, verbose=F)
            norm <- meffil.normalize.samples(norm.objects, just.beta=F, cpglist.remove=qc.summary$bad.cpgs$name)
            beta <- meffil.get.beta(norm$M, norm$U)
            pcs <- meffil.methylation.pcs(beta)
            norm.summary <- meffil.normalization.summary(norm.objects, pcs=pcs)
            meffil.normalization.report(norm.summary, output.file=norm.report.fname)
            return(beta)}""")(qc_list, n_pcs, normalization_report_fname,n_cores)

        self.beta_final = set_missing(self.beta_final, pval_beadnum)

        if rm_sex:
            self.beta_final = remove_sex(self.beta_final)

        try:
            grdevice = importr("grDevices")
            geneplotter = importr("geneplotter")
            qr_report_fname=os.path.abspath(qc_report_fname).split('/')
            qr_report_fname[-1]='beta_dist.jpg'
            grdevice.jpeg('/'.join(qr_report_fname),height=900,width=600)
            base.par(mfrow=robjects.vectors.IntVector([1,2]))
            self.enmix.multidensity(self.beta_final, main="Multidensity")
            self.enmix.multifreqpoly(self.beta_final, xlab="Beta value")
            grdevice.dev_off()
        except:
            pass

    def load_idats(self):
        """For minfi pipeline, load IDATs from specified idat_dir."""
        targets = self.minfi.read_metharray_sheet(self.idat_dir)
        self.RGset = self.minfi.read_metharray_exp(targets=targets, extended=True)
        return self.RGset

    def plot_qc_metrics(self, output_dir):
        """Plot QC results from ENmix pipeline and possible minfi. Still experimental.

        Parameters
        ----------
        output_dir
            Where to store plots."""
        self.enmix.plotCtrl(self.RGset)
        grdevice = importr("grDevices")
        geneplotter = importr("geneplotter")
        base = importr('base')
        anno=self.minfi.getAnnotation(self.RGset)
        anno_py = pandas2ri.ri2py(robjects.r['as'](anno,'data.frame'))
        beta_py = pandas2ri.ri2py(self.beta)
        beta1=numpy2ri.py2ri(beta_py[anno_py["Type"]=="I"])
        beta2=numpy2ri.py2ri(beta_py[anno_py["Type"]=="II"])
        grdevice.jpeg(output_dir+'/dist.jpg',height=900,width=600)
        base.par(mfrow=robjects.vectors.IntVector([3,2]))
        self.enmix.multidensity(self.beta, main="Multidensity")
        self.enmix.multifreqpoly(self.beta, xlab="Beta value")
        self.enmix.multidensity(beta1, main="Multidensity: Infinium I")
        self.enmix.multifreqpoly(beta1, main="Multidensity: Infinium I", xlab="Beta value")
        self.enmix.multidensity(beta2, main="Multidensity: Infinium II")
        self.enmix.multifreqpoly(beta2, main="Multidensity: Infinium II", xlab="Beta value")
        grdevice.dev_off()
        self.minfi.qcReport(self.RGset, pdf = "{}/qcReport.pdf".format(output_dir))
        self.minfi.mdsPlot(self.RGset)
        self.minfi.densityPlot(self.RGset, main='Beta', xlab='Beta')

    def return_beta(self):
        """Return minfi RSet after having created MSet."""
        self.RSet = self.minfi.ratioConvert(self.MSet, what = "both", keepCN = True)
        return self.RSet

    def get_beta(self):
        """Get beta value matrix from minfi after finding RSet."""
        self.beta = self.minfi.getBeta(self.RSet)
        return self.beta

    def filter_beta(self):
        """After creating beta, filter out outliers."""
        self.beta_final=self.enmix.rm_outlier(self.beta,qcscore=self.qcinfo)
        return self.beta_final

    def get_meth(self):
        """Get methylation intensity matrix from MSet"""
        return self.minfi.getMeth(self.MSet)

    def get_unmeth(self):
        """Get unmethylated intensity matrix from MSet"""
        return self.minfi.getUnmeth(self.MSet)

    def extract_pheno_data(self, methylset=False):
        """Extract pheno data from MSet or RGSet, minfi.

        Parameters
        ----------
        methylset
            If MSet has beenn created, set to True, else extract from original RGSet."""
        self.pheno = robjects.r("pData")(self.MSet) if methylset else robjects.r("pData")(self.RGset)
        return self.pheno

    def extract_manifest(self):
        """Get manifest from RGSet."""
        self.manifest = self.minfi.getManifest(self.RGset)
        return self.manifest

    def preprocess_enmix_pipeline(self, n_cores=6, pipeline='enmix', noob=False, qc_only=False, use_cache=False):
        """Run complete ENmix or minfi preprocessing pipeline.

        Parameters
        ----------
        n_cores
            Number CPUs.
        pipeline
            Run enmix or minfi
        noob
            Noob norm or RAW if minfi running.
        qc_only
            Save and quit after only running QC?
        use_cache
            Load preexisting RGSet instead of running QC again."""
        cache_storage_path = os.path.join(self.idat_dir,'RGSet.rds')
        if use_cache:
            self.RGset=robjects.r('readRDS')(cache_storage_path)
        else:
            self.load_idats()
        if qc_only:
            robjects.r('saveRDS')(self.RGset,cache_storage_path)
            exit()
        if not os.path.exists(cache_storage_path):
            robjects.r('saveRDS')(self.RGset,cache_storage_path)
        if pipeline =='enmix':
            self.preprocessENmix(n_cores)
        else:
            if noob:
                self.preprocessNoob()
            else:
                self.preprocessRAW()
        self.return_beta()
        self.get_beta()
        self.filter_beta()
        self.extract_pheno_data(methylset=True)
        return self.pheno, self.beta_final

    def plot_original_qc(self, output_dir):
        """Plot QC results from ENmix pipeline and possible minfi. Still experimental.

        Parameters
        ----------
        output_dir
            Where to store plots."""
        self.preprocessRAW()
        self.return_beta()
        self.get_beta()
        self.plot_qc_metrics(output_dir)

    def output_pheno_beta(self, meffil=False):
        """Get pheno and beta dataframe objects stored as attributes for input to MethylationArray object.

        Parameters
        ----------
        meffil
            True if ran meffil pipeline."""
        self.pheno_py=pandas2ri.ri2py(robjects.r['as'](self.pheno,'data.frame'))
        if not meffil:
            self.beta_py=pd.DataFrame(pandas2ri.ri2py(self.beta_final),index=numpy2ri.ri2py(robjects.r("featureNames")(self.RSet)),columns=numpy2ri.ri2py(robjects.r("sampleNames")(self.RSet))).transpose()
            self.pheno_py['Sample_Name']=np.vectorize(lambda x: x.split('/')[-1])(self.pheno_py['Basename'])
            self.pheno_py = self.pheno_py.set_index('Sample_Name').loc[self.beta_py.index,:]
        else:
            self.beta_py=pd.DataFrame(pandas2ri.ri2py(self.beta_final),index=robjects.r("rownames")(self.beta_final),columns=robjects.r("colnames")(self.beta_final)).transpose()
            print(self.beta_py)
            print(self.beta_py.index)
            print(self.pheno_py)
            self.pheno_py = self.pheno_py.set_index('Sample_Name').loc[self.beta_py.index,:]

    def export_pickle(self, output_pickle, disease=''):
        """Export pheno and beta dataframes to pickle, stored in python dict that can be loaded into MethylationArray

        Parameters
        ----------
        output_pickle
            Where to store MethylationArray.
        disease
            Custom naming scheme for data."""
        output_dict = {}
        if os.path.exists(output_pickle):
            output_dict = pickle.load(open(output_pickle,'rb'))
        output_dict['pheno' if not disease else 'pheno_{}'.format(disease)] = self.pheno_py
        output_dict['beta' if not disease else 'beta_{}'.format(disease)] = self.beta_py
        pickle.dump(output_dict, open(output_pickle,'wb'),protocol=4)

    def export_sql(self, output_db, disease=''):
        """Export pheno and beta dataframes to SQL

        Parameters
        ----------
        output_db
            Where to store data, sqlite db.
        disease
            Custom naming scheme for data."""
        conn = sqlite3.connect(output_db)
        self.pheno_py.to_sql('pheno' if not disease else 'pheno_{}'.format(disease), con=conn, if_exists='replace')
        self.beta_py.to_sql('beta' if not disease else 'beta_{}'.format(disease), con=conn, if_exists='replace')
        conn.close()

    def export_csv(self, output_dir):
        """Export pheno and beta dataframes to CSVs

        Parameters
        ----------
        output_dir
            Where to store csvs."""
        self.pheno_py.to_csv('{}/pheno.csv'.format(output_dir))
        self.beta_py.to_csv('{}/beta.csv'.format(output_dir))

    def to_methyl_array(self,disease=''):
        """Convert results from preprocessing into MethylationArray, and directly return MethylationArray object.

        Parameters
        ----------
        disease
            Custom naming scheme for data."""
        return MethylationArray(self.pheno_py,self.beta_py, disease)
