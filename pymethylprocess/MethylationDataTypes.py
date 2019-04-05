"""
MethylationDataTypes.py
=======================
Contains datatypes core to storing beta and phenotype methylation data, and imputation.
"""

import pandas as pd, numpy as np
import pickle, os

class MethylationArray:
    """Stores beta and phenotype information and performs various operations.
    Initialize MethylationArray object by inputting dataframe of phenotypes
    and dataframe of beta values with samples as index.

    pheno_df
        Phenotype dataframe (samples x covariates)
    beta_df
        Beta Values Dataframe (samples x cpgs)"""
    def __init__(self, pheno_df, beta_df, name=''):
        """"""
        self.pheno=pheno_df
        self.beta=beta_df
        self.name=name

    def export(self, output_pickle):
        pass

    def write_csvs(self, output_dir):
        """Write phenotype data and beta values to csvs.

        Parameters
        ----------
        output_dir
            Directory to output csv files.
        """
        self.pheno.to_csv('{}/pheno.csv'.format(output_dir))
        self.beta.to_csv('{}/beta.csv'.format(output_dir))

    def write_pickle(self, output_pickle, disease=''):
        """Store phenotype data and beta values in pickle file.
        Is default file format for storing MethylationArray objects.

        Parameters
        ----------
        output_pickle
            Pickle file to store MethylationArray data.
        """
        output_dict = {}
        if 0 and os.path.exists(output_pickle):
            output_dict = pickle.load(open(output_pickle,'rb'))
        output_dict['pheno'] = self.pheno #  if not disease else 'pheno_{}'.format(disease)
        output_dict['beta'] = self.beta #  if not disease else 'beta_{}'.format(disease)
        pickle.dump(output_dict, open(output_pickle,'wb'),protocol=4)

    def write_db(self, conn, disease=''):
        """Store phenotype data and beta values in SQL database.

        Parameters
        ----------
        conn
            SQLite connection.
        disease
            Create new tables in db that are related to disease state by this name.
        """
        self.pheno.to_sql('pheno' if not disease else 'pheno_{}'.format(disease), con=conn, if_exists='replace')
        self.beta.to_sql('beta' if not disease else 'beta_{}'.format(disease), con=conn, if_exists='replace')

    def impute(self, imputer):
        """Perform imputation on NaN beta vaues. Input imputater returned from ImputerObject.

        Parameters
        ----------
        imputer
            Type of imputer object, in sklearn type interface.
        """
        self.beta = pd.DataFrame(imputer.fit_transform(self.beta),index=self.beta.index,columns=list(self.beta))

    def return_shape(self):
        """Return dimensionality and number of samples of beta matrix."""
        return self.beta.shape

    def bin_column(self, col, n_bins):
        """Turn continuous variable/covariate into categorical bins.
        Returns name of new column and updates phenotype matrix to reflect this change.

        Parameters
        ----------
        col
            Continuous column of phenotype array to bin.
        n_bins
            Number of bins to create.
        """
        new_col_name='{}_{}'.format(col,'binned')
        self.pheno[new_col_name]=np.vectorize(lambda x: x.replace(' ',''))(pd.cut(self.pheno[col],n_bins).astype(str))
        return new_col_name

    def split_train_test(self, train_p=0.8, stratified=True, disease_only=False, key='disease', subtype_delimiter=',', val_p=0.): # https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html
        """Split MethylationArray into training and test sets, with option to stratify by categorical covariate.

        Parameters
        ----------
        train_p
            Fraction of methylation array to use as training set.
        stratified
            Whether to stratify by categorical variable.
        disease_only
            Consider disease superclass by some delimiter. For instance if disease is S1,s2, superclass would be S1.
        key
            Column to stratify on.
        subtype_delimiter
            How to split disease column into super/subclass.
        val_p
            If set greater than 0, will create additional validation set, fraction of which is broken off from training set.
        """
        np.random.seed(42)
        if stratified:
            train_pheno=self.pheno.groupby(self.split_key(key,subtype_delimiter) if disease_only else key,group_keys=False).apply(lambda x: x.sample(frac=train_p))
        else:
            train_pheno = self.pheno.sample(frac=train_p)
        test_pheno = self.pheno.drop(train_pheno.index)
        train_arr, test_arr = MethylationArray(train_pheno,self.beta.loc[train_pheno.index.values,:],'train'),MethylationArray(test_pheno,self.beta.loc[test_pheno.index.values,:],'test')
        if val_p > 0:
            train_arr, val_arr = train_arr.split_train_test(1.-val_p,stratified,disease_only,key,subtype_delimiter)
            return train_arr, test_arr, val_arr
        else:
            return train_arr, test_arr

    def subsample(self, key='disease', n_samples=None, frac=None, categorical=False):
        """Subsample MethylationArray, make the set randomly smaller.

        Parameters
        ----------
        key
            If stratifying, use this column of pheno array.
        n_samples
            Number of samples to consider overall, or per stratum.
        frac
            Alternative to n_samples, where x frac of array or stratum is considered.
        categorical
            Whether to stratify by column.
        """
        np.random.seed(42)
        if categorical:
            pheno=self.pheno.groupby(key,group_keys=False).apply(lambda x: x.sample(frac=frac, n=min(x.shape[0],n_samples)))
        else:
            pheno=self.pheno.sample(frac=frac, n=min(n_samples,self.pheno.shape[0]))
        return MethylationArray(pheno,self.beta.loc[pheno.index.values,:],'subsampled')

    def split_key(self, key, subtype_delimiter):
        """Manipulate an entire phenotype column, splitting each element up by some delimiter.

        Parameters
        ----------
        key
            Phenotype column.
        subtype_delimiter
            How to break up strings in columns. S1,s2 -> S1 for instance.
        """
        new_key = '{}_only'.format(key)
        self.pheno.loc[:,new_key] = self.pheno[key].map(lambda x: x.split(subtype_delimiter)[0])
        return new_key

    def remove_whitespace(self, key):
        """Remove whitespaces from phenotype column.

        Parameters
        ----------
        key
            Phenotype column.
        """
        self.pheno.loc[:,key]=self.pheno[key].map(lambda x: x.replace(' ',''))

    def split_by_subtype(self, disease_only, subtype_delimiter):
        """Split MethylationArray into generator of MethylationArrays by phenotype column. Much akin to groupby. Only splits from disease column.

        Parameters
        ----------
        disease_only
            Consider disease superclass.
        subtype_delimiter
            How to break up disease column if using disease_only.
        """
        for disease, pheno_df in self.pheno.groupby(self.split_key('disease',subtype_delimiter) if disease_only else 'disease'):
            new_disease_name = disease.replace(' ','')
            beta_df = self.beta.loc[pheno_df.index,:]
            yield MethylationArray(pheno_df,beta_df,new_disease_name)

    def feature_select(self, n_top_cpgs, feature_selection_method='mad', metric='correlation', nn=10):
        """Perform unsupervised feature selection on MethylationArray.

        Parameters
        ----------
        n_top_cpgs
            Number of CpGs to retain.
        feature_selection_method
            Method to perform selection.
        metric
            If considering structural feature selection like SPEC, use this distance metric.
        nn
            Number of nearest neighbors.
        """
        if feature_selection_method=='mad':
            mad_cpgs = self.beta.mad(axis=0).sort_values(ascending=False)
            top_mad_cpgs = np.array(list(mad_cpgs.iloc[:n_top_cpgs].index))
            self.beta = self.beta.loc[:, top_mad_cpgs]
        elif feature_selection_method=='spectral':
            from pynndescent.pynndescent_ import PyNNDescentTransformer
            from skfeature.function.similarity_based.SPEC import spec
            if nn:
                W=PyNNDescentTransformer(n_neighbors=nn,metric=metric).fit_transform(self.beta.values)
                print('W computed...')
                f_weights = spec(self.beta.values,W=W)
                print('weights',f_weights)
            else:
                f_weights = spec(self.beta.values)
            self.beta = self.beta.iloc[:, np.argsort(f_weights)[::-1][:n_top_cpgs]]

    def remove_missingness(self, cpg_threshold=None, sample_threshold=None):
        """Remove samples and CpGs with certain level of missingness..

        Parameters
        ----------
        cpg_threshold
            If more than fraction of Samples for this CpG are missing, remove cpg.
        sample_threshold
            If more than fraction of CpGs for this sample are missing, remove sample.
        """
        cpgs_remove = np.array([])
        samples_remove = np.array([])
        cpgs=self.return_cpgs()
        samples=self.return_idx()
        na_frame = pd.isna(self.beta).astype(int)
        if cpg_threshold != None:
            cpgs_remove=cpgs[na_frame.sum(axis=0).apply(lambda x: x/float(na_frame.shape[0])).values >= cpg_threshold]
        if sample_threshold != None:
            samples_remove=samples[na_frame.sum(axis=1).apply(lambda x: x/float(na_frame.shape[1])) >= sample_threshold]
        self.beta = self.beta.loc[~np.isin(samples,samples_remove),~np.isin(cpgs,cpgs_remove)]
        self.pheno = self.pheno.loc[self.beta.index,:]

    def remove_na_samples(self, outcome_cols):
        """Remove samples of MethylationArray who have missing values in phenotype column.

        Parameters
        ----------
        outcome_cols
            Phenotype columns, if any rows contain missing values, samples are removed.
        """
        vals=self.pheno[outcome_cols].values
        if len(vals.shape) < 2:
            vals = vals[:,np.newaxis]
        remove_bool = ~np.isnan(vals).any(axis=1).flatten()
        self.pheno=self.pheno[remove_bool]
        self.beta=self.beta[remove_bool]

    def subset_index(self,index):
        """Subset MethylationArray by samples.

        Parameters
        ----------
        index
            Sample names to subset by.
        """
        return MethylationArray(self.pheno.loc[index,:],self.beta.loc[index,:])

    def return_idx(self):
        """Return sample names of MethylationArray.
        """
        return np.array(list(self.pheno.index))

    def return_cpgs(self):
        """Return list of cpgs of MethylationArray"""
        return np.array(list(self.beta))

    def return_raw_beta_array(self):
        """Return numpy array of methylation beta vaues."""
        return self.beta.values

    def subset_cpgs(self,cpgs):
        """Subset beta matrix by list of Cpgs.
            Parameters
            ----------
            cpgs
                Cpgs to subset by.
        """
        return MethylationArray(self.pheno,self.beta.loc[:,cpgs])

    def categorical_breakdown(self, key):
        """Print categorical distribution, counts for each unique value in phenotype column.

        Parameters
        ----------
        key
            Phenotype Column.
        """
        from collections import Counter
        print('\n'.join(['{}:{}'.format(k,v) for k,v in Counter(self.pheno[key].values).items()]))

    def overwrite_pheno_data(self, preprocess_sample_df):
        """Feed in another phenotype dataframe that will overwrite overlapping keys of existing phenotype array.

        Parameters
        ----------
        preprocess_sample_df
            New phenotype dataframe.
        """
        preprocess_sample_df=preprocess_sample_df.loc[self.pheno.index,:]
        for col in list(preprocess_sample_df):
            if col in list(self.pheno):
                self.pheno.loc[:,col]=preprocess_sample_df[col]
            else:
                self.pheno[col]=preprocess_sample_df[col]

    def merge_preprocess_sheet(self, preprocess_sample_df):
        """Feed in another phenotype dataframe that will be merged with existing phenotype array.

        Parameters
        ----------
        preprocess_sample_df
            New phenotype dataframe.
        """
        self.pheno=self.pheno.merge(preprocess_sample_df,on=['Basename'],how='inner')
        if 'disease_x' in list(self.pheno):
            self.pheno = self.pheno.rename(columns={'disease_x':'disease'})
        self.pheno = self.pheno[[col for col in list(self.pheno) if not col.startswith('Unnamed:')]]
        self.pheno=self.pheno.set_index([np.vectorize(lambda x: x.split('/')[-1])(self.pheno['Basename'])],drop=False)

    def groupby(self, key):
        """Groupby for Methylation Array. Returns generator of methylation arrays grouped by key.

        Parameters
        ----------
        preprocess_sample_df
            New phenotype dataframe.
        """
        for name, df in self.pheno.groupby(key):
            yield (str(name), MethylationArray(pheno_df=self.pheno.loc[df.index],beta_df=self.beta.loc[df.index]))

    @classmethod
    def from_pickle(self, input_pickle):
        """Load MethylationArray stored in pickle.

        Usage: MethylationArray.from_pickle([input_pickle])

        Parameters
        ----------
        input_pickle
            Stored MethylationArray pickle.
        """
        return MethylationArray(*extract_pheno_beta_df_from_pickle_dict(pickle.load(open(input_pickle,'rb'))))

class MethylationArrays:
    """Literally a list of methylation arrays, with methods operate on these arrays that is memory efficient.
    Initialize with list of methylation arrays. Can optionally leave list empty or with one element.

    list_methylation_arrays
        List of methylation arrays."""
    def __init__(self, list_methylation_arrays):
        """
        """
        self.methylation_arrays = list_methylation_arrays

    def __len__(self):
        """Number of stored methylation arrays."""
        return len(self.methylation_arrays)

    def combine(self, array_generator=None): # FIXME add sort based on samples
        """Combine the list of methylation arrays into one array via concatenation of beta matrices and phenotype arrays.

        Parameters
        ----------
        array_generator
            Generator of additional methylation arrays for computational memory minimization.
        """
        if len(self.methylation_arrays)> 1:
            pheno_df=pd.concat([methylArr.pheno for methylArr in self.methylation_arrays], join='inner')#.sort()
            beta_df=pd.concat([methylArr.beta for methylArr in self.methylation_arrays], join='inner')#.sort()
        else:
            pheno_df=self.methylation_arrays[0].pheno
            beta_df=self.methylation_arrays[0].beta
        if array_generator != None:
            for methylArr in array_generator:
                pheno_df=pd.concat([pheno_df,methylArr.pheno], join='inner')
                beta_df=pd.concat([beta_df,methylArr.beta], join='inner')
        return MethylationArray(pheno_df,beta_df)

    def write_dbs(self, conn):
        """Write list of methylation arrays to SQL database. Recommend naming MethylationArray.

        Parameters
        ----------
        conn
            SQL connection.
        """
        for methyl_arr in self.methylation_arrays:
            methyl_arr.write_db(conn, methyl_arr.name)

    def write_pkls(self,pkl):
        """Write list of methylation arrays to single pickle. Recommend naming each MethylationArray.

        Parameters
        ----------
        pkl
            Pickle file to write to.
        """
        for methyl_arr in self.methylation_arrays:
            methyl_arr.write_pickle(pkl, methyl_arr.name)

    def impute(self, imputer):
        """Impute all methylation arrays.

        Parameters
        ----------
        imputer
            Type of imputation, sklearn-like.
        """
        for i in range(len(self.methylation_arrays)):
            self.methylation_arrays[i].impute(imputer)

class ImputerObject:
    """Class that stores and accesses different types of imputers.
    Construct sklearn-like imputer given certain input arguments.

    solver
        Library for imputation, eg. sklearn, fancyimpute.
    method
        Imputation method in library, named.
    opts
        Additional options to assign to imputer."""
    def __init__(self, solver, method, opts={}):
        """
        """
        from fancyimpute import KNN, NuclearNormMinimization, SoftImpute, IterativeImputer, BiScaler
        from sklearn.impute import SimpleImputer
        import inspect
        imputers = {'fancyimpute':dict(KNN=KNN,MICE=IterativeImputer,BiScaler=BiScaler,Soft=SoftImpute),
                    'impyute':dict(),
                    'sklearn':dict(Mean=SimpleImputer(strategy='mean'),Zero=SimpleImputer(strategy='constant'))}
        try:
            from sklearn.impute import IterativeImputer
            f=IterativeImputer
            opts['n_nearest_features']=opts['k']
            opts={key: opts[key] for key in opts if key in inspect.getargspec(f.__init__)[0]}
            opts.update(dict(min_value=0.,sample_posterior=True,max_value=1.))
            imputers['sklearn']['MICE']=f(**opts)
        except:
            print('Please install development branch of iterative imputer from sklearn. Defaulting to mean imputation.')
            imputers['sklearn']['MICE']=imputers['sklearn']['Mean']
        try:
            if solver == 'fancyimpute':
                f=imputers[solver][method]
                opts={key: opts[key] for key in opts if key in inspect.getargspec(f.__init__)[0]}
                self.imputer=f(**opts)
            else:
                self.imputer = imputers[solver][method]
        except:
            print('{} {} not a valid combination.\nValid combinations:{}'.format(
                solver, method, '\n'.join('{}:{}'.format(solver,','.join(imputers[solver].keys())) for solver in imputers)))
            exit()

    def return_imputer(self):
        """Return initialized sklearn-like imputer."""
        return self.imputer

#### FUNCTIONS ####

def extract_pheno_beta_df_from_folder(folder):
    """Return phenotype and beta dataframes from specified folder with csv.

    Parameters
    ----------
    folder
        Input folder.
    """
    return pd.read_csv(folder+'/pheno.csv'), pd.read_csv(folder+'/beta.csv')

def extract_pheno_beta_df_from_pickle_dict(input_dict, disease=''):
    """Return phenotype and beta dataframes from specified dictionary storing MethylationArray python dictionary.

    Parameters
    ----------
    input_dict
        Python disctionary storing pheno/beta information.
    """
    if disease:
        return input_dict['pheno_{}'.format(disease)], input_dict['beta_{}'.format(disease)]
    else:
        return input_dict['pheno'], input_dict['beta']

def extract_pheno_beta_df_from_sql(conn, disease=''):
    """Return phenotype and beta dataframes from SQL tables storing MethylationArray info.

    Parameters
    ----------
    conn
        SQL connection.
    """
    if disease:
        return pd.read_sql('select * from {};'.format('pheno_{}'.format(disease)),conn), pd.read_sql('select * from {};'.format('beta_{}'.format(disease)),conn)
    else:
        return pd.read_sql('select * from {};'.format('pheno'),conn), pd.read_sql('select * from {};'.format('beta'),conn)
