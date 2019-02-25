import pandas as pd, numpy as np
import pickle, os

class MethylationArray:
    def __init__(self, pheno_df, beta_df, name=''):
        self.pheno=pheno_df
        self.beta=beta_df
        self.name=name

    def export(self, output_pickle):
        pass

    def write_csvs(self, output_dir):
        self.pheno.to_csv('{}/pheno.csv'.format(output_dir))
        self.beta.to_csv('{}/beta.csv'.format(output_dir))

    def write_pickle(self, output_pickle, disease=''):
        output_dict = {}
        if 0 and os.path.exists(output_pickle):
            output_dict = pickle.load(open(output_pickle,'rb'))
        output_dict['pheno'] = self.pheno #  if not disease else 'pheno_{}'.format(disease)
        output_dict['beta'] = self.beta #  if not disease else 'beta_{}'.format(disease)
        pickle.dump(output_dict, open(output_pickle,'wb'),protocol=4)

    def write_db(self, conn, disease=''):
        self.pheno.to_sql('pheno' if not disease else 'pheno_{}'.format(disease), con=conn, if_exists='replace')
        self.beta.to_sql('beta' if not disease else 'beta_{}'.format(disease), con=conn, if_exists='replace')

    def impute(self, imputer):
        self.beta = pd.DataFrame(imputer.fit_transform(self.beta),index=self.beta.index,columns=list(self.beta))

    def return_shape(self):
        return self.beta.shape

    def bin_column(self, col, n_bins):
        new_col_name='{}_{}'.format(col,'binned')
        self.pheno[new_col_name]=np.vectorize(lambda x: x.replace(' ',''))(pd.cut(self.pheno[col],n_bins).astype(str))
        return new_col_name

    def split_train_test(self, train_p=0.8, stratified=True, disease_only=False, key='disease', subtype_delimiter=',', val_p=0.): # https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html
        np.random.seed(42)
        if stratified:
            train_pheno=self.pheno.groupby(self.split_key('disease',subtype_delimiter) if disease_only else key,group_keys=False).apply(lambda x: x.sample(frac=train_p))
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
        np.random.seed(42)
        if categorical:
            pheno=self.pheno.groupby(key,group_keys=False).apply(lambda x: x.sample(frac=frac, n=min(x.shape[0],n_samples)))
        else:
            pheno=self.pheno.sample(frac=frac, n=min(n_samples,self.pheno.shape[0]))
        return MethylationArray(pheno,self.beta.loc[pheno.index.values,:],'subsampled')

    def split_key(self, key, subtype_delimiter):
        new_key = '{}_only'.format(key)
        self.pheno.loc[:,new_key] = self.pheno[key].map(lambda x: x.split(subtype_delimiter)[0])
        return new_key

    def remove_whitespace(self, key):
        self.pheno.loc[:,key]=self.pheno[key].map(lambda x: x.replace(' ',''))

    def split_by_subtype(self, disease_only, subtype_delimiter):
        for disease, pheno_df in self.pheno.groupby(self.split_key('disease',subtype_delimiter) if disease_only else 'disease'):
            new_disease_name = disease.replace(' ','')
            beta_df = self.beta.loc[pheno_df.index,:]
            yield MethylationArray(pheno_df,beta_df,new_disease_name)

    def feature_select(self, n_top_cpgs, feature_selection_method='mad', metric='correlation', nn=10):
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
        vals=self.pheno[outcome_cols].values
        if vals.shape[0] == 1:
            vals = vals[:,np.newaxis]
        remove_bool = ~np.isnan(vals).any(axis=1)
        self.pheno=self.pheno[remove_bool]
        self.beta=self.beta[remove_bool]

    def subset_index(self,index):
        return MethylationArray(self.pheno.loc[index,:],self.beta.loc[index,:])

    def return_idx(self):
        return np.array(list(self.pheno.index))

    def return_cpgs(self):
        return np.array(list(self.beta))

    def return_raw_beta_array(self):
        return self.beta.values

    def subset_cpgs(self,cpgs):
        return MethylationArray(self.pheno,self.beta.loc[:,cpgs])

    def categorical_breakdown(self, key):
        from collections import Counter
        print('\n'.join(['{}:{}'.format(k,v) for k,v in Counter(self.pheno[key].values).items()]))

    def overwrite_pheno_data(self, preprocess_sample_df):
        preprocess_sample_df=preprocess_sample_df.loc[self.pheno.index,:]
        for col in list(preprocess_sample_df):
            if col in list(self.pheno):
                self.pheno.loc[:,col]=preprocess_sample_df[col]
            else:
                self.pheno[col]=preprocess_sample_df[col]

    def merge_preprocess_sheet(self, preprocess_sample_df):
        self.pheno=self.pheno.merge(preprocess_sample_df,on=['Basename'],how='inner')
        if 'disease_x' in list(self.pheno):
            self.pheno = self.pheno.rename(columns={'disease_x':'disease'})
        self.pheno = self.pheno[[col for col in list(self.pheno) if not col.startswith('Unnamed:')]]
        self.pheno=self.pheno.set_index([np.vectorize(lambda x: x.split('/')[-1])(self.pheno['Basename'])],drop=False)

    def groupby(self, key):
        for name, df in self.pheno.groupby(key):
            yield (str(name), MethylationArray(pheno_df=self.pheno.loc[df.index],beta_df=self.beta.loc[df.index]))

    @classmethod
    def from_pickle(self, input_pickle):
        return MethylationArray(*extract_pheno_beta_df_from_pickle_dict(pickle.load(open(input_pickle,'rb'))))

class MethylationArrays:
    def __init__(self, list_methylation_arrays):
        self.methylation_arrays = list_methylation_arrays

    def __len__(self):
        return len(self.methylation_arrays)

    def combine(self, array_generator=None): # FIXME add sort based on samples
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
        for methyl_arr in self.methylation_arrays:
            methyl_arr.write_db(conn, methyl_arr.name)

    def write_pkls(self,pkl):
        for methyl_arr in self.methylation_arrays:
            methyl_arr.write_pickle(pkl, methyl_arr.name)

    def impute(self, imputer):
        for i in range(len(self.methylation_arrays)):
            self.methylation_arrays[i].impute(imputer)

class ImputerObject:
    def __init__(self, solver, method, opts={}):
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
            imputers['sklearn']['MICE']=imputers['sklearn']['SimpleImputer']
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
        return self.imputer

#### FUNCTIONS ####

def extract_pheno_beta_df_from_folder(folder):
    return pd.read_csv(folder+'/pheno.csv'), pd.read_csv(folder+'/beta.csv')

def extract_pheno_beta_df_from_pickle_dict(input_dict, disease=''):
    if disease:
        return input_dict['pheno_{}'.format(disease)], input_dict['beta_{}'.format(disease)]
    else:
        return input_dict['pheno'], input_dict['beta']

def extract_pheno_beta_df_from_sql(conn, disease=''):
    if disease:
        return pd.read_sql('select * from {};'.format('pheno_{}'.format(disease)),conn), pd.read_sql('select * from {};'.format('beta_{}'.format(disease)),conn)
    else:
        return pd.read_sql('select * from {};'.format('pheno'),conn), pd.read_sql('select * from {};'.format('beta'),conn)
