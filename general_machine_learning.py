from MethylationDataTypes import MethylationArray
import pickle, numpy as np

class MachineLearning:
    def __init__(self, model, options):
        self.model=model(**options)

    def fit(self, train_methyl_array, outcome_cols=None):
        if self.outcome_cols != None:
            self.model.fit(train_methyl_array.beta,train_methyl_array.pheno[outcome_cols])
        else:
            self.model.fit(train_methyl_array.beta)
        return self.model

    def transform(self, test_methyl_array):
        self.results=self.model.transform(test_methyl_array)
        return self.results

    def fit_transform(self, train_methyl_array, outcome_cols=None):
        self.results=self.fit(train_methyl_array, outcome_cols=None).transform(train_methyl_array)
        return self.results

    def predict(self, test_methyl_array):
        self.results=self.model.predict(test_methyl_array)
        return self.results

    def fit_predict(self, train_methyl_array, outcome_cols=None):
        self.results=self.fit(train_methyl_array, outcome_cols).predict(train_methyl_array)
        return self.results

    def store_results(self, output_pkl, results_dict={}):
        if not results_dict:
            results_dict=dict(results=self.results)
        pickle.dump(results_dict,open(results_dict,'wb'))

    def assign_results_to_pheno_col(self, methyl_array, new_col, output_pkl):
        methyl_array.pheno[new_col]=self.results
        methyl_array.write_pickle(output_pkl)

    def transform_results_to_beta(self, methyl_array, output_pkl):
        methyl_array.beta = pd.DataFrame(self.results,index=self.beta.index)
        methyl_array.write_pickle(output_pkl)

    def return_outcome_metric(self, methyl_array, outcome_cols, metric, run_bootstrap=False):
        y_true = methyl_array.pheno[outcome_cols]
        y_pred = self.results
        if not bootstrap:
            return metric(y_true,y_pred)
        else:
            from sklearn.utils import resample
            boot_results=np.array([metric(*resample(y_true,y_pred,random_state=123)) for i in range(n_bootstrap)])
            original = metric(y_true,y_pred)
            std_err = np.std(boot_results)
            boot_results=np.sort(boot_results)
            ci=0.95
            bound=(1-ci)/2.
            # BORROWED FROM MLXTEND
            def quantile(x, q):
                rank = round(q * x.shape[0]) - 1
                if rank >= x.shape[0]:
                    rank = x.shape[0]
                elif rank <= 0:
                    rank = 0
                rank = int(round(rank))
                return x[rank]
            high_ci = quantile(boot_results, q=(ci + bound))
            low_ci = quantile(boot_results, q=bound)
            return original, std_err, (low_ci,high_ci)
