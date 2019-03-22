"""
general_machine_learning.py
===========================
Contains a machine learning class to perform scikit-learn like operations, along with held-out hyperparameter grid search.
"""
from pymethylprocess.MethylationDataTypes import MethylationArray
import pickle, numpy as np
from hypopt import GridSearch
from sklearn.preprocessing import LabelEncoder

class MachineLearning:
    """Machine learning class to run sklearn-like pipeline on MethylationArray data.
    Initialize object with scikit-learn model, and optionally supply a hyperparameter search grid.

    model
        Scikit-learn-like model, classification, regression, dimensionality reduction, clustering etc.
    options
        Options to supply model in form of dictionary.
    grid
        Alternatively, supply search grid to search for bets hyperparameters.
    labelencode
        T/F encode string labels.
    n_eval
        Number of evaluations for randomized grid search, if set to 0, perform exhaustive grid search
    """
    def __init__(self, model, options, grid={}, labelencode=False, n_eval=0):
        if grid:
            self.model = GridSearch(model = model(), param_grid=grid, num_random_search=None if not n_eval else n_eval)
            self.param_grid_exists=True
            self.grid=grid
        else:
            self.model=model(**options)
            self.param_grid_exists=False
        if labelencode:
            self.encoder=LabelEncoder()
        else:
            self.encoder=None

    def fit(self, train_methyl_array, val_methyl_array=None, outcome_cols=None):
        """Fit data to model.

        Parameters
        ----------
        train_methyl_array
            Training MethylationArray.
        val_methyl_array
            Validation MethylationArray. Can set to None.
        outcome_cols
            Set to none if not needed, but phenotype column to train on, can be multiple.
        """
        if outcome_cols != None:
            if self.encoder != None:
                self.encoder.fit(train_methyl_array.pheno[outcome_cols])
            if self.param_grid_exists:
                self.model.fit(train_methyl_array.beta,self.encoder.transform(train_methyl_array.pheno[outcome_cols]) if self.encoder != None else train_methyl_array.pheno[outcome_cols], val_methyl_array.beta,self.encoder.transform(val_methyl_array.pheno[outcome_cols]) if self.encoder != None else val_methyl_array.pheno[outcome_cols], scoring='accuracy' if self.encoder!=None else 'r2')
            else:
                self.model.fit(train_methyl_array.beta,self.encoder.transform(train_methyl_array.pheno[outcome_cols]) if self.encoder != None else train_methyl_array.pheno[outcome_cols])
        else:
            self.model.fit(train_methyl_array.beta)
        return self.model

    def transform(self, test_methyl_array):
        """Transform test methylation array.

        Parameters
        ----------
        test_methyl_array
            Testing MethylationArray.
        """
        self.results=self.model.transform(test_methyl_array.beta)
        return self.results

    def fit_transform(self, train_methyl_array, outcome_cols=None):
        """Fit and transform to training data.

        Parameters
        ----------
        train_methyl_array
            Training MethylationArray.
        outcome_cols
            Set to none if not needed, but phenotype column to train on, can be multiple.
        """
        self.results=self.fit(train_methyl_array, outcome_cols=None).transform(train_methyl_array)
        return self.results

    def predict(self, test_methyl_array):
        """Make new predictions on test methylation array.

        Parameters
        ----------
        test_methyl_array
            Testing MethylationArray.
        """
        self.results=self.model.predict(test_methyl_array.beta)
        if self.encoder != None:
            self.results=self.encoder.inverse_transform(self.results)
        return self.results

    def fit_predict(self, train_methyl_array, outcome_cols=None):
        """Fit and predict training data.

        Parameters
        ----------
        train_methyl_array
            Training MethylationArray.
        outcome_cols
            Set to none if not needed, but phenotype column to train on, can be multiple.
        """
        self.results=self.fit(train_methyl_array, outcome_cols).predict(train_methyl_array)
        return self.results

    def store_results(self, output_pkl, results_dict={}):
        """Store results in pickle file.

        Parameters
        ----------
        output_pkl
            Output pickle to dump results to.
        results_dict
            Supply own results dict to be dumped.
        """
        if not results_dict:
            results_dict=dict(results=self.results)
        pickle.dump(results_dict,open(results_dict,'wb'))

    def assign_results_to_pheno_col(self, methyl_array, new_col, output_pkl):
        """Assign results to new phenotype column.

        Parameters
        ----------
        methyl_array
            MethylationArray.
        new_col
            New column name.
        output_pkl
            Output pickle to dump MethylationArray to.
        """
        methyl_array.pheno[new_col]=self.results
        methyl_array.write_pickle(output_pkl)

    def transform_results_to_beta(self, methyl_array, output_pkl):
        """Transform beta matrix into reduced beta matrix and store.

        Parameters
        ----------
        methyl_array
            MethylationArray.
        output_pkl
            Output pickle to dump MethylationArray to.
        """
        methyl_array.beta = pd.DataFrame(self.results,index=self.beta.index)
        methyl_array.write_pickle(output_pkl)

    def return_outcome_metric(self, methyl_array, outcome_cols, metric, run_bootstrap=False):
        """Supply metric to evaluate results.

        Parameters
        ----------
        methyl_array
            MethylationArray to evaluate.
        outcome_cols
            Outcome phenotype columns.
        metric
            Sklearn evaluation metric.
        run_bootstrap
            Make 95% CI from 1k bootstraps.
        """
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
