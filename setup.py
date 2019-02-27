from setuptools import setup
setup(name='pymethylprocess',
      version='0.1',
      description='What the module does',
      url='https://github.com/Christensen-Lab-Dartmouth/PyMethylProcess',
      author='Joshua Levy',
      author_email='joshualevy44@berkeley.edu',
      license='MIT',
      scripts=['bin/run_parallel','bin/pymethyl-install_r_dependencies'],
      entry_points={
            'console_scripts':['pymethyl-basic-install=pymethylprocess.basic_installer:main',
                               'pymethyl-install=pymethylprocess.installer:install',
                               'pymethyl-visualize=pymethylprocess.visualizations:visualize',
                               'pymethyl-preprocess=pymethylprocess.preprocess:preprocess',
                               'pymethyl-utils=pymethylprocess.utils:util']
      },
      packages=['pymethylprocess'],
      install_requires=['rpy2',
                        'numpy',
                        'scs @ git+https://github.com/bodono/scs-python.git@bb45c69ce57b1fbb5ab23e02b30549a7e0b801e3',
                        'kneed',
                        'Cython',
                        'pathos',
                        'nevergrad',
                        'umap-learn>=0.3.7',
                        'plotly>=3.4.2',
                        'fancyimpute>=0.4.2',
                        'pandas>=0.23.4',
                        'scikit-learn>=0.20.1',
                        'shap',
                        'matplotlib',
                        'seaborn',
                        'mlxtend',
                        'click==6.7',
                        'scs @ git+https://github.com/jlevy44/hypopt.git@27aefef62483174736bd6d5a1b3983dbaf4184dc'])
