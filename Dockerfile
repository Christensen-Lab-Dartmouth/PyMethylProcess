FROM rpy2/rpy2:latest

RUN apt-get update -y

RUN apt-get install -y vim nano

RUN apt-get remove curl && apt-get install -y curl

RUN apt-get remove libssl-dev && apt-get install -y libssl-dev

RUN apt-get install -y libgtk2.0-dev xvfb xauth xfonts-base libcairo2-dev libcurl4-openssl-dev # libxaw xlib

RUN apt-get install -y libx11-dev openssl libxt-dev #libxaw xlib

RUN pip install --upgrade rpy2

RUN pip install --upgrade pip

RUN pip install numpy kneed Cython pathos nevergrad

RUN pip install umap-learn>=0.3.7 plotly>=3.4.2 fancyimpute>=0.4.2 pandas>=0.23.4 scikit-learn>=0.20.1

RUN pip install shap matplotlib seaborn mlxtend click==6.7

RUN pip install git+https://github.com/bodono/scs-python.git@bb45c69ce57b1fbb5ab23e02b30549a7e0b801e3

RUN pip install git+https://github.com/jlevy44/hypopt.git@27aefef62483174736bd6d5a1b3983dbaf4184dc

RUN pip install -i https://test.pypi.org/simple/ pymethylprocess --no-deps

RUN pymethyl-install_r_dependencies
