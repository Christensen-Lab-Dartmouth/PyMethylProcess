FROM rpy2/rpy2:latest

RUN apt-get update -y

RUN apt-get install -y vim nano

RUN apt-get remove curl && apt-get install -y curl

RUN apt-get remove libssl-dev && apt-get install -y libssl-dev

RUN apt-get install -y libgtk2.0-dev xvfb xauth xfonts-base libcairo2-dev libcurl4-openssl-dev # libxaw xlib

RUN apt-get install -y libx11-dev openssl libxt-dev #libxaw xlib

RUN pip install --upgrade rpy2

RUN pip install pymethylprocess

RUN pymethyl-install_r_dependencies
