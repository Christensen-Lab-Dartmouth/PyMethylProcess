FROM joshualevy44/pymethylprocess:0.1.3

RUN mkdir /pymethyl_repo/

COPY . /pymethyl_repo/

RUN apt-get -y install python3-setuptools

RUN ls /pymethyl_repo

RUN cd /pymethyl_repo/ && python3 setup.py sdist bdist_wheel && pip3 install dist/pymethylprocess-0.1.4.tar.gz --no-deps && cd -

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
