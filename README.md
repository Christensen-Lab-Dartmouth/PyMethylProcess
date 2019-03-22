## PyMethylProcess

![Overview](https://raw.githubusercontent.com/Christensen-Lab-Dartmouth/PyMethylProcess/master/docs/yimages/pymethylprocess_overview.jpeg)

https://github.com/Christensen-Lab-Dartmouth/PyMethylProcess  

Help documentation: https://christensen-lab-dartmouth.github.io/PyMethylProcess/  

Alternatively, you can access the pdf: PyMethylProcess.pdf

What is it:
* Preprocess 450k and 850k methylation IDAT files in parallel using Minfi, ENmix, and meffil  
* Convenient and scalable implementation  
* Imputation and Feature Selection  
* Preparation for machine learning pipelines    

Why:
* Make DNAm accessible to python developers and more machine learning oriented researchers  
* Streamlined analysis makes processing easy  

Getting Started:  
* Installation:    
    * pip install pymethylprocess && pymethyl-install_r_dependencies  
    * docker pull joshualevy44/pymethylprocess  
    * Or see example scripts for usage, install from github
* Example Usage Scripts (in github repo): Located in ./example_scripts/  
* Help docs (in github repo): https://christensen-lab-dartmouth.github.io/PyMethylProcess/

*PyMethyProcess* is pending submission and review, and link to paper will be posted shortly.

![Download](https://raw.githubusercontent.com/Christensen-Lab-Dartmouth/PyMethylProcess/master/docs/yimages/pipeline-download.jpeg?token=ASyRZ6VozGl-fk1XitRFkeYYZJjArIe7ks5cnkqcwA%3D%3D)

![Format](https://raw.githubusercontent.com/Christensen-Lab-Dartmouth/PyMethylProcess/master/docs/yimages/pipeline-format.jpeg)

![PreProcess](https://raw.githubusercontent.com/Christensen-Lab-Dartmouth/PyMethylProcess/master/docs/yimages/pipeline-preprocess.jpeg)

![Visualize](https://raw.githubusercontent.com/Christensen-Lab-Dartmouth/PyMethylProcess/master/docs/yimages/pipeline-visualize.jpeg)

![Split](https://raw.githubusercontent.com/Christensen-Lab-Dartmouth/PyMethylProcess/master/docs/yimages/pipeline-train-test-split.jpeg)

Note: May need to prefix pip install with MACOSX_DEPLOYMENT_TARGET=10.9 CC=clang CXX=clang++ for Mac OS install
