from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import click
import os, subprocess


CONTEXT_SETTINGS = dict(help_option_names=['-h','--help'], max_content_width=90)

@click.group(context_settings= CONTEXT_SETTINGS)
@click.version_option(version='0.1')
def install():
    pass

class PackageInstaller:
    def __init__(self):
        self.lib_xml_location=os.popen('type xml2-config').read().split()[-1]
        self.lib_xml_location=robjects.r('function (xml.config) {Sys.setenv(XML_CONFIG=xml.config)}')(self.lib_xml_location)

    def install_bioconductor(self):
        base = importr('base')
        base.source("http://www.bioconductor.org/biocLite.R")

    def install_tcga_biolinks(self):
        biocinstaller = importr("BiocInstaller")
        biocinstaller.biocLite("TCGAbiolinks")

    def install_minfi_others(self):
        biocinstaller = importr("BiocInstaller")
        biocinstaller.biocLite(robjects.vectors.StrVector(["minfi","ENmix",
                                "minfiData","sva","GEOquery","geneplotter"]))

    def install_custom(self, custom, manager):
        if not manager:
            biocinstaller = importr("BiocInstaller")
            biocinstaller.biocLite(robjects.vectors.StrVector(custom),suppressUpdates=True)
        else:
            biocinstaller = importr("BiocManager")
            for c in custom:
                if '=' in c:
                    pkg,version= tuple(c.split('='))
                    biocinstaller.install(pkg,ask=False,version=version)
                else:
                    biocinstaller.install(c,ask=False)

    def install_devtools(self):
        subprocess.call('conda install -y -c r r-cairo=1.5_9 r-devtools=1.13.6',shell=True)
        robjects.r('install.packages')('devtools')

    def install_r_packages(self, custom):
        robjects.r["options"](repos=robjects.r('structure(c(CRAN="http://cran.wustl.edu/"))'))
        robjects.r('install.packages')(robjects.vectors.StrVector(custom))

    def install_meffil(self, git=False):
        if git:
            remotes=importr('remotes')
            remotes.install_github('perishky/meffil')
        else:
            subprocess.call("wget https://github.com/perishky/meffil/archive/master.zip && unzip master.zip && mv meffil-master meffil && R CMD INSTALL meffil",shell=True)

@install.command()
def change_gcc_path():
    """Change GCC and G++ paths if don't have version 7.2.0. [Experimental]"""
    bin_path = os.path.join(os.popen('conda list | grep "packages in environment at" | awk "{print $6}"').read().split()[-1].replace(':',''),'bin')
    subprocess.call('export CC={}'.format(os.path.join(bin_path,'x86_64-conda_cos6-linux-gnu-gcc')),shell=True)
    subprocess.call('export CXX={}'.format(os.path.join(bin_path,'x86_64-conda_cos6-linux-gnu-g++')),shell=True)

## Install ##
@install.command()
def install_bioconductor():
    """Installs bioconductor."""
    installer = PackageInstaller()
    installer.install_bioconductor()

@install.command()
@click.option('-p', '--package', multiple=True, default=['ENmix'], help='Custom packages.', type=click.Path(exists=False), show_default=True)
@click.option('-m', '--manager', is_flag=True, help='Use BiocManager (recommended).')
def install_custom(package,manager):
    """Installs bioconductor packages."""
    installer = PackageInstaller()
    installer.install_custom(package,manager)

@install.command()
@click.option('-p', '--package', multiple=True, default=[''], help='Custom packages.', type=click.Path(exists=False), show_default=True)
def install_r_packages(package):
    """Installs r packages."""
    installer = PackageInstaller()
    installer.install_r_packages(package)

@install.command()
def install_minfi_others():
    """Installs minfi and other dependencies."""
    installer = PackageInstaller()
    installer.install_minfi_others()

@install.command()
def install_tcga_biolinks():
    """Installs tcga biolinks."""
    installer = PackageInstaller()
    installer.install_tcga_biolinks()

@install.command()
def install_meffil():
    """Installs meffil (update!)."""
    installer = PackageInstaller()
    installer.install_meffil()

@install.command()
def install_some_deps():
    """Installs bioconductor, minfi, enmix, tcga biolinks, and meffil."""
    installer = PackageInstaller()
    installer.install_bioconductor()
    installer.install_minfi_others()
    installer.install_tcga_biolinks()
    installer.install_meffil()

# install.packages("png", "/home/user/anaconda3/lib/R/library")



if __name__ == '__main__':
    install()
