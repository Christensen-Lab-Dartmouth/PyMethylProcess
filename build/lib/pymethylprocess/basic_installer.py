import argparse
import subprocess

def install_r_packages(packages,bioconductor=False, github=False):
    if github:
        for package in packages:
            subprocess.call('Rscript -e "library(\'devtools\'); {}"'.format("install_github('{}')".format(package)),shell=True)
    else:
        if bioconductor:
            for package in packages:
                subprocess.call('Rscript -e "library(\'BiocManager\'); {}"'.format("install('{}',ask=F)".format(package)),shell=True)
        else:
            packages = 'c({})'.format(','.join(["'{}'".format(package) for package in packages]))
            subprocess.call('Rscript -e "options(repos=structure(c(CRAN=\'http://cran.wustl.edu/\'))); install.packages({})"'.format(packages),shell=True)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-g', '--github', action='store_true', help='Install from github.')
    p.add_argument('-b', '--bioconductor', action='store_true', help='Install from bioconductor.')
    p.add_argument('-p','--packages', nargs='+', type=str, help='List of packages.', required=True)
    args=p.parse_args()
    github=args.github
    bioconductor=args.bioconductor
    packages = args.packages
    install_r_packages(packages,bioconductor, github)


if __name__=='__main__':
    main()
