from rpy2.robjects.packages import importr


def importr_(package, r_lib_file='rlibs.txt'):
    try:
        return importr(package)
    except:
        with open('rlibs.txt') as f:
            rlibs=f.read().splitlines()
        for rlib in rlibs:
            try:
                return importr(package, lib_loc=rlib)
            except:
                pass
        return None
