"""
meffil_functions.py
===================
Contains a few R functions that interact with meffil and minfi.
"""

import rpy2.robjects as robjects

def load_detection_p_values_beadnum(qc_list, n_cores):
    """Return list of detection p-value matrix and bead number matrix.

    Parameters
    ----------
    qc_list
        R list containing qc objects.
    n_cores
        Number of cores to use in computation.
    """
    pval_beadnum = robjects.r("""function(qc.list,mc.cores=1,
                                            max.bytes=2^30-1,
                                            verbose=F,
                                            ...) {
        qc.objects <- qc.list$qc.objects
        options(mc.cores=mc.cores)
        stopifnot(all(sapply(qc.objects, meffil:::is.qc.object)))

        featuresets <- sapply(qc.objects, function(qc.object) qc.object$featureset)
        featureset <- featuresets[1]

        if (is.list(featuresets)) ## backwards compatibility
            featureset <- featuresets <- "450k"

        if (any(featuresets != featureset))
            stop("Multiple feature sets were used to create these QC objects:",
                 paste(unique(featuresets), collapse=", "))

        feature.names <- meffil.get.features(featureset)$name

        if (!all(sapply(qc.objects, function(qc.object) meffil:::exists.rg(qc.object$basename))))
             stop("IDAT files are not accessible for all QC objects")

        ret.pvalue <- meffil:::mcsapply.safe(qc.objects, function(qc.object) {
            if (is.null(qc.object$featureset)) ## backwards compatibility
                qc.object$chip <- "450k"

            rg <- meffil:::read.rg(qc.object$basename, verbose=verbose)
            probes <- meffil.probe.info(qc.object$chip)
            pvalues <- meffil:::extract.detection.pvalues(rg, probes, verbose=verbose)
            unname(pvalues[feature.names])
        }, ..., max.bytes=max.bytes)

        ret.beadnum <- meffil:::mcsapply.safe(qc.objects, function(qc.object) {
            if (is.null(qc.object$featureset)) ## backwards compatibility
                qc.object$chip <- "450k"
            rg <- meffil:::read.rg(qc.object$basename, verbose=verbose)
            probes <- meffil.probe.info(qc.object$chip)
            beadnum <- meffil:::extract.beadnum(rg, probes, verbose=verbose)
            unname(beadnum[feature.names])
        }, ..., max.bytes=max.bytes)

        dimnames(ret.pvalue) <- list(feature.names, names(qc.objects))
        dimnames(ret.beadnum) <- list(feature.names, names(qc.objects))
        return(list(p.values=ret.pvalue, beadnum=ret.beadnum))
        }""")(qc_list, n_cores)
    return pval_beadnum

def set_missing(beta, pval_beadnum, detection_val=1e-6):
    """Set missing beta values to NA, taking into account detection values and bead number thesholds.

    Parameters
    ----------
    pval_beadnum
        Detection pvalues and number of beads per cpg/samples
    detection_val
        If threshold to set site to missingness based on p-value detection.
    """
    beta = robjects.r("""function (beta, pval.beadnum, detection.p=1e-6){
        p.values <- pval.beadnum$p.values[rownames(beta),colnames(beta)]
        beadnum <- pval.beadnum$beadnum[rownames(beta),colnames(beta)]
        beta[((p.values >= detection.p)+(beadnum<3))>0]<-NA
        return(beta)
        }""")(beta, pval_beadnum, detection_val)
    return beta

def remove_sex(beta, array_type='450k'):
    """Remove non-autosomal cpgs from beta matrix.

    Parameters
    ----------
    array_type
        450k/850k array?
    """
    beta = robjects.r("""function (beta,array.type){
        featureset<-array.type
        autosomal.sites <- meffil.get.autosomal.sites(featureset)
        autosomal.sites <- intersect(autosomal.sites, rownames(beta))
        beta <- beta[autosomal.sites,]
        return(beta)
        }""")(beta,array_type)
    return beta

def r_autosomal_cpgs(array_type='450k'):
    """Return list of autosomal cpg probes per platform.

    Parameters
    ----------
    array_type
        450k/850k array?
    """
    robjects.r('library(meffil)')
    cpgs = robjects.r("""meffil.get.autosomal.sites('{}')""".format(array_type))
    return cpgs

def r_snp_cpgs(array_type='450k'):
    """Return list of SNP cpg probes per platform.

    Parameters
    ----------
    array_type
        450k/850k array?
    """
    robjects.r('library(meffil)')
    cpgs = robjects.r("""meffil.snp.names('{}')""".format(array_type))
    return cpgs

def est_cell_counts_meffil(qc_list, cell_type_reference):
    """Given QCObject list R object, estimate cell counts using reference approach via meffil.

    Parameters
    ----------
    qc_list
        R list containing qc objects.
    cell_type_reference
        Reference blood/tissue set."""
    cell_count_estimates = robjects.r("""function (qc.list, cell.type.reference) {
        qc.objects <- qc.list$qc.objects
        cc<-t(sapply(qc.objects, function(obj) meffil.estimate.cell.counts(obj,cell.type.reference)))
        cc<-data.frame(IID=row.names(cc),cc)
        return(cc)
        }""")(qc_list,cell_type_reference)
    return cell_count_estimates

def est_cell_counts_minfi(rgset):
    """Given RGSet object, estimate cell counts using reference approach via minfi.

    Parameters
    ----------
    rgset
        RGSet object stored in python via rpy2"""
    robjects.r('library(FlowSorted.Blood.450k)')
    cell_count_estimates = robjects.r("""function (RGset) {
        cellCounts <- as.table(estimateCellCounts(RGset))
        return(cellCounts)
        }""")(rgset)
    return cell_count_estimates

def est_cell_counts_IDOL(rgset,library):
    """Given RGSet object, estimate cell counts for 450k/850k using reference approach via IDOL library.

    Parameters
    ----------
    rgset
        RGSet object stored in python via rpy2
    library
        What type of CpG library to use."""
    robjects.r('library(FlowSorted.Blood.EPIC)')
    cell_count_estimates = robjects.r("""function (RGset) as.table(estimateCellCounts2(RGset,IDOLOptimizedCpGs={})$counts)""".format(library))(rgset)
    return cell_count_estimates


def bmiq_mc(beta, nCores, nfit):
	# Credits to Lucas A. Salas
	from rpy2.robjects.packages import importr
	enmix=importr('ENmix')
	meffil=importr('meffil')
	array_types={'IlluminaHumanMethylation450k':'450k','IlluminaHumanMethylationEPIC':'850k'}
	try:
		beta=remove_sex(beta, array_type=array_types['IlluminaHumanMethylation450k'])
	except:
		beta = remove_sex(beta, array_type=array_types['IlluminaHumanMethylationEPIC'])
	imputation_args=dict(k = 10,rowmax=0.5,colmax=0.8,maxp=1500)
	beta=enmix.rm_outlier(beta,impute=True,rmcr=True, **imputation_args)
	beta=robjects.r("""bmiq.mc2<-function (mdat, nCores = 1, ...)
{
 if (!is(mdat, "matrix")) {
   stop("object needs to be of class 'matrix'")
 }
 if (nCores > detectCores()) {
   nCores <- detectCores()
   cat("Only ", nCores, " cores are available in your computer,",
	   "argument nCores was reset to nCores=", nCores, "\n")
 }
 library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
 anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
 #mdat<-betasronnerblad
 anno<-anno[rownames(mdat),]
 cat("Analysis is running, please wait...!", "\n")
 beta.b <- mdat
 rm(mdat)
 ncpgs<-nrow(beta.b)
 beta.b[beta.b <= 0] <- 1e-06
 design.v <- as.vector(anno$Type)
 design.v[design.v == "I"] = 1
 design.v[design.v == "II"] = 2
 design.v <- as.numeric(design.v)
 coln = colnames(beta.b)
 N = ceiling(ncol(beta.b)/(nCores * 10))
 parts = rep(1:N, each = ceiling(ncol(beta.b)/N))[1:ncol(beta.b)]
 c1 <- makeCluster(nCores)
 registerDoParallel(c1)
 library(wateRmelon)
 #writeLines(c(""), "log.txt")
 for (i in 1:N) {
   id = which(parts == i)
   beta.b1 = beta.b[, id]
   #print(id)
   #print(nrow(beta.b1),ncol(beta.b1))
   beta.b1 <- foreach(s = 1:ncol(beta.b1), .combine = cbind,
					  .export = c("BMIQ")) %dopar% {
						s = s
						#sink("log.txt", append=TRUE)
						#print(s)
						#print(beta.b1[, s])
						out <- tryCatch(BMIQ(beta.b1[, s], design.v = design.v, plots = FALSE)$nbeta,error=function(e){rep(-1,ncpgs)})
						out
					  }
   beta.b[, id] = beta.b1
 }
 stopCluster(c1)
 #print(as(beta.b,'matrix'))
 if (is.matrix(beta.b)) {
   if (sum(is.na(beta.b)) > 0) {
	 stop("BMIQ estimates \n    encountered error, try to run it again")
   }
 }
 else {
   stop("BMIQ estimates encountered error, try to run it again")
 }
 colnames(beta.b) <- coln
 return(beta.b)
}""")(beta, nCores, nfit)
	return beta
