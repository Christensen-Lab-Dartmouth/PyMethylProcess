import rpy2.robjects as robjects

def load_detection_p_values_beadnum(qc_list, n_cores):
    pval_beadnum = robjects.r("""function(qc.list,mc.cores=1,
                                            max.bytes=2^30-1,
                                            verbose=F,
                                            ...) {
        #library(meffil)
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
    beta = robjects.r("""function (beta, pval.beadnum, detection.p=1e-6){
        p.values <- pval.beadnum$p.values[rownames(beta),colnames(beta)]
        beadnum <- pval.beadnum$beadnum[rownames(beta),colnames(beta)]
        beta[((p.values >= detection.p)+(beadnum<3))>0]<-NA
        return(beta)
        }""")(beta, pval_beadnum, detection_val)
    return beta
