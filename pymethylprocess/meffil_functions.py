import rpy2.robjects as robjects

def load_detection_p_values_beadnum(qc_list, n_cores):
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
    beta = robjects.r("""function (beta, pval.beadnum, detection.p=1e-6){
        p.values <- pval.beadnum$p.values[rownames(beta),colnames(beta)]
        beadnum <- pval.beadnum$beadnum[rownames(beta),colnames(beta)]
        beta[((p.values >= detection.p)+(beadnum<3))>0]<-NA
        return(beta)
        }""")(beta, pval_beadnum, detection_val)
    return beta

def remove_sex(beta, array_type='450k'):
    beta = robjects.r("""function (beta,array.type){
        featureset<-array.type
        autosomal.sites <- meffil.get.autosomal.sites(featureset)
        autosomal.sites <- intersect(autosomal.sites, rownames(norm.beta))
        norm.beta <- norm.beta[autosomal.sites,]
        return(beta)
        }""")(beta,array_type)
    return beta

def r_autosomal_cpgs(array_type='450k'):
    robjects.r('library(meffil)')
    cpgs = robjects.r("""meffil.get.autosomal.sites('{}')""".format(array_type))
    return cpgs

def r_snp_cpgs(array_type='450k'):
    robjects.r('library(meffil)')
    cpgs = robjects.r("""meffil.snp.names('{}')""".format(array_type))
    return cpgs

def est_cell_counts_meffil(qc_list, cell_type_reference):
    cell_count_estimates = robjects.r("""function (qc.list, cell.type.reference) {
        qc.objects <- qc.list$qc.objects
        cc<-t(sapply(qc.objects, function(obj) meffil.estimate.cell.counts(obj,cell.type.reference)))
        cc<-data.frame(IID=row.names(cc),cc)
        return(cc)
        }""")(qc_list,cell_type_reference)
    return cell_count_estimates

def est_cell_counts_minfi(rgset):
    robjects.r('library(FlowSorted.Blood.450k)')
    cell_count_estimates = robjects.r("""function (RGset) {
        cellCounts <- as.table(estimateCellCounts(RGset))
        return(cellCounts)
        }""")(rgset)
    return cell_count_estimates

def est_cell_counts_IDOL(rgset,library):
    robjects.r('library(FlowSorted.Blood.EPIC)')
    cell_count_estimates = robjects.r("""function (RGset) as.table(estimateCellCounts2(RGset,IDOLOptimizedCpGs={})$counts)""".format(library))(rgset)
    return cell_count_estimates
