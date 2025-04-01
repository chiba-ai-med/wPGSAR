#' Weighted Parametric Gene Set Analysis (wPGSA)
#'
#' Weighted Parametric Gene Set Analysis (wPGSA) to infer the activity of transcription factors by evaluating the expression variance of their target genes using weighted t-statistics. When group labels are provided, the package can also perform statistical tests to compare transcription factor activity between groups.
#' @param network_data Output object of retrieveNetworkData()
#' @param expression_data Gene expression matrix specified by user
#' @param group Group label vector to perform statistical tests to compare TF activity between groups. (optional)
#' @returns A list containing the result of wPGSA:
#' \describe{
#'   \item{tstat}{TF activity scores (t-statistics) computed from weighted expression and background.}
#'   \item{tstat_se}{Standard error used in tstat calculation.}
#'   \item{group_mean}{Group-wise averages of t-statistics, estimated by linear regression.}
#'   \item{group_se}{Standard errors of group means.}
#'   \item{group_t}{t-values of group means (group_mean / group_se).}
#'   \item{pval}{P-values from linear regression testing whether group means differ significantly from 0.}
#'   \item{fdr}{Adjusted p-values (FDR) using Benjamini-Hochberg correction.}
#' }
#' @examples
#' library("wPGSAR")
#' data(GSE143371_count)
#' data(GSE143371_group)
#' network_data <- retrieveNetworkData()
#' res1 <- wPGSA(network_data, GSE143371_count)
#' res2 <- wPGSA(network_data, GSE143371_count, GSE143371_group)
#' @importFrom Matrix sparseMatrix colMeans rowSums
#' @importFrom stats coef lm p.adjust var setNames
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom data.table as.data.table
#' @importFrom methods is
#' @export
wPGSA <- function(network_data, expression_data, group=NULL){
	# Argument Check
	commonGeneName <- .checkwPGSA(network_data, expression_data)
	# Background statistics
	cat("(1/4) Background statistics are calculated.\n")
	global_mean <- .global_mean(expression_data)
	global_var <- apply(expression_data, 2, var)
	# Filtering
	res.filter <- .filterwPGSA(
		network_data, expression_data, commonGeneName)
	# TF-specific statistics by Sparse Matrix Product
	cat("(2/4) TF-specific statistics are calculated.\n")
	network_sparse <- .sparseMatrix(
		res.filter$network_data, res.filter$TFs, commonGeneName)
	total_sum <- rowSums(network_sparse)
	target_mean <- .targetMean(network_sparse,
			res.filter$expression_data, total_sum)
	target_var <- .targetVar(network_sparse,
			res.filter$expression_data, total_sum, target_mean)
	# t-statistics
	cat("(3/4) TF enrichment scores are calculated.\n")
	tstat_se <- .tstat_se(global_var, target_var, total_sum)
	tstat <- .tstat(global_mean, target_mean, tstat_se)
	# Statistical tests
	res.tests <- .tests(tstat, group)
	# Output
	list(tstat=tstat, tstat_se=tstat_se,
		group_mean=res.tests$group_mean,
		group_se=res.tests$group_se, group_t=res.tests$group_t,
		pval=res.tests$pval, fdr=res.tests$fdr)
}

.checkwPGSA <- function(network_data, expression_data){
	commonGeneName <- intersect(
		unique(network_data$GeneName),
		rownames(expression_data))
	if(length(commonGeneName) == 0){
		msg1 <- "There is no common gene name "
		msg2 <- "between network_data and expression_data!"
		stop(paste0(msg1, msg2))
	}
	commonGeneName
}

.global_mean <- function(expression_data){
	if("dgCMatrix" %in% is(expression_data)){
		Matrix::colMeans(expression_data)
	}else{
		base::colMeans(expression_data)
	}
}

.filterwPGSA <- function(network_data, expression_data, commonGeneName){
	network_data <- network_data[
		network_data$GeneName %in% commonGeneName, ]
	TFs <- unique(network_data$TF_name)
	expression_data <- expression_data[
		commonGeneName, , drop = FALSE]	
	list(network_data=network_data,
		TFs=TFs,
		expression_data=expression_data)
}

.sparseMatrix <- function(network_data, TFs, commonGeneName) {
  dt <- as.data.table(network_data)
  tf_index <- setNames(seq_along(TFs), TFs)
  gene_index <- setNames(seq_along(commonGeneName), commonGeneName)
  dt[, tf_idx := tf_index[TF_name]]
  dt[, gene_idx := gene_index[GeneName]]
  sparseMatrix(
    i = dt$tf_idx,
    j = dt$gene_idx,
    x = dt$weight,
    dims = c(length(TFs), length(commonGeneName)),
    dimnames = list(TFs, commonGeneName)
  )
}

.targetMean <- function(network_sparse, expression_data, total_sum){
	(network_sparse %*% expression_data) / total_sum
}

.targetVar <- function(network_sparse, expression_data, total_sum, target_mean){
	(network_sparse %*% expression_data^2) / total_sum - target_mean^2
}

.tstat_se <- function(global_var, target_var, total_sum){
	sqrt(Matrix::t(global_var + Matrix::t(target_var)) / total_sum)
}

.tstat <- function(global_mean, target_mean, tstat_se){
	Matrix::t(Matrix::t(target_mean) - global_mean) / tstat_se
}

.extract_col <- function(tmp, i, group_levels){
	out <- do.call(rbind, lapply(tmp, function(x){x[, i]}))
	colnames(out) <- group_levels
	out
}

.apply_pb <- function(X, MARGIN, FUN, ...)
{
  env <- environment()
  pb_Total <- sum(dim(X)[MARGIN])
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total,
                       style = 3)

  wrapper <- function(...)
  {
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir= env)
    setTxtProgressBar(get("pb", envir= env),
                           curVal +1)
    FUN(...)
  }
  res <- apply(X, MARGIN, wrapper, ...)
  close(pb)
  res
}

.lmTestCore <- function(x, group){
  df <- data.frame(val = x, group = group)
  fit <- lm(val ~ 0 + group, data = df)
  coef(summary(fit))
}

.lmTest <- function(tstat, group){
  .apply_pb(tstat, 1, function(x) .lmTestCore(x, group), simplify=FALSE)
}

.tests <- function(tstat, group){
	if(!is.null(group)){
		cat("(4/4) Statistical tests between groups are performed.\n")
		group <- as.factor(group)
		group_levels <- levels(group)
		tmp <- .lmTest(tstat, group)
		group_mean <- .extract_col(tmp, 1, group_levels)
		group_se <- .extract_col(tmp, 2, group_levels)
		group_t <- .extract_col(tmp, 3, group_levels)
		pval <- .extract_col(tmp, 4, group_levels)
		fdr <- apply(pval, 2, p.adjust)
	}else{
		cat("(4/4) Statistical tests between groups are skipped.\n")
		group_mean <- NULL
		group_se <- NULL
		group_t <- NULL
		pval <- NULL
		fdr <- NULL
	}
	list(group_mean=group_mean, group_se=group_se,
		group_t=group_t, pval=pval, fdr=fdr)
}