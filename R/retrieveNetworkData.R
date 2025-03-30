#' Retrieval of TF - target gene table and TF metadata from ChIP-Atlas
#'
#' Two tab-delimited files are retrieved from ChIP-Atlas.
#' To check the local cache path, type the code below:
#' tools::R_user_dir("BiocFileCache", which = "cache")
#' @param assembly Genome Assembly (hg19 or mm10)
#' @param out.dir Output Directory (Default: tempdir())
#' @returns A corresponding table between TF (SRX ID) and target gene (Gene Symbol) and TF metadata (SRX_ID, TF, celltype1, celltype2)
#' @examples
#' library("wPGSAR")
#' network_data <- retrieveNetworkData()
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcrpath bfcneedsupdate bfcinfo
#' @importFrom tools R_user_dir
#' @importFrom dplyr left_join
#' @importFrom utils read.delim
#' @export
retrieveNetworkData <- function(assembly=c("hg19", "mm10"),
	out.dir=tempdir()){
	# Argument Check
	assembly <- match.arg(assembly)
	if(assembly == "hg19"){
		# URL
		url1 <- "https://chip-atlas.dbcls.jp/data/hg38/wpgsa/wpgsa_GRN.hg38.5.tsv"
		url2 <- "https://chip-atlas.dbcls.jp/data/hg38/wpgsa/wpgsa_hash.hg38.5.tsv"
	}else{
		stop("Mouse Data (mm10) is not available for now...")
	}
	# File Cache
	bfc <- BiocFileCache()
	# Get Cache ID
	res1 <- bfcquery(bfc, url1, field = "rname")
	res2 <- bfcquery(bfc, url2, field = "rname")
	# Check whether exist
	rid1 <- res1$rid
	rid2 <- res2$rid
	is_new1 <- length(rid1) == 0
	is_new2 <- length(rid2) == 0
	.checkCache(bfc, is_new1, rid1)
	.checkCache(bfc, is_new2, rid2)
	# Download (TF - target)
	cat("(1/3) TF - target table is being loaded...\n")
	outfile1 <- bfcrpath(bfc, url1)
	# Download (TF metadata)
	cat("(2/3) TF metadata is being loaded...\n")
	outfile2 <- bfcrpath(bfc, url2)
	# Loading
	tmp1 <- read.delim(outfile1, header=FALSE)
	tmp2 <- read.delim(outfile2, header=FALSE)
	colnames(tmp1) <- c("SRX_id", "GeneName")
	colnames(tmp2) <- c("SRX_id", "TF", "celltype1", "celltype2")
	# Merge
	cat("(3/3) Two files are being merged...\n")
	network_data <- left_join(tmp1, tmp2, by="SRX_id")
	network_data$TF_name <- paste(
		network_data$TF,
		network_data$celltype1,
		network_data$celltype2,
		network_data$SRX_id, sep="_")
	network_data$weight <- 1
	network_data[, c("TF_name", "GeneName", "weight")]
}

.checkCache <- function(bfc, is_new, rid){
	if(is_new){
		cat("Data file is being retrieved from remote server.\n")
	}else{
		needs_update <- bfcneedsupdate(bfc, rid)
		cache_path <- R_user_dir("BiocFileCache",
			which = "cache")
		if(needs_update){
			msg <- paste0("Local cache file in ",
				cache_path, " is old. ",
				"Data is retrived from remote server again.\n")
		}else{
			msg <- paste0("Local cache file in ",
				cache_path, " will be used.\n")
		}
		cat(msg)
	}
}
