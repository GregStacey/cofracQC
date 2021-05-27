# calculate QC metrics for all datasets

source("R/functions.R")
exp = read_csv("data/files.csv") %>% as.data.frame

# get corum
corum = read_tsv("data/allComplexes.txt") %>% as.data.frame()
if (!file.exists("data/corum_edges.txt")) {
  corum.edges = make.corumedgelist.by.organism("data/allComplexes.txt", "subunits(UniProt IDs)")
  write_tsv(corum.edges, file = "data/corum_edges.txt")
} else {
  corum.edges = read_tsv("data/corum_edges.txt") %>% as.data.frame()
}
corum.prots = unique(c(corum.edges[,1], corum.edges[,2]))

# calculate metrics
exp$Nid = exp$Nquant = exp$Noutliers =
  exp$mean.ngauss = exp$mean.r2 = exp$mean.sigma.maxA = exp$sd.sigma.maxA = 
  exp$N.corum.overlap = exp$corum.corr = NA
for (ii in 1:nrow(exp)) {
  chroms = get.chroms(exp$fn[ii])
  colnames(chroms) = paste0("fraction", 1:ncol(chroms))
  
  exp0 = tryCatch({
    # IDs and quantification
    exp$Nid[ii] = nrow(chroms)
    i.well.quantified = which(rowSums(!is.na(chroms))>=5)
    exp$Nquant[ii] = length(i.well.quantified)
    
    
    # modern
    autocor = calculate_autocorrelation(chroms, 10, "pearson", 1)
    zz = calculate_z_scores(autocor, bins = 1)
    exp$Noutliers[ii] = sum(abs(zz) > 2)
    chroms.modern = chroms
    chroms.modern[abs(zz)>5] = NA
    
    
    # gaussian parameters
    gauss = build_gaussians(chroms.modern, max_gaussians = 3)
    exp$mean.ngauss[ii] = mean(unlist(sapply(gauss[i.well.quantified], function(x) x$n_gaussians)))
    exp$mean.r2[ii] = mean(unlist(sapply(gauss[i.well.quantified], function(x) x$R2)))
    exp$mean.sigma.maxA[ii] = mean(unlist(sapply(gauss[i.well.quantified], function(x) {
      unlist(x$coefs$sigma[which.max(x$coefs$A)])
    })))
    exp$sd.sigma.maxA[ii] = sd(unlist(sapply(gauss[i.well.quantified], function(x) {
      unlist(x$coefs$sigma[which.max(x$coefs$A)])
    })))
    
    # corum overlap
    exp$N.corum.overlap[ii] = sum(rownames(chroms) %in% corum.prots)
    adj = corum.edges %>% dplyr::select(1:2) %>% 
      filter(.[[1]]%in%rownames(chroms) & .[[2]]%in%rownames(chroms)) %>%
      adjacency_matrix_from_data_frame()
    dd = 1 - amap::Dist(chroms, method = "pearson") %>% as.matrix()
    dd = dd[match(rownames(adj), rownames(dd)), match(rownames(adj), rownames(dd))]
    exp$corum.corr[ii] = mean(dd[adj==1])
    
    return(exp)
  },
  error=function(cond) {
    message(paste0("Failed to do QC for ", exp$fn[ii]))
    message("original error message:")
    message(cond)
    return(NA)
  })
  
  if (!is.na(exp0)) exp = exp0
  
  # write
  write_csv(exp, file = "data/files_cofracQC.csv")
}
