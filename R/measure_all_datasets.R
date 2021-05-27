# calculate QC metrics for all datasets

source("R/functions.R")
dir("R/modern/", full.names = T) %>% sapply(.,source,.GlobalEnv)
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
exp$Nid = exp$Nquant = exp$peak.width = 
  exp$mean.ngauss = exp$mean.r2[ii] = exp$mean.sigma.maxA = exp$sd.sigma.maxA = 
  exp$N.corum.overlap = exp$corum.corr = NA
for (ii in 1:nrow(exp)) {
  chroms = get.chroms(exp$fn[ii])
  
  # IDs and quantification
  exp$Nid[ii] = nrow(chroms)
  i.well.quantified = which(rowSums(!is.na(chroms))>=5)
  exp$Nquant[ii] = length(i.well.quantified)
  
  
  # modern
  zz = detect_outliers(chroms)
  exp$Noutliers[ii] = sum(abs(zz) > 5)
  
  # gaussian parameters
  gauss = build_gaussians(chroms.modern, max_gaussians = 3)
  exp$mean.ngauss[ii] = mean(unlist(sapply(gauss[i.well.quantified], function(x) x$n_gaussians)))
  exp$mean.r2[ii] = mean(unlist(sapply(gauss[i.well.quantified], function(x) x$R2)))
  exp$mean.sigma.maxA[ii] = mean(unlist(sapply(gauss[i.well.quantified], function(x) {
    unlist(x$coefs$sigma[which.max(x$coefs$A)])
  })))
  exp$sd.sigma.maxA[ii] = mean(unlist(sapply(gauss[i.well.quantified], function(x) {
    unlist(x$coefs$sigma[which.max(x$coefs$A)])
  })))
  
  
  # corum overlap
  exp$N.corum.overlap[ii] = sum(rownames(chroms) %in% corum.prots)
  
}
