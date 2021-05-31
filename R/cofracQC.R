# calculate QC metrics for a single dataset and compare to others
source("R/functions.R")

cofraqcQC = function(fn = NULL, fn.out = NULL) {
  # read results
  fn.results = "data/cofrac_benchmark.csv"
  if (!file.exists(fn.results)) {
    exp = dir("data/files_cofracQC/", 
              pattern = "files_cofracQC_", full.names = T) %>%
      map(., read_csv) %>% do.call(rbind, .) %>% as.data.frame()
    write_csv(exp, file = fn.results)
  } else exp = read_csv(fn.results) %>% as.data.frame()
  
  #if (is.null(fn.out)) fn.out = paste0(file_path_sans_ext(fn), "_benchmark.png")
  fn.out = paste0("figures/", paste(unique(get.group.from.fn(fn)),collapse = "-"), "_benchmark.png")
  
  # 
  if (!is.null(fn)) {
    fn.group = get.group.from.fn(fn)
    fn.group[is.na(fn.group)] = basename(fn[is.na(fn.group)])
    ia = which(get.group.from.fn(exp$fn) %in% fn.group)
    if (length(ia)>0) {
      # is fn already in the results?
      this.results = exp[ia,]
      exp = exp[-ia,]
    } else {
      # make results for fn
      this.results = calculate_measures(fn)
    }
  }
  
  # plot
  this.theme = theme_minimal() + 
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_blank())
  pp.invalid = scale_fill_grey()
  pp.valid = scale_fill_viridis_d()
  corum.valid = pp.valid
  if (this.results$N.corum.overlap==0) corum.valid = pp.invalid
  
  # ids, quantification, outlier
  p1 = ggplot(exp, aes(x=Nid, fill=source)) + geom_density(alpha = .5) +
    geom_vline(data=this.results, aes(xintercept=Nid), color="red", size=2, alpha = .5) +
    theme_minimal() + scale_fill_viridis_d() + ggtitle("Number of IDs") + this.theme
  p2 = ggplot(exp, aes(x=Nquant/Nid, fill=source)) + geom_density(alpha = .5) +
    geom_vline(data=this.results, aes(xintercept=Nquant/Nid), color="red", size=2, alpha = .5) +
    theme_minimal() + scale_fill_viridis_d() + ggtitle("fraction quantified IDs") + this.theme
  p3 = ggplot(exp, aes(x=Noutliers/Nid, fill=source)) + geom_density(alpha = .5) +
    geom_vline(data=this.results, aes(xintercept=Noutliers/Nid), color="red", size=2, alpha = .5) +
    theme_minimal() + scale_fill_viridis_d() + ggtitle("fraction modern outliers") + this.theme
  
  # corum 
  p4 = ggplot(exp, aes(x=corum.corr, fill=source)) + geom_density(alpha = .5) +
    geom_vline(data=this.results, aes(xintercept=corum.corr), color="red", size=2, alpha = .5) +
    theme_minimal() + corum.valid + ggtitle("Average CORUM pairwise correlation") + this.theme +
    theme(legend.position = "right")
  p9 = ggplot(exp, aes(x=N.corum.overlap/Nid, fill=source)) + geom_density(alpha = .5) +
    geom_vline(data=this.results, aes(xintercept=N.corum.overlap/Nid), color="red", size=2, alpha = .5) +
    theme_minimal() + corum.valid + ggtitle("Fraction of proteins in CORUM") + this.theme
  
  # gauss parameters
  p5 = ggplot(exp, aes(x=mean.ngauss, fill=source)) + geom_density(alpha = .5) +
    geom_vline(data=this.results, aes(xintercept=mean.ngauss), color="red", size=2, alpha = .5) +
    theme_minimal() + scale_fill_viridis_d() + ggtitle("avg number of gaussians") + this.theme
  p6 = ggplot(exp, aes(x=mean.r2, fill=source)) + geom_density(alpha = .5) +
    geom_vline(data=this.results, aes(xintercept=mean.r2), color="red", size=2, alpha = .5) +
    theme_minimal() + scale_fill_viridis_d() + ggtitle("avg model fit (R^2)") + this.theme
  p7 = ggplot(exp, aes(x=mean.sigma.maxA, fill=source)) + geom_density(alpha = .5) +
    geom_vline(data=this.results, aes(xintercept=mean.sigma.maxA), color="red", size=2, alpha = .5) +
    theme_minimal() + scale_fill_viridis_d() + ggtitle("mean peak width (fractions)") + this.theme
  p8 = ggplot(exp, aes(x=sd.sigma.maxA, fill=source)) + geom_density(alpha = .5) +
    geom_vline(data=this.results, aes(xintercept=sd.sigma.maxA), color="red", size=2, alpha = .5) +
    scale_fill_viridis_d() + ggtitle("peak width variability (fractions)")  + this.theme
  
  if (!is.null(fn.out)) ggsave(fn.out, 
                               plot = (p1 | p2 | p3 | p4 ) / (p5 | p6 |p7 |p8 | p9), 
                               width = 16)
  return((p1 | p2 | p3 | p4 ) / (p5 | p6 |p7 |p8 | p9))
}


