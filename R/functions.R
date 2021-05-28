options(stringsAsFactors = F)
library(argparse)
library(tidyverse)
library(magrittr)
library(PrInCE)
library(flavin)
library(dplyr)
library(magrittr)
require(reshape2)
require(igraph)
require(LiblineaR)
require(amap)
require(patchwork)
require(tools)
require(gtools)
dir("R/modern/", full.names = T) %>% sapply(.,source,.GlobalEnv)

taxa = c("human"=9606, "mouse"=10090)
fn.proteome = c("human"="../no-reference/data/proteomes/uniprot-human-9606.txt",
                "mouse"="../no-reference/data/proteomes/uniprot-mouse-10090.txt",
                "bee"=NULL)
gold.standard = c("human"="../no-reference/data/complex/CORUM_human_edgelist.txt",
                  "mouse"="../no-reference/data/complex/CORUM_mouse_edgelist.txt",
                  "yeast"="../no-reference/data/complex/benschop_yeast_edgelist.txt",
                  "og"="../no-reference/data/complex/CORUM_edgelist_OG.txt",
                  "bee"="../no-reference/data/complex/CORUM_edgelist_OG.txt",
                  "fake01"="data/complex/fakechroms_01_edgelist.txt")
# submit a job to the appropriate cluster
submit_job = function(jobs0, script, allocation,
                      system = c('cedar', 'sockeye', 'elasti'),
                      jobs_per_array = 100) {
  system = match.arg(system)
  
  if (is.character(jobs0)){
    jobs = as.data.frame(read_csv(jobs0))
    Njobs = nrow(jobs)+1
  } else if(is.numeric(jobs0)) {
    Njobs = jobs0+1
  }
  
  
  if (system == 'cedar') {
    system(paste0("cd ~/project; ",
                  "sbatch --account=", allocation, " --array=1-", Njobs, 
                  " ", script))
  } else if (system == 'elasti') {
    system(paste0("cd ~/project; ",
                  "sbatch --account=", allocation, " --array=1-", Njobs,
                  " ", script))
  } else if (system == 'sockeye') {
    n_jobs = Njobs
    ## Sockeye only lets you run 1,000 jobs at a time
    n_submissions = ifelse(n_jobs > jobs_per_array,
                           ceiling(n_jobs / jobs_per_array), 1)
    for (submission_idx in seq_len(n_submissions)) {
      job_start = (submission_idx - 1) * jobs_per_array + 1
      job_end = ifelse(submission_idx == n_submissions,
                       ifelse(n_jobs %% jobs_per_array == 0,
                              submission_idx * jobs_per_array,
                              job_start - 1 + n_jobs %% jobs_per_array),
                       submission_idx * jobs_per_array)
      system(paste0("qsub -A ", allocation, " -J ", job_start, "-", job_end,
                    " ", script))
    }
  }
}


get.chroms = function(fn) {
  chroms = as.data.frame(read_tsv(fn)) %>% filter(!grepl("__", .[,1]))
  
  # remove replicate column if needed
  i.rep = which(tolower(names(chroms)) %in% "replicate")
  if (length(i.rep)>0){
    if (length(unique(chroms[,i.rep]))>1) {
      stop("this function can't handle multiple replicates.")
    }
    print(paste0(fn, " has a replicate column. Removing it..."))
    chroms = chroms[,-i.rep]
  }
  
  prots = chroms[,1]
  prots = unlist(sapply(strsplit(prots, ";"), "[", 1))
  chroms = as.matrix(chroms[,2:ncol(chroms)])
  rownames(chroms) = prots
  # log-transform
  #if (mean(chroms, na.rm=T)>1000) {
  #  chroms = log10(chroms)
  #  chroms[is.infinite(chroms)] = NA
  #}
  return(chroms)
}


# return a dataframe of corum edges with organism annotation
make.corumedgelist.by.organism = function(fn, column.name) {
  complexes = read_tsv(fn) %>% dplyr::select(c(column.name, "Organism")) %>% as.data.frame() %>%
    `colnames<-`(c("complex","organism")) %>%
    filter(grepl(";", complex))
  edgelist = sapply(1:nrow(complexes), function(x) {
    x = unlist(strsplit(complexes$complex[x], ";")) %>% combn(., 2) %>% t()
  }) %>% lapply(., as.data.frame)
  edgelist = lapply(1:length(edgelist), function(x) {edgelist[[x]] %<>% mutate(organism = complexes$organism[x])}) %>%
    do.call(rbind, .)
  return(edgelist)
}

