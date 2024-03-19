
dir.create("chr_h2_rintra")

template_lines <- readLines("7.5_Chromosome_partitioning_template.R")

for(i in c(1:29)){
  
  writeLines(c(paste0("h = ", i),template_lines), paste0("chr_h2_rintra/GWAS_", i, "_rintra.R")) 
  
}