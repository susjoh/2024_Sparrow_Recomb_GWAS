
create.dir("chr_h2_intra")

template_lines <- readLines("7.2_Chromosome_partitioning_template.R")

for(i in c(1:29)){
  
  writeLines(c(paste0("h = ", i),template_lines), paste0("chr_h2_intra/GWAS_", i, ".R")) 
  
}