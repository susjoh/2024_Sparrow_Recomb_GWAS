
create.dir("chr_h2_perm")

template_lines <- readLines("8.2_Chromosome_partitioning_template.R")

for(i in c(1:29)){
  
  writeLines(c(paste0("h = ", i),template_lines), paste0("chr_h2_perm/Perm_", i, ".R")) 
  
}