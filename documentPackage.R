{
  rm(list=ls()); options(error = NULL)
  # install.packages("~/Documents/lpme-software/lpme",repos = NULL, type = "source",force = F)
  #data(package = "lpme")
  data("KnowledgeVoteDuty",package="lpme")
  
  # set path and specify package name
  setwd(sprintf("~/Documents/%s-software", package_name <- "lpme"))

  # document package
  tools::add_datalist(package_path <- sprintf("~/Documents/%s-software/%s",package_name,package_name),
                      force = TRUE, small.size = 1L)
  devtools::document(package_path)

  # remove old PDF
  try(file.remove(sprintf("./%s.pdf",package_name)),T)
 
  # create new PDF
  system(sprintf("R CMD Rd2pdf %s",package_path))

  # Check package to ensure it meets CRAN standards.
  # devtools::check( package_path )
}

