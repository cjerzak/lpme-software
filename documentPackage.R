{
  rm(list=ls()); options(error = NULL)

  # set path and specify package name
  setwd(sprintf("~/Documents/%s-software", package_name <- "lpmec"))

  # get version number from DESCRIPTION
  package_path <- sprintf("~/Documents/%s-software/%s", package_name, package_name)
  versionNumber <- read.dcf(file.path(package_path, "DESCRIPTION"), fields = "Version")[1, 1]

  # document package
  tools::add_datalist(package_path <- sprintf("~/Documents/%s-software/%s",package_name,package_name),
                      force = TRUE, small.size = 1L)
  devtools::document(package_path)

  # run tests (stop on failure)
  test_results <- devtools::test(package_path)
  if (any(as.data.frame(test_results)$failed > 0)) {
    stop("Tests failed! Stopping build process.")
  }
  cat("\n✓ All tests passed!\n\n") 

  # remove old PDF
  try(file.remove(sprintf("./%s.pdf",package_name)),T)
 
  # create new PDF
  system(sprintf("R CMD Rd2pdf %s",package_path))
  
  # build tar
  system( paste(shQuote(file.path(R.home("bin"), "R")),
                "R CMD build --resave-data", shQuote(package_path)) )
  
  # check as cran
  system( paste(shQuote(file.path(R.home("bin"), "R")),
                "R CMD check --as-cran",
                shQuote(
                  paste(package_name, "_", versionNumber, ".tar.gz", sep = "")
                ))  )

  # Check package to ensure it meets CRAN standards.
  # devtools::check( package_path )
  
  # check data integrity
  # devtools::install_github(repo = "cjerzak/lpmec-software/lpmec")
  # data(package = "lpmec")
  # data("KnowledgeVoteDuty",package="lpmec")

  # install current local build
  install.packages( "~/Documents/lpmec-software/lpmec",repos = NULL, type = "source",force = F) # install from local
}

