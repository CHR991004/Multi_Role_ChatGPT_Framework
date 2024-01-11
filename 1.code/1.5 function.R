getDirListing <- function(url) {
  # Takes a URL and returns a character vector of filenames
  while(TRUE){
    a <-  try(xml2::read_html(url), silent=F)
    if(!is(a, 'try-error')) break
  }
  fnames = grep('^G',xml_text(xml_find_all(a,'//a/@href')),value=TRUE)
  return(fnames)
}
download_fun <- function(a,i){
  while(TRUE){
    mess <-  try(download.file(a[[i]],destfile = paste0("./download_dir/",basename(a[[i]]))), silent=T)
    if(!is(mess, 'try-error')) break
  }
}
url<- function(GEO,merge,getDirListing){
  require(doParallel)
  library(stringr)
  library(GEOquery)
  library(xml2)
  library(parallel)
  stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/"
  b = getDirListing(sprintf(gdsurl, stub, GEO))
  ret <- list()
  for (x in 1:length(b))
  {
    ret[[x]]<-sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s", 
                      stub, GEO, b[x])
  }
  merge<-c(merge,unlist(ret))
  return(merge)
}