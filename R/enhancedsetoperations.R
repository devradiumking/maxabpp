#'@export
union <- function(x,y,...){
  dat <- list(x,y,...)
  if(length(dat) < 2) return(unlist(dat))
  u <- as.vector(dat[[1]])
  for(i in 2:length(dat)){
    u <- unique(c(u,as.vector(dat[[i]])))
  }
  u
}
#'@export
intersect <- function(x,y,...){
  dat <- list(x,y,...)
  if(length(dat) < 2) return(unlist(dat))
  common <- as.vector(dat[[1]])
  for(i in 2:length(dat)){
    common <- unique(common[match(as.vector(dat[[i]]), common, 0L)])
    if(length(common) == 0)
      break
  }
  common
}
#'@export
#list all possible intersections
quick_venn_table=function(x){
  #return Venn diagram entry sizes
  #x: a list of sets
  if(!is.list(x)) stop('Input x must be list\n')
  nL=length(x)
  if(nL<2) stop('Input x should have at least two entries\n')
  allE=unique(unlist(x))
  venn_areas=rep('',length(allE))
  for(i in 1:nL){
    venn_areas=paste(venn_areas,ifelse(allE %in% x[[i]],'1','0'),sep='')
  }
  Res=data.frame(Entry=allE,Venn.area=venn_areas,stringsAsFactors=FALSE)
 return(Res)
}
#' Perform grepl on multiple patterns; it's like  AND-ing or OR-ing successive
#'  grepl statements.
#'
#' @export
#' @param pattern character vector of patterns
#' @param x the character vector to search
#' @param op logical vector operator back quoted, defaults to `|`
#' @param ... further arguments for \code{grepl} like \code{fixed} etc.
#' @return logical vector
mgrepl <- function(pattern, x, op = `|`, ... ){
  Reduce(op, lapply(pattern, grepl, x, ...))
}
