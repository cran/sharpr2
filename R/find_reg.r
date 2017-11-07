#' find_reg
#'
#' Given an object from sharpr2 and a position, this function finds the tiled region containing the position.
#' @param re an object obtained from sharpr2.
#' @param pos a position for which the tiled region is searched.  
#' @keywords sharpr2 HiDRA
#' @return ind: the index of the tiled region in the object from sharpr2. If no such tile region is found, NA is returned.
#' @export
#' @examples
#' # find_reg(re, 100000)


find_reg <- function(re, pos)
{
	if(class(re)!="sharpr2")
	{
		stop('The first argument must be an object obtained from sharpr2.')
	}
	
	s_e <- unlist(strsplit(unlist(re$region),'-'))
	start <- as.numeric(s_e[seq(from=1,to=length(s_e),by=2)])
	end <- as.numeric(s_e[seq(from=2,to=length(s_e),by=2)])
	ind_s <- pos>=start
	ind_e <- pos<=end
	ind <- which((ind_s+ind_e)==2)
	if(length(ind)==0)	
	{
		ind <- NA
	}
	
	return(ind)
}
