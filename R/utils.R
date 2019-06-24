#' A utility function to convert a \code{runjags} or \code{runjags.list} object
#' to a tibble
#' @param obj the object to be converted
#' @param paramPrefix the name of any indexed model parameters, including the
#' opening square bracket. eg 'prob['
#' @param ... arbitrary arguments that were passed to runjags when \code{obj}
#' was created
#' @return the corresponding tibble
runjagsResultsToTibble <- function(obj, paramPrefix, ...)
{
  if (!("runjags" %in% class(obj)) & !("runjags.list" %in% class(obj)))
    stop("obj must be a runjags or runjags.list object")
  rv <- list(length(obj$mcmc))
  for (i in 1:length(obj$mcmc))
  {
    rv[[i]] <- tibble::as_tibble(obj$mcmc[[i]]) %>%
      tibble::add_column(Chain=i,
                         Sample=seq(from=stats::start(obj$mcmc[[i]]), 
                                    to=stats::end(obj$mcmc[[i]]),
                                    #deltat(mcmc) returns 1 even when thin != 1.
                                    by=ifelse(is.null(list(...)[["thin"]]),
                                              1, 
                                              list(...)[["thin"]])))
    indexedTibble <- rv[[i]] %>%
      dplyr::select(Chain, Sample, tidyselect::starts_with(paramPrefix)) %>% 
      tidyr::gather(key=Column, value=Value, tidyselect::starts_with(paramPrefix)) %>%
      tidyr::separate(Column,
                      into=c("Parameter", "Index", "Dummy"),
                      sep="[\\[,\\]]",
                      convert=TRUE,
                      retain=FALSE) %>%
      dplyr::select(-Dummy) %>% 
      dplyr::arrange(Chain, Sample, Parameter)
    unindexedTibble <- rv[[i]] %>%
      dplyr::select(-tidyselect::starts_with(paramPrefix))
    if (ncol(unindexedTibble) > 2)
    {
      unindexedTibble <- unindexedTibble %>% 
        tidyr::gather(key=Parameter, value=Value, -Chain, -Sample) %>%
        dplyr::mutate(Index=NA) %>% 
        dplyr::arrange(Chain, Sample, Parameter) %>% 
        dplyr::select(Chain, Sample, Parameter, Index, Value)
      rv[[i]] <- indexedTibble %>% 
        dplyr::full_join(unindexedTibble, 
                         by=c("Chain", "Sample", "Parameter", "Index", "Value"))
    } else rv[[i]] <- indexedTibble
  }
  return (dplyr::bind_rows(rv))
}
#' Parse a model string, repacing tokens as required
#'
#' Parsing is recursive.
#'
#' @param s the string to be parsed
#' @param tokens a list whose names define the tokens to be replaced and whose
#' values define the replacements
#' @return the parsed string
#' @seealso \code{\link{getModelString}}, \code{\link{fitBinomialModel}}
#' @examples
#' #Do nothing
#' s <- parseModelString(s)
#' #Edit the supplied string so that a (binomial) model is fitted to the tibble
#' \{myTibble} in which groups are defined by the variable \code{Centre},
#' observation counts by \code{N} and event counts by \code{R}.
#' s <- parseModelString(s, list="<!DATA!>="myTibble", "<!GROUP!>="Centre",
#' "<!COUNT!>="N", "<!EVENTS!>!="R")
# parseModelString <- function(s, tokens=list())
# {
#   oldString <- s
#   #s <- mapply(stringr::str_replace_all, names(tokens), tokens)
#   for (key in names(tokens))
#     s <- stringr::str_replace_all(s, key, tokens[[key]])
#   #Recurse until no change
#   if (oldString != s) s <- parseModelString(s, tokens)
#   return(s)
# }
