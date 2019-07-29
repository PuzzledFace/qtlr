#' Fit a Normal model to subject-level data
#'
#' @param x a vector of observed data 
#' @param g a vector indicating the group membership for each datum.  Either a 
#' factor defining group membership, or a vector of integers between 1 and the
#' number of groups.
#' @param inits a list containing the initial values of the hyperparameters.
#' See Usage Notes below.
#' @param modelString The JAGS model string that defines the model to be fitted
#' @param autoRun Logical.  If \code{TRUE}, use \code{autorun.jags}.  Otherwise,
#' use \code{run.jags}.
#' @param raw logical if \code{TRUE}, the \code{mcmc} object created by
#'  \code{run.jags()} is returned.  Otherwise, a \code{tibble} containing the
#' concatenation of elements of the \code{mcmc.list} created by JAGS
#' is returned
#' @param ... passed to JAGS
#' @return Either the \code{mcmc} (or \code{mcmc.list}) object returned by
#' JAGS or a tibble containing the MCMC samples from the posterior
#' distribution
#' @section Usage notes:
#' If \code{modelString == NULL}, the model string is obtained by calling
#' \code{getModelString("tte")}.\cr\cr  
#' If \code{raw == FALSE}, the chain from which each observation is drawn is
#' indicated by \code{Chain} and the dataset is transformed into tidy format,
#' with the model parameter indicated by \code{Parameter}.\cr\cr  
#' The \code{inits} parameter can be used to define the number of chains created
#' my JAGS.  If a list of lists, the number of elements in the
#' outer list defines the number of chains and the elements of each sub-list
#' define the initial value for each hyperparameter.  For example, the default
#' value of \code{inits} requests two chains.  The initial values of \code{shape}
#'  and \code{scale} in the first chain are both 0.0001.  In the second
#' chain, they are both 1.  If the MCMC model has converged and is stationary, 
#' the initial values of the hyperparameters will be irrelevant.   To check for 
#' convergence, it is necessary - but not sufficient - to obtain more than one
#' chain and to use different initial values for each chain.\cr\cr
#' Note that some parameters to \code{run.jags} are incompatible with 
#' \code{autorun.jags} and will therefor cause an error when passed in 
#' \code{...} unless \code{autoRun} is set to \code{FALSE}.  These parameters
#' include \code{samples} and \code{burnin} and will result in an error similar
#' to "Error in \code{extend.jags(runjags.object, add.monitor = add.monitor,
#' drop.monitor = drop.monitor,}  : formal argument "burnin" matched by multiple
#' actual arguments".
#' @seealso \code{\link{fitBinomialModel}}, \code{\link{fitPoissonModel}},
#' \code{\link{fitBinaryModel}}, \code{\link{fitTteModel}}
#' @examples
#' #Fit a Normal QTL model to artifical data
#' nCentres <- 6
#' centreSizes <- ceiling(runif(nCentres, min=8, max=25))
#' group <- rep(1:nCentres, times=centreSizes)
#' centreMeans <- rnorm(nCentres, mean=5, sd=1.5)
#' means <- rep(centreMeans, times=centreSizes)
#' x <- rnorm(length(group), means, 3)
#' fitNormalModel(x, group) 
#' @export
fitNormalModel <- function(x,
                           g,
                           inits=NULL,
                           modelString=NULL,
                           autoRun=TRUE,
                           raw=FALSE,
                           ...)
{
  #Validate
  if (length(x) != length(x)) stop ("x and g must be of the same length.")
  if (!is.factor(g)) 
  {
    if (any(g > length(g))) stop("Invalid group index found in g")
    if (any(g < 0)) stop("Invalid group index found in g")
  }
  
  #Initialise
  ###Prepare model input
  d <- list()
  d$x <- c(x, NA)
  d$g <- c(g, max(g, na.rm=TRUE)+1)
  d$n <- length(d$x)
  d$k <- length(unique(d$g))
  
  if (is.null((inits)))
  {
    inits <- list()
    inits1 <- list(mu=rep(0, d$k), invTau=1e06)
    inits2 <- list(mu=stats::rnorm(d$k, 
                                   mean(0), 
                                   sd=sqrt(stats::var(x))),
                   invTau=1e03)
    inits <- list(inits1, inits2)
  }
  if (is.null(modelString)) modelString <- getModelString("normal")
  
  #Execute
  if (autoRun) f <- runjags::autorun.jags
  else f <- runjags::run.jags
  results <- f(model=modelString,
               monitor=c("mu", "tau"),
               data=d,
               inits=inits,
               ...)
  if (raw) return (results)
  results <- runjagsResultsToTibble(results, "mu[", ...)
  return (results)
}

#' Fit a time-to-event QTL model
#' 
#' @param time a vector containing the event/censoring times for each datum.
#' @param status a vector containing the censoring status for each datum.
#' @param censorFlag a scalar or vector containing the values of \code{status}
#' that indicate that the corresponding event time is censored.
#' @param g a vector indicating the group membership for each datum.  Either a 
#' factor defining group membership, or a vector of integers between 1 and the
#' number of groups.
#' @param inits a list containing the initial values of the hyperparameters.
#' See Usage Notes below.
#' @param modelString The JAGS model string that defines the model to be fitted
#' @param autoRun Logical.  If \code{TRUE}, use \code{autorun.jags}.  Otherwise,
#' use \code{run.jags}.
#' @param raw logical if \code{TRUE}, the \code{mcmc} object created by
#'  \code{run.jags()} is returned.  Otherwise, a \code{tibble} containing the
#' concatenation of elements of the \code{mcmc.list} created by JAGS
#' is returned
#' @param ... passed to JAGS
#' @return Either the \code{mcmc} (or \code{mcmc.list}) object returned by
#' JAGS or a tibble containing the MCMC samples from the posterior
#' distribution
#' @section Usage notes:
#' If \code{modelString == NULL}, the model string is obtained by calling
#' \code{getModelString("tte")}.\cr\cr  
#' If \code{raw == FALSE}, the chain from which each observation is drawn is
#' indicated by \code{Chain} and the dataset is transformed into tidy format,
#' with the model parameter indicated by \code{Parameter}.\cr\cr  
#' The \code{inits} parameter can be used to define the number of chains created
#' my JAGS.  If a list of lists, the number of elements in the
#' outer list defines the number of chains and the elements of each sub-list
#' define the initial value for each hyperparameter.  For example, the default
#' value of \code{inits} requests two chains.  The initial values of \code{shape}
#'  and \code{scale} in the first chain are both 0.0001.  In the second
#' chain, they are both 1.  If the MCMC model has converged and is stationary, 
#' the initial values of the hyperparameters will be irrelevant.   To check for 
#' convergence, it is necessary - but not sufficient - to obtain more than one
#' chain and to use different initial values for each chain.\cr\cr
#' Note that some parameters to \code{run.jags} are incompatible with 
#' \code{autorun.jags} and will therefor cause an error when passed in 
#' \code{...} unless \code{autoRun} is set to \code{FALSE}.  These parameters
#' include \code{samples} and \code{burnin} and will result in an error similar
#' to "Error in \code{extend.jags(runjags.object, add.monitor = add.monitor,
#' drop.monitor = drop.monitor,}  : formal argument "burnin" matched by multiple
#' actual arguments".
#' @seealso \code{\link{fitBinomialModel}}, \code{\link{fitPoissonModel}},
#' \code{\link{fitBinaryModel}}, \code{\link{fitNormalModel}}
#' @examples 
#' #Fit a QTL model to the VA lung dataset reported by Kalbfleisch and 
#' #Prentice (1980), using cell type to define group membership
#' fitTteModel(time=vaLung$SurvivalTime,
#'             status=vaLung$Status,
#'             g=vaLung$CellType)
#' @export
fitTteModel <- function(time,
                        status,
                        censorFlag=0,
                        g,
                        inits=NULL,
                        modelString=NULL,
                        autoRun=TRUE,
                        raw=FALSE,
                        ...)
{
  #Validate
  if (length(time) != length(status)) stop ("time and status must be of the same length.")
  if (length(time) != length(g)) stop ("time and g must be of the same length.")
  if (!is.factor(g)) 
  {
    if (any(g > length(g))) stop("Invalid group index found in g")
    if (any(g < 0)) stop("Invalid group index found in g")
  }

  #Initialise
  ##Calculate by-group summary statistics
  inputData <- tibble::tibble(Group=g,
                              Time=time,
                              Censored=(status %in% censorFlag))

  summaryData <- inputData %>% 
                   dplyr::group_by(Group) %>% 
                   dplyr::summarise(N=dplyr::n(),
                                    Events=N-sum(Censored),
                                    Censored=sum(Censored),
                                    TotalTime=sum(Time),
                                    MeanTime=TotalTime/Events)
  ##Prepare model input
  d <- list()
  d$logMean <- c(log(summaryData$MeanTime), NA)
  d$n <- c(summaryData$N, 1)
  d$k <- length(d$n)

  if (is.null((inits)))
  {
    inits <- list()
    inits1 <- list(mu=rep(mean(summaryData$MeanTime), d$k), invTau=1e06)
    inits2 <- list(mu=stats::rnorm(d$k, 
                                   mean(summaryData$MeanTime), 
                                   sd=sqrt(stats::var(summaryData$MeanTime))),
                   invTau=1e03)
    inits <- list(inits1, inits2)
  }
  if (is.null(modelString)) modelString <- getModelString("tte")
  
  #Execute
  if (autoRun) f <- runjags::autorun.jags
  else f <- runjags::run.jags
  results <- f(model=modelString,
               monitor=c("mu", "tau"),
               data=d,
               inits=inits,
               ...)
  if (raw)
  {
    for(i in 1:length(x))
      for(j in 1:d$k) 
        x[[i]][ ,j] <- exp(x[[i]][ ,j])
    return (results)
  }
  results <- runjagsResultsToTibble(results, "mu[", ...)
  results <- results %>% dplyr::mutate(Value=exp(Value))
  return (results)
}

#' Fit a binary QTL model.
#' 
#' @param g the vector holding the group ID for this observation
#' @param r the vector holding the response status of this observation
#' @param inits a list containing the initial values of the hyperparameters.
#' See Usage Notes below.
#' @param modelString The JAGS model string that defines the model to be fitted
#' @param autoRun Logical.  If \code{TRUE}, use \code{autorun.jags}.  Otherwise,
#' use \code{run.jags}.
#' @param raw logical if \code{TRUE}, the \code{mcmc} object created by
#'  JAGS is returned.  Otherwise, a \code{tibble} containing the
#' concatenation of elements of the \code{mcmc.list} created by JAGS
#' is returned
#' @param ... passed to JAGS
#' @return Either the \code{mcmc} (or \code{mcmc.list}) object returned by
#' JAGS or a tibble containing the MCMC samples from the posterior
#' distribution
#' @examples
#' #Simple use
#' b <- createBerryData(binary=TRUE) 
#' m <- fitBinaryModel(b$Event, b$Study) %>% dplyr::filter(Index == 10)
#' #Passing parameters to run.jags()
#' inits1 <- list(a=4, b=2)
#' inits2 <- list(a=1, b=1)
#' inits3 <- list(a=2, b=10)
#' m <- fitBinaryModel(b$Event, b$Study,
#'                     inits=list(inits1, inits2, inits3),
#'                     thin=2)
#' @section Usage notes:
#' If \code{modelString == NULL}, the model string is obtained by calling
#' \code{getModelString("binary")}.\cr\cr  
#' If \code{raw == FALSE}, the chain from which each observation is drawn is
#' indicated by \code{Chain} and the dataset is transformed into tidy format,
#' with the model parameter indicated by \code{Parameter}.\cr\cr  
#' The \code{inits} parameter can be used to define the number of chains created
#' by JAGS.  If a list of lists, the number of elements in the
#' outer list defines the number of chains and the elements of each sub-list
#' define the initial value for each hyperparameter.  For example, the default
#' value of \code{inits} requests two chains.  The initial values of \code{a}
#'  and \code{b} in the first chain are 2 and 4, respectively.  In the second
#' chain, the corresponding values are 9 and 1.  If the MCMC model has converged
#' and is stationary, the initial values of the hyperparameters will be
#' irrelevant.   To check for convergence, it is necessary - but not sufficient
#' - to obtain more than one chain and to use different initial values for each
#' chain.\cr\cr
#' Note that some parameters to \code{run.jags} are incompatible with 
#' \code{autorun.jags} and will therefor cause an error when passed in 
#' \code{...} unless \code{autoRun} is set to \code{FALSE}.  These parameters
#' include \code{samples} and \code{burnin} and will result in an error similar
#' to "Error in \code{extend.jags(runjags.object, add.monitor = add.monitor,
#' drop.monitor = drop.monitor,}  : formal argument "burnin" matched by multiple
#' actual arguments".
#' @seealso \code{\link{fitBinomialModel}}, \code{\link{fitPoissonModel}},
#' \code{\link{fitTteModel}}, \code{\link{fitNormalModel}}
#' @export
fitBinaryModel <- function(r,
                           g,
                           inits=list(list("a"=2, "b"=4), list("a"=9, "b"=1)),
                           modelString=NULL,
                           autoRun=TRUE,
                           raw=FALSE,
                           ...)
{
  #Validation
  if (length(r) != length(g)) stop ("r and g must be of the same length.")
  if (!all.equal(as.integer(g), g)) stop("Non-integer values found in g.")
  if (is.logical(r)) r <- as.numeric(r)
  if (!all.equal(as.integer(r), r)) stop("Non-integer values found in r.")
  if (any(g > length(g))) stop("Invalid group index found in g")
  if (any(g < 0)) stop("Invalid group index found in g")
  if (any(r < 0)) stop("Negative value(s) found in r.")
  if (any(r > 1)) stop("Invalid response status found in r")

  #Initialise
  if (is.null(modelString) | length(modelString) == 0)
    modelString <- getModelString("binary")

  #Begin
  #Append dummy, "posterior", observation
  d <- list()
  d$group <- c(g, length(unique(g))+1)
  d$r <- c(as.numeric(r), NA)
  d$k <- length(d$group)
  d$m <- length(unique(d$group))
  #Execute
  if (autoRun) f <- runjags::autorun.jags
  else f <- runjags::run.jags
  results <- runjags::autorun.jags(model=modelString,
                                   monitor=c(c("a", "b"),
                                             paste0("p[", 1:d$m , "]")),
                                   data=d,
                                   inits=inits,
                                   ...)
  if (raw) return (results)
  results <- runjagsResultsToTibble(results, "p[", ...)
  return (results)
}

#' Fit a Poisson QTL model to rate or event count data
#' @param e a vector containing the number of events in this group
#' @param t a vector containing the total exposure in this group
#' @param modelString the string containing the JAGS model definition
#' @param inits a list containing the initial values of the hyperparameters.
#' See Usage Notes below.
#' @param autoRun Logical.  If \code{TRUE}, use \code{autorun.jags}.  Otherwise,
#' use \code{run.jags}.
#' @param raw logical if \code{TRUE}, the \code{mcmc} object created by
#'  \code{run.jags()} is returned.  Otherwise, a \code{tibble} containing the
#' concatenation of elements of the \code{mcmc.list} created by JAGS
#' is returned
#' @param ... passed to JAGS
#' @return Either the \code{mcmc} (or \code{mcmc.list}) object returned by
#' JAGS or a tibble containing the MCMC samples from the posterior
#' distribution
#' @section Usage notes:
#' If \code{modelString == NULL}, the model string is obtained by calling
#' \code{getModelString("binary")}.\cr\cr  
#' If \code{raw == FALSE}, the chain from which each observation is drawn is
#' indicated by \code{Chain} and the dataset is transformed into tidy format,
#' with the model parameter indicated by \code{Parameter}.\cr\cr  
#' The \code{inits} parameter can be used to define the number of chains created
#' by JAGS.  If a list of lists, the number of elements in the
#' outer list defines the number of chains and the elements of each sub-list
#' define the initial value for each hyperparameter.  For example, the default
#' value of \code{inits} requests two chains.  The initial values of \code{a}
#'  and \code{b} in the first chain are 2 and 4, respectively.  In the second
#' chain, the corresponding values are 9 and 1.  If the MCMC model has converged
#' and is stationary, the initial values of the hyperparameters will be
#' irrelevant.   To check for convergence, it is necessary - but not sufficient
#' - to obtain more than one chain and to use different initial values for each
#' chain.\cr\cr
#' Note that some parameters to \code{run.jags} are incompatible with 
#' \code{autorun.jags} and will therefor cause an error when passed in 
#' \code{...} unless \code{autoRun} is set to \code{FALSE}.  These parameters
#' include \code{samples} and \code{burnin} and will result in an error similar
#' to "Error in \code{extend.jags(runjags.object, add.monitor = add.monitor,
#' drop.monitor = drop.monitor,}  : formal argument "burnin" matched by multiple
#' actual arguments".
#' @seealso \code{\link{fitBinomialModel}}, \code{\link{fitBinaryModel}},
#' \code{\link{fitTteModel}}, \code{\link{fitNormalModel}}
#' @export
fitPoissonModel <- function(e, 
                            t, 
                            modelString, 
                            inits=list(list("shape"=1, "scale"=1),
                                       list("shape"=2, "scale"=10)),
                            autoRun=TRUE,
                            raw=FALSE,
                            ...)
{
  #Validate
  if (length(e) != length(t)) stop ("e and t must be of the same length.")
  if (!all.equal(as.integer(e), e)) stop("Non-integer values found in e.")
  if (any(t < 0)) stop("Negative value(s) found in t.")
  if (any(is.na(t))) stop("NA(s) found in t.")
  if (any(is.na(e))) stop("NA(s) found in e.")
  if (any(is.null(t))) stop("NULL(s) found in t.")
  if (any(is.null(e))) stop("NULL(s) found in e.")
  #Initialise
  #Add additional, "posterior" observation to each input vector
  d <- list()
  d$events <- c(e, NA)
  d$exposure <- c(t, 1)
  d$k <- length(d$events)
  #Execute
  modelString <- getModelString("poisson")
  if (autoRun) f <- runjags::autorun.jags
  else f <- runjags::run.jags
  results <- f(modelString,
               monitor=c(c("scale", "shape"), 
                         paste0("mu[", 1:d$k , "]")),
               data=d,
               inits=inits,
               ...)
  if (raw) return(results)
  else return (runjagsResultsToTibble(results, "mu["))
}

#'Obtain the quantile based on a sampled posterior
#'@param probs a vector containing the required quantiles
#'@param data a tibble containing the MCMC samples from which the quantile
#' should be estimated
#'@param type the quantile type to be returned.  See the help for \code{quantile()} for
#' more information.  Default 4, to ensure quality with sampled values.
#'@return a vector containing the quantiles obtained from the Value column in 
#' \code{data}
#'@examples
#'b <- createBerryData()
#'m <- fitBinomialModel(b$Subjects, b$Events) %>% dplyr::filter(Index==10)
#'qtlFromQuantile(m, c(0.1, 0.2, 0.8, 0.9))
#'@export
qtlFromQuantile <- function(data, probs, type=4)
{
  rv <- data %>% 
          dplyr::summarise(QTL=list(stats::quantile(Value, 
                                         probs=probs,
                                         type=type))
                          ) %>% 
          tidyr::unnest() %>% 
          tibble::add_column(Quantile=probs, .before=1)
  return (rv)
}

#' Create the QTL plot
#'@param mcmcData the tibble containing the posterior MCMC samples of the metric
#'@param mcmcVar the column in \code{mcmcData} containing the sampled values of
#'the metric
#'@param groupData the tibble containing the observed group data
#'@param groupResponse the column in \code{groupData} containing the observed
#' by-group metric values
#'@param groupSize the column in the \code{groupData} containing the group sizes
#'@param groupType The column in \code{groupData} containing the group type for
#' for the observation.  If \code{NULL}, group types are not identified in the
#' plot.  Otherwise, \code{groupType} is converted to a factor and used as an
#' aesthetic.
#' @param basicTheme one of the themes provided by \code{ggplot2}, for example
#' \code{theme_light} or \code{theme_dark}.
#'@param qtls a list containing the quantiles that define the QTL thresholds.
#'If \code{NULL} or of zero length, the QTL thresholds are not shown.
#'@param postScaleFactor scales the posterior density relative to the height of 
#'the \code{groupSize} bars.  If \code{NULL}, \code{postScaleFactor} is chosen 
#'so that the mode of the posterior is twice as high as the tallest group size 
#'bar.  See Usage Notes below.  
#'@param mcmcAlpha the transparency of the fill used to indicate QTL bands
#'@param nDensity the number of bins used to estimate the density.  Must be a 
#'power of two.  Normally, this can be left unchanged, but occasionally a
#' different value is needed to correct visual artefacts related to the shading
#' of the plot
#'@param xAxisRange the range of the x-axis, as defined by a call to 
#'\code{coord_cartesian()}.  NULL uses the default axis range.
#'@return the QTL plot
#' @examples
#' b <- createBerryData() %>% 
#'        tibble::add_column(Region=as.factor(c(1,1,1,2,2,2,3,3,3)))
#' m <- fitBinomialModel(b$Subjects, b$Events) %>% dplyr::filter(Index==10)
#' createQtlPlot(groupData=b,
#'               groupResponse=ObservedResponse,
#'               groupSize=Subjects,
#'               groupType=Region,
#'               mcmcData=m,
#'               mcmcVar=Value) +
#'               ggplot2::scale_fill_manual(name="QTL",
#'                                          labels=c("Action low",
#'                                                   "Warn low",
#'                                                   "",
#'                                                   "Warn high",
#'                                                   "Action high"),
#'                                          values=c("red",
#'                                                   "yellow",
#'                                                   "white",
#'                                                   "yellow",
#'                                                   "red")) +
#'               ggplot2::scale_colour_manual(name="Region",
#'                                            values=1:3,
#'                                            labels=c("US",
#'                                                     "EMEA",
#'                                                     "Other")) +
#'               ggplot2::labs(x="Observed response",
#'                             title="Example QTL plot")
#'@section Usage Notes:
#'The plot returned by the function is deliberately crude: the x axis is 
#'unlabelled, the legends are unfomatted, and so on.  However, it can be easily 
#'customised to include more appropriate titles, legends and so forth, as shown
#'in the second example.
#'
#'The scaled kernel estimate of the density is multiplied by
#' \code{postScaleFactor} to obtain the curve shown in the plot.  The height of 
#' the scaled kernel density at its mode is 1.  Therefore, by default, 
#' \code{postScaleFactor} is set to twice the largest group size.
#'@export
createQtlPlot <- function(mcmcData,
                          mcmcVar=Value,
                          groupData,
                          groupResponse=ObservedResponse,
                          groupSize=N,
                          groupType=NULL,
                          basicTheme=ggplot2::theme_light,
                          qtls=c(0.1, 0.2, 0.8, 0.9),
                          postScaleFactor=NULL,
                          mcmcAlpha=0.2,
                          nDensity=512,
                          xAxisRange=NULL)
{
  #Validate
  if (is.null(groupData)) stop("groupData cannot be null")
  if (is.null(quote(groupResponse))) stop("groupResponse cannot be null")
  if (!tibble::is_tibble(groupData)) stop("groupData must be a tibble")
  if (is.null(mcmcData)) stop("mcmcData cannot be null")
  if (!tibble::is_tibble(mcmcData)) stop("mcmcData must either be null or a tibble.")
  if (is.null(quote(mcmcVar))) stop("mcmcVar cannot be null")
  if (!is.null(xAxisRange)) if (length(xAxisRange) != 2) stop("xAxisRange must either be NULL or a numeric vector of length 2")
  #Begin
  qResp <- dplyr::enquo(groupResponse)
  qSize <- dplyr::enquo(groupSize)
  if (!is.null(quote(groupType))) qGroup <- dplyr::enquo(groupType)
  qVar <- dplyr::enquo(mcmcVar)
  plot <- groupData %>% ggplot2::ggplot()
  #Plot the observed data
  if (is.null(qGroup)) {
    plot <- plot + ggplot2::geom_linerange(ggplot2::aes_q(x=qResp, ymin=0, ymax=qSize))
  } else {
    plot <- plot + ggplot2::geom_linerange(ggplot2::aes_q(x=qResp, ymin=0, ymax=qSize, colour=qGroup))
  }
  #Plot the posterior density.  Could use geom_density directly, but shading
  #the action limits requires manipulation...
  if (is.null(postScaleFactor)) postScaleFactor <- 2 * (groupData %>%  dplyr::summarise(Max=max(!! qSize)))$Max[1]
  if (!is.null(qtls) & length(qtls) > 0)
  {
    #Shading required
    qtlTibble <- qtlFromQuantile(mcmcData, qtls)
    d <- ggplot2::ggplot_build(mcmcData %>% 
                    ggplot2::ggplot() +
                    ggplot2::geom_density(ggplot2::aes_q(qVar, y=~..scaled..*postScaleFactor), 
                                          n=nDensity)
                      )$data[[1]]
    d <- d %>% dplyr::mutate(Area=cut(x, 
                               breaks=c(-Inf, qtlTibble$Quantile, Inf),
                               labels=1:(nrow(qtlTibble)+1)))
    plot <- plot +
              ggplot2::geom_line(data=d, ggplot2::aes(x=x, y=y))
    for (a in unique(d$Area))
      plot <- plot + ggplot2::geom_area(data=d %>% dplyr::filter(Area == a),
                               ggplot2::aes(x=x, y=y, fill=Area),
                               alpha=mcmcAlpha)
  } else {
    #No shading required
    plot <- plot + ggplot2::geom_density(data=mcmcData,
                                         ggplot2::aes_q(qVar, y=~..scaled..*postScaleFactor), 
                                         n=nDensity)
  }
  plot <- plot + 
           basicTheme() +
           ggplot2::theme(axis.ticks.y=ggplot2::element_blank(),
                          axis.text.y=ggplot2::element_blank(),
                          axis.title.y=ggplot2::element_blank())
  if (!is.null(xAxisRange)) plot <- plot + ggplot2::coord_cartesian(xlim=xAxisRange)
  return(plot)
}

#' Utility function to create the data on page 63 of Berry et al
#' @param binary if \code{FALSE} return data at the study level.  Otherwise,
#'  return data at the subject level 
#' @return if \code{binary==FALSE}, a tibble with columns Study, Subjects and
#'  Events.  Otherwise, a tibble containing columns Study, SubjectID and Event
#' @examples
#' data <- createBerryData()
#' @export
createBerryData <- function(binary=FALSE)
{
  t <- tibble::tibble(Study=1:9,
                      Subjects=c(20, 10, 16, 19, 14, 46, 10, 9, 6),
                      Events=  c(20,  4, 11, 10,  5, 36,  9, 7, 4)) %>% 
         dplyr::mutate(ObservedResponse=Events/Subjects)
  if (!binary) return (t)
  t <- t %>% 
    dplyr::group_by(Study) %>% 
    dplyr::mutate(SubjectID=list(1:Subjects)) %>% 
    tidyr::unnest() %>% 
    dplyr::mutate(Event=SubjectID <= Events, 
                  SubjectID=SubjectID + 100*Study) %>% 
    dplyr::select(Study, SubjectID, Event)
  return (t)
}

#' Function to obtain the string that defines the default JAGS model for each
#' data type
#'
#' @param label the type of model
#' @return a string defining the required model
#' @examples
#' #To fit a binary model
#' binModel <- getModelString("binary")
#' @section Usage notes:
#' \code{label} is determined case insensitively
#' @export
getModelString <- function(label=c("binary", "binomial", "normal", "poisson", "tte"))
{
  label <- match.arg(label)
  s <- ""
  s <- switch(label,
         binary="model
                 {
                   for (j in 1:m)
                   {
                     p[j] ~ dbeta(a, b)
                   }
                   for (i in 1:k)
                   {
                     r[i] ~ dbern(p[group[i]])
                   }
                   a ~ dunif(0, 10)
                   b ~ dunif(0, 10)
                 }",
         binomial="model
                   {
                      for (i in 1:k)
                      {
                         r[i] ~ dbin(p[i], n[i])
                         p[i] ~ dbeta(a, b)
                      }
                      a ~ dunif(0, 10)
                      b ~ dunif(0, 10)
                   }",
         poisson="model 
               {
                 for (i in 1:k)
                 {
                   events[i] ~ dpois(mu[i])
                   mu[i] <- lambda[i]*exposure[i]
                   lambda[i] ~ dgamma(shape, 1/scale)
                 }
                 scale ~ dgamma(1, 1)
                 shape ~ dgamma(1, 1)
               }",
         tte="model
                {
                  #Likelihood
                  for (i in 1:k)
                  {
                     logMean[i] ~ dnorm(mu[i], n[i] / tau)
                  }
                  #Prior
                  for (i in 1:k)
                  {
                     mu[i] ~ dnorm(mu0, tau)
                  }
                  mu0 ~ dnorm(0, 1e-06)
                  invTau ~ dgamma(1e-06, 1e-06)
                  #Transform
                  tau <- 1/invTau  # Because inverse gamma isn't directly supported
                }",
          normal="model
                {
                  #Likelihood
                  for (i in 1:n)
                  {
                     x[i] ~ dnorm(mu[g[i]], tau)
                  }
                  #Prior
                  for (i in 1:k)
                  {
                     mu[i] ~ dnorm(mu0, tau)
                  }
                  mu0 ~ dnorm(0, 1e-06)
                  invTau ~ dgamma(1e-06, 1e-06)
                  #Transform
                  tau <- 1/invTau  # Because inverse gamma isn't directly supported
                }")
  if (length(s) == 0) stop(paste0("Unsupported data type [", label, "]."))
  return (s)
}


#' Fit a model to binomial data
#' @param n the vector holding the number of observations
#' @param r the vector holding the number of responders
#' @param inits a list containing the initial values of the hyperparameters.
#' See Usage Notes below.
#' @param modelString The JAGS model string that defines the model to be fitted
#' @param autoRun Logical.  If \code{TRUE}, use \code{autorun.jags}.  Otherwise,
#' use \code{run.jags}.
#' @param raw logical if \code{TRUE}, the \code{mcmc} object created by
#' JAGS is returned.  Otherwise, a \code{tibble} containing the
#' concatenation of elements of the \code{mcmc.list} created by JAGS
#' is returned
#' @param ... passed to JAGS
#' @return Either the \code{mcmc} (or \code{mcmc.list}) object returned by
#' JAGS or a tibble containing the MCMC samples from the posterior
#' distribution
#' @examples
#' #Simple use
#' b <- createBerryData() 
#' m <- fitBinomialModel(b$Subjects, b$Events) %>% dplyr::filter(Index == 10)
#' #Passing parameters to run.jags()
#' inits1 <- list(a=4, b=2)
#' inits2 <- list(a=1, b=1)
#' inits3 <- list(a=2, b=10)
#' m <- fitBinomialModel(b$Subjects, b$Events,
#'                       inits=list(inits1, inits2, inits3),
#'                       thin=2)
#' @section Usage notes:
#' If \code{modelString == NULL}, the model string is obtained by calling
#' \code{getModelString("binomial")}.\cr\cr  
#' If \code{raw == FALSE}, the chain from which each observation is drawn is
#' indicated by \code{Chain} and the dataset is transformed into tidy format,
#' with the model parameter indicated by \code{Parameter}.\cr\cr  
#' The \code{inits} parameter can be used to define the number of chains created
#' my JAGS.  If a list of lists, the number of elements in the
#' outer list defines the number of chains and the elements of each sub-list
#' define the initial value for each hyperparameter.  For example, the default
#' value of \code{inits} requests two chains.  The initial values of \code{a}
#'  and \code{b} in the first chain are 2 and 4, respectively.  In the second
#' chain, the corresponding values are 9 and 1.  If the MCMC model has converged
#' and is stationary, the initial values of the hyperparameters will be
#' irrelevant.   To check for convergence, it is necessary - but not sufficient
#' - to obtain more than one chain and to use different initial values for each
#' chain.
#' @seealso \code{\link{fitPoissonModel}},  \code{\link{fitBinaryModel}},
#' \code{\link{fitTteModel}}
#' @export
fitBinomialModel <- function(n,
                             r,
                             inits=list(list("a"=2, "b"=4), list("a"=9, "b"=1)),
                             modelString=NULL,
                             autoRun=TRUE,
                             raw=FALSE,
                             ...)
{
  #Validation
  if (length(n) != length(r)) stop ("n and r must be of the same length.")
  if (!all.equal(as.integer(n), n)) stop("Non-integer values found in n.")
  if (!all.equal(as.integer(r), r)) stop("Non-integer values found in r.")
  if (any(n < 1)) stop("Non positive-definite value(s) found in n.")
  if (any(r < 0)) stop("Negative value(s) found in r.")
  if (any(r > n)) stop("More successes than trials in r and n.")

  #Initialise
  if (is.null(modelString) | length(modelString) == 0)
    modelString <- getModelString("binomial")

  #Begin
  #Append dummy, "posterior", observation
  tempData <- list()
  tempData$n <- c(n, 1)
  tempData$r <- c(r, NA)
  tempData$k <- length(tempData$n)
  if (autoRun) f <- runjags::autorun.jags
  else f <- runjags::run.jags
  results <- f(model=modelString,
               monitor=c(c("a", "b"), 
                         paste0("p[", 1:tempData$k , "]")),
               data=tempData,
               inits=inits,
               ...)
  if (raw) return (results)
  results <- runjagsResultsToTibble(results, "p[", ...)
  return (results)
}
