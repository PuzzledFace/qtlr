## R CMD check results
There were no ERRORs or NOTEs. 

There was one WARNING:
WARNING

*  ‘qpdf’ is needed for checks on size reduction of PDFs

  I am unable to install qpdf in my corporate environment.

There was 1 NOTE:

* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'R6'

  R6 is a build-time dependency.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Examples
Several examples run for more than 5 seconds elapsed time.  The underlying 
methodology is MCMC, so this seems unavoidable.  I have copied a typical summary
table below.

   Examples with CPU or elapsed time > 5s
                      user system elapsed
   fitBinomialModel 15.661  0.544  16.206
   fitBinaryModel   15.112  0.420  15.532
   createQtlPlot     6.786  0.355   7.149
   qtlFromQuantile   5.150  0.233   5.383
   fitTteModel       4.807  0.329   5.136
   
## Other comments
This is my first CRAN submission.
