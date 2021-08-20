#' Check for required Library and Python Packages
#'
#' The function will download all required Python Packages that
#' are rquired for the e-shotgun to run properly.
#'
#' @import reticulate
#'
#' @param method method of installation
#' @param conda environment
#'
#' @export
#' @examples checkLibraries()
checkLibraries <- function(method = "auto", conda = "auto") {

  # check if a python environment exists
  pythonEnv <- reticulate::virtualenv_exists()

  if(pythonEnv) {
    # Conda doesn't include all needed imports
    reticulate::py_install("numpy", method = method, conda = conda, pip=TRUE)
    reticulate::py_install("GPy", method = method, conda = conda, pip=TRUE)
    reticulate::py_install("pygmo",method = method, conda = conda)
    reticulate::py_install("scipy", method = method, conda = conda, pip=TRUE)
    reticulate::py_install("cma", method = method, conda = conda, pip=TRUE)
    reticulate::py_install("nlopt", method = method, conda = conda, pip=TRUE)
    reticulate::py_install("pyDOE2", method = method, conda = conda, pip=TRUE)
    reticulate::py_install("numpy-stl", method = method, conda = conda, pip=TRUE)
    reticulate::py_install("matplotlib", method = method, conda = conda, pip=TRUE)
  }else {
    print("Python environment is missing")
  }
}



#' Call the e-shotgun Version 1 with Pareto front selection function in python
#'
#' The function checks the passed parameter and than calls the e-shotgun Python
#' implementation and returns a matrix with the evaluated points.
#'
#' George De Ath, Richard M. Everson, Jonathan E. Fieldsend, and Alma A. M. Rahat. 2020.
#' e-shotgun : e-greedy Batch Bayesian Optimisation. In Genetic and Evolutionary Computation Conference (GECCO ’20), July 8–12, 2020, Cancún, Mexico.
#' ACM, New York, NY, USA, 9 pages.
#' https://doi.org/10.1145/3377930.3390154
#' https://github.com/georgedeath/eshotgun
#'
#' @param Xtr a matrix containing the initial points
#' @param Ytr a matrix containing the evaluation of Xtr with a given function
#' @param f_lb a vector with the values of the lower bounds
#' @param f_ub a vector with the values of the upper bounds
#' @param q the amount if points that the e-shotgun should evaluate
#' @param epsilon the epsilon value of the e-shotgun
#' @param pf boolean that decides if pareto front is used
#'
#' @return a matrix or a vector
#' @export
#'
#' @examples
callEshotgun <- function(Xtr, Ytr, f_lb, f_ub, q=10L, epsilon=0.1, pf=FALSE) {
  py_run_file(system.file("EshotgunPy.py", package="CallEshotgun"))
  np <- import("numpy", convert = FALSE, delay_load = TRUE)

  Xnew <- tryCatch({
    errorMsg <- "Unkown"
    xrow <- nrow(Xtr)
    xcol <- ncol(Xtr)

    yrow <- nrow(Ytr)

    if(is.null(yrow)) {
      yrow <- length(Ytr)
    }

    dimLb <- length(f_lb)
    dimUb <- length(f_ub)


    # check for equal dimensions for Xtr and Ytr
    if(xrow != yrow) {
      errorMsg <- paste("Xtr and Ytr have unequal rows\nXtr is "
                        , xrow, " and Ytr is ", yrow, sep = "")
      stop()
    }

    #check lower and upper bound for a fitting dimension
    if(!(dimLb == dimUb && dimLb == xcol)) {

      if(dimLb != dimUb) {
        errorMsg <- paste("Dimension of bounds don't match!",
                      paste("\nDimension of Lower bound: ", dimLb, sep=""),
                      paste("\nDimension of Upper bound: ", dimUb, sep=""), sep="")
      }else {
        errorMsg <- paste("Dimension of bounds and Xtr don't match!",
                      paste("\nDimension of Lower bound: ", dimLb, sep=""),
                      paste("\nDimension of Xtr: ", xcol, sep=""), sep="")
      }

      stop()
    }

    #check that every number at the corresponding index of lower bound is
    #smaller than the upper bound
    #All values of the lower bound have to be strictly smaller
    if(any(f_lb > f_ub)){
      errorMsg <- "All Values of the lower bound have to be strictly smaller."

      stop()
    }


    #check epsilon between 0.0 and 1.0
    if(!(epsilon >= 0.0 && epsilon <= 1.0)) {
      errorMsg <- paste("Epsilon has to be between 0.0 and 1.0\n","Passed Epsilon is ",
                        epsilon, sep ="")
      stop()
    }

    # if the column is 1 choose special case
    if(xcol >= 2) {
      py$callShotgun(np$array(Xtr), np$array(Ytr), np$array(f_lb), np$array(f_ub), q, epsilon, pf)
    }else {
      py$callShotgun(np$array(Xtr), np$array(Ytr), f_lb, f_ub, q, epsilon, pf)
    }

  }, error = function(e) {
    cat(paste("Error:\n", errorMsg, "\n", sep=""))
  })

  return(Xnew)
}



#' Call the e-shotgun version 2 with random selection function in python
#'
#' The function checks the passed parameter and than calls the e-shotgun Python
#' implementation and returns a matrix with the evaluated points.
#'
#' George De Ath, Richard M. Everson, Jonathan E. Fieldsend, and Alma A. M. Rahat. 2020.
#' e-shotgun : e-greedy Batch Bayesian Optimisation. In Genetic and Evolutionary Computation Conference (GECCO ’20), July 8–12, 2020, Cancún, Mexico.
#' ACM, New York, NY, USA, 9 pages.
#' https://doi.org/10.1145/3377930.3390154
#' https://github.com/georgedeath/eshotgun
#'
#' @param Xtr a matrix containing the initial points
#' @param Ytr a matrix containing the evaluation of Xtr with a given function
#' @param f_lb a vector with the values of the lower bounds
#' @param f_ub a vector with the values of the upper bounds
#' @param q the amount if points that the e-shotgun should evaluate
#' @param epsilon the epsilon value of the e-shotgun
#' @param pf boolean that decides if pareto front is used
#'
#' @return a matrix or a vector
#' @export
#'
#' @examples
callEshotgunV2 <- function(Xtr, Ytr, f_lb, f_ub, q=10L, epsilon=0.1, pf=FALSE) {
  py_run_file(system.file("EshotgunPy.py", package="CallEshotgun"))
  np <- import("numpy", convert = FALSE, delay_load = TRUE)

  Xnew <- tryCatch({
    errorMsg <- "Unkown"
    xrow <- nrow(Xtr)
    xcol <- ncol(Xtr)

    yrow <- nrow(Ytr)

    if(is.null(yrow)) {
      yrow <- length(Ytr)
    }

    dimLb <- length(f_lb)
    dimUb <- length(f_ub)


    # check for equal dimensions for Xtr and Ytr
    if(xrow != yrow) {
      errorMsg <- paste("Xtr and Ytr have unequal rows\nXtr is "
                        , xrow, " and Ytr is ", yrow, sep = "")
      stop()
    }

    #check lower and upper bound for a fitting dimension
    if(!(dimLb == dimUb && dimLb == xcol)) {

      if(dimLb != dimUb) {
        errorMsg <- paste("Dimension of bounds don't match!",
                          paste("\nDimension of Lower bound: ", dimLb, sep=""),
                          paste("\nDimension of Upper bound: ", dimUb, sep=""), sep="")
      }else {
        errorMsg <- paste("Dimension of bounds and Xtr don't match!",
                          paste("\nDimension of Lower bound: ", dimLb, sep=""),
                          paste("\nDimension of Xtr: ", xcol, sep=""), sep="")
      }

      stop()
    }

    #check that every number at the corresponding index of lower bound is
    #smaller than the upper bound
    #All values of the lower bound have to be strictly smaller
    if(any(f_lb > f_ub)){
      errorMsg <- "All Values of the lower bound have to be strictly smaller."

      stop()
    }


    #check epsilon between 0.0 and 1.0
    if(!(epsilon >= 0.0 && epsilon <= 1.0)) {
      errorMsg <- paste("Epsilon has to be between 0.0 and 1.0\n","Passed Epsilon is ",
                        epsilon, sep ="")
      stop()
    }

    # if the column is 1 choose special case
    if(xcol >= 2) {
      py$callShotgunV2(np$array(Xtr), np$array(Ytr), np$array(f_lb), np$array(f_ub), q, epsilon, pf)
    }else {
      py$callShotgunV2(np$array(Xtr), np$array(Ytr), f_lb, f_ub, q, epsilon, pf)
    }

  }, error = function(e) {
    cat(paste("Error:\n", errorMsg, "\n", sep=""))
  })

  return(Xnew)
}
