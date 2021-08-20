#' Check that all Python imports exist
#'
#' @return b True or False
#'
#' @examples
checkForPackages <- function() {
  return(
      reticulate::py_module_available("numpy") &
      reticulate::py_module_available("GPy") &
      reticulate::py_module_available("pygmo") &
      reticulate::py_module_available("scipy") &
      reticulate::py_module_available("cma") &
      reticulate::py_module_available("nplot") &
      reticulate::py_module_available("pyDOE2") &
      reticulate::py_module_available("numpy-stl") &
      reticulate::py_module_available("matplotlib")
  )
}

#checking for the right dimensions and legal values of bounds
test_that("test bounds", {

  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- sphere(Xtr)

  # equal rows
  testthat::equals(
    nrow(Xtr), nrow(Ytr)
  )

  if (!reticulate::virtualenv_exists())
    skip("missing python environment")

  if(!checkForPackages()) {
    skip("Packages Missing")
  }

  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- sphere(Xtr)

  # right dimension, legal values
  testthat::expect_type(
    callEshotgun(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, 0.1, pf=TRUE),
    typeof(Xtr)
  )

  # right dimension, legal values
  testthat::expect_type(
    callEshotgun(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, 0.1),
    typeof(Xtr)
  )

  # right dimension, legal values
  testthat::expect_type(
    callEshotgun(Xtr, Ytr, c(-4,-5), c(-3,1), 10L, 0.1)
    , typeof(Xtr)
  )

  #wrong dimension
  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(-5.12,-1), c(5.12, 0, 1), 10L, 0.1)
    , NULL
  )

  #wrong dimension
  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(-5.12,-1, 0), c(5.12, 0, 1), 10L, 0.1)
    , NULL
  )

  #values are not strictly smaller
  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(1,1), c(1,1), 10L, 0.1)
    , NULL
  )

  #values are not strictly smaller
  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(0,1), c(2,1), 10L, 0.1)
    , NULL
  )

  #value of lower bound is bigger than the value of the upper bound
  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(0,3), c(2,2), 10L, 0.1)
    , NULL
  )

  #value of upper bound is lower than the value of the lower bound
  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(0,0), c(2,-1), 10L, 0.1)
    , NULL
  )

  #value of upper bound is lower than the value of the lower bound
  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(0,0), c(-3,-1), 10L, 0.1)
    , NULL
  )

  Xtr <- matrix(runif(20),ncol=3)
  Ytr <- sphere(Xtr)

  testthat::expect_type(
    callEshotgun(Xtr, Ytr, c(5,0,200), c(10,10,210), 10L, 0.1)
    , typeof(Xtr)
  )

  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(5,0,200), c(10,10,190), 10L, 0.1)
    , NULL
  )

  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- sphere(Xtr)

  # right dimension, legal values
  testthat::expect_type(
    callEshotgunV2(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, 0.1, pf=TRUE),
    typeof(Xtr)
  )

  # right dimension, legal values
  testthat::expect_type(
    callEshotgunV2(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, 0.1),
    typeof(Xtr)
  )

  # right dimension, legal values
  testthat::expect_type(
    callEshotgunV2(Xtr, Ytr, c(-4,-5), c(-3,1), 10L, 0.1)
    , typeof(Xtr)
  )

  #wrong dimension
  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(-5.12,-1), c(5.12, 0, 1), 10L, 0.1)
    , NULL
  )

  #wrong dimension
  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(-5.12,-1, 0), c(5.12, 0, 1), 10L, 0.1)
    , NULL
  )

  #values are not strictly smaller
  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(1,1), c(1,1), 10L, 0.1)
    , NULL
  )

  #values are not strictly smaller
  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(0,1), c(2,1), 10L, 0.1)
    , NULL
  )

  #value of lower bound is bigger than the value of the upper bound
  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(0,3), c(2,2), 10L, 0.1)
    , NULL
  )

  #value of upper bound is lower than the value of the lower bound
  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(0,0), c(2,-1), 10L, 0.1)
    , NULL
  )

  #value of upper bound is lower than the value of the lower bound
  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(0,0), c(-3,-1), 10L, 0.1)
    , NULL
  )

  Xtr <- matrix(runif(20),ncol=3)
  Ytr <- sphere(Xtr)

  testthat::expect_type(
    callEshotgunV2(Xtr, Ytr, c(5,0,200), c(10,10,210), 10L, 0.1)
    , typeof(Xtr)
  )

  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(5,0,200), c(10,10,190), 10L, 0.1)
    , NULL
  )


})

#check for the right shape of the return value
test_that("test dimensions", {

  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- sphere(Xtr)

  # equal rows
  testthat::equals(
    nrow(Xtr), nrow(Ytr)
  )

  if (!reticulate::virtualenv_exists())
    skip("missing python environment")

  if (!checkForPackages()) {
    skip("Packages Missing")
  }

  #working
  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- sphere(Xtr)

  #working
  newX <- callEshotgun(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, 0.1)
  testthat::expect_equal(
    nrow(newX), nrow(Xtr)
  )

  #working
  newX <- callEshotgun(Xtr, Ytr, c(-4,-5), c(-3,1), 10L, 0.1)
  testthat::expect_equal(
    nrow(newX), nrow(Xtr)
  )


  Xtr <- matrix(runif(20),ncol=1)
  Ytr <- sphere(Xtr)

  #working
  newX <-callEshotgun(Xtr, Ytr, c(-5.12), c(5.12), 20L, 0.2)
  testthat::expect_equal(
    ncol(newX), ncol(Xtr)
  )

  #10 points get evaluaded
  testthat::expect_equal(
    nrow(newX), nrow(Xtr)
  )


  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- modifiedLevy(Xtr)

  #working
  newX <- callEshotgun(Xtr, Ytr, c(-5,-1), c(5,1), 10L, 0.1)
  testthat::expect_equal(
    nrow(newX),
    nrow(Xtr)
  )

  #working
  testthat::expect_equal(
    ncol(newX),
    ncol(Xtr)
  )

  # Edit the Ytr-Matrix to provoke a error
  Xtr <- matrix(runif(20),ncol=1)
  Ytr <- modifiedLevy(Xtr)[1:8]
  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(-5), c(5), 10L, 0.1),
    NULL
  )


  # Edit the Ytr-Matrix to provoke a error
  Xtr <- matrix(runif(20),ncol=1)
  Ytr <- modifiedLevy(Xtr)[1:15]
  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(-5), c(5), 10L, 0.1),
    NULL
  )


  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- modifiedBranin(Xtr)

  #working
  newX <- callEshotgun(Xtr, Ytr, c(-5,-1), c(5,1), 10L, 0.1)
  testthat::expect_equal(
    nrow(newX),
    nrow(Xtr)
  )

  #working
  testthat::expect_equal(
    ncol(newX),
    ncol(Xtr)
  )

  Xtr <- matrix(runif(20),ncol=3)
  Ytr <- modifiedBranin(Xtr)
  newX <- callEshotgun(Xtr, Ytr, c(-5,-4,-3), c(5,6,7), 10L, 0.1)

  #working
  testthat::expect_equal(
    ncol(newX),
    ncol(Xtr)
  )

  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- sphere(Xtr)

  #working
  newX <- callEshotgunV2(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, 0.1)
  testthat::expect_equal(
    nrow(newX), nrow(Xtr)
  )

  #working
  newX <- callEshotgunV2(Xtr, Ytr, c(-4,-5), c(-3,1), 10L, 0.1)
  testthat::expect_equal(
    nrow(newX), nrow(Xtr)
  )


  Xtr <- matrix(runif(20),ncol=1)
  Ytr <- sphere(Xtr)

  #working
  newX <-callEshotgunV2(Xtr, Ytr, c(-5.12), c(5.12), 20L, 0.2)
  testthat::expect_equal(
    ncol(newX), ncol(Xtr)
  )

  #10 points get evaluaded
  testthat::expect_equal(
    nrow(newX), nrow(Xtr)
  )


  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- modifiedLevy(Xtr)

  #working
  newX <- callEshotgunV2(Xtr, Ytr, c(-5,-1), c(5,1), 10L, 0.1)
  testthat::expect_equal(
    nrow(newX),
    nrow(Xtr)
  )

  #working
  testthat::expect_equal(
    ncol(newX),
    ncol(Xtr)
  )

  # Edit the Ytr-Matrix to provoke a error
  Xtr <- matrix(runif(20),ncol=1)
  Ytr <- modifiedLevy(Xtr)[1:8]
  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(-5), c(5), 10L, 0.1),
    NULL
  )


  # Edit the Ytr-Matrix to provoke a error
  Xtr <- matrix(runif(20),ncol=1)
  Ytr <- modifiedLevy(Xtr)[1:15]
  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(-5), c(5), 10L, 0.1),
    NULL
  )


  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- modifiedBranin(Xtr)

  #working
  newX <- callEshotgunV2(Xtr, Ytr, c(-5,-1), c(5,1), 10L, 0.1)
  testthat::expect_equal(
    nrow(newX),
    nrow(Xtr)
  )

  #working
  testthat::expect_equal(
    ncol(newX),
    ncol(Xtr)
  )

  Xtr <- matrix(runif(20),ncol=3)
  Ytr <- modifiedBranin(Xtr)
  newX <- callEshotgunV2(Xtr, Ytr, c(-5,-4,-3), c(5,6,7), 10L, 0.1)

  #working
  testthat::expect_equal(
    ncol(newX),
    ncol(Xtr)
  )

})

#For v1 and v2 alike:
#2 tests for epsilon lower than 0.0
#2 tests for epsilon higher than 1.0
#2 tests for epsilon equal to 0.0 and 1.0
#1 test for different dimension
test_that("test epsilon", {

  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- modifiedBranin(Xtr)

  # equal rows
  testthat::equals(
    nrow(Xtr), nrow(Ytr)
  )


  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- modifiedBranin(Xtr)

  if(!reticulate::virtualenv_exists())
    skip("missing python environment")

  if (!checkForPackages()) {
    skip("Packages Missing")
  }

  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, 2.1)
    , NULL
  )

  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, 1.1)
    , NULL
  )

  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, -2.1)
    , NULL
  )

  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, -1.1)
    , NULL
  )

  Xtr <- matrix(runif(20),ncol=3)
  Ytr <- modifiedBranin(Xtr)
  testthat::expect_type(
    callEshotgun(Xtr, Ytr, c(-5,-4,-3), c(5,6,7), 10L, 0.0),
    typeof(Xtr)
  )

  testthat::expect_type(
    callEshotgun(Xtr, Ytr, c(-5,-4,-3), c(5,6,7), 10L, 1.0),
    typeof(Xtr)
  )

  Xtr <- matrix(runif(20),ncol=3)
  Ytr <- modifiedBranin(Xtr)
  testthat::expect_equal(
    callEshotgun(Xtr, Ytr, c(-5,-4,-3), c(5,6,7), 10L, 2.0),
    NULL
  )

  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, 2.1)
    , NULL
  )

  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, 1.1)
    , NULL
  )

  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, -2.1)
    , NULL
  )

  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(-5.12,-1), c(5.12, 0), 10L, -1.1)
    , NULL
  )

  Xtr <- matrix(runif(20),ncol=3)
  Ytr <- modifiedBranin(Xtr)
  testthat::expect_type(
    callEshotgunV2(Xtr, Ytr, c(-5,-4,-3), c(5,6,7), 10L, 0.0),
    typeof(Xtr)
  )

  testthat::expect_type(
    callEshotgunV2(Xtr, Ytr, c(-5,-4,-3), c(5,6,7), 10L, 1.0),
    typeof(Xtr)
  )

  Xtr <- matrix(runif(20),ncol=3)
  Ytr <- modifiedBranin(Xtr)
  testthat::expect_equal(
    callEshotgunV2(Xtr, Ytr, c(-5,-4,-3), c(5,6,7), 10L, 2.0),
    NULL
  )
})
