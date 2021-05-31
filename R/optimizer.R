#' runSampleOpt
#'
#' Optimizing a function with the e-shotgun
#'
#' @param fn function of the testproblem
#' @param budget budget for the run
#'
#' @return
#' @export
#'
#' @examples
runSampleOpt <- function(fn, budget = 100){
  initBudget <- budget
  Xtr <- matrix(runif(20),ncol=2)
  Ytr <- fn(Xtr)

  budget <- initBudget - nrow(Xtr)
  while(budget > 0){
    newX <- callEshotgun(Xtr, Ytr, c(-5,-5), c(5,5))

    newY <- fn(newX)
    Xtr <- rbind(Xtr, newX)

    Ytr <- rbind(Ytr, newY)

    budget <- initBudget - nrow(Xtr)
  }

  print(Xtr)
  print(Ytr)
}

#runSampleOpt(sphere)
#runSampleOpt(modifiedBranin)
#runsampleOpt(modifiedEgg)
#runSampleOpt(modifiedLevy)
#runSampleOpt(modifiedschwef)
