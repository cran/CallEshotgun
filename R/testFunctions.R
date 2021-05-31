#' sphere
#'
#' ##########################################################################
#'#
#'# SPHERE FUNCTION, MODIFIED
#'#
#'# Authors: Sonja Surjanovic, Simon Fraser University
#'#          Derek Bingham, Simon Fraser University
#'# Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
#'#
#'# Copyright 2013. Derek Bingham, Simon Fraser University.
#'#
#'# THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
#'# FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
#'# derivative works, such modified software should be clearly marked.
#'# Additionally, this program is free software; you can redistribute it
#'# and/or modify it under the terms of the GNU General Public License as
#'# published by the Free Software Foundation; version 2.0 of the License.
#'# Accordingly, this program is distributed in the hope that it will be
#'# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#'# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#'# General Public License for more details.
#'#
#'# For function details and reference information, see:
#'# http://www.sfu.ca/~ssurjano/
#'#
#'##########################################################################
#'#
#'# INPUT:
#'#
#'# xx = c(x1, x2, x3, x4, x5, x6)
#'#
#'##########################################################################
#'
#' @param x Vector x = c(x1, x2, x3, x4, x5, x6)
#'
#' @return a matrix with the sphere function applied to each row
#' @export
#'
#' @examples
sphere <- function (x) {
  matrix(apply(x, # matrix
               1, # margin (apply over rows)
               function(x) {
                 sum(x ^ 2)  # objective function
               }),
         , 1) # number of columns
}


#' egg
#'
#'##########################################################################
#'#
#'# EGGHOLDER FUNCTION
#'#
#'# Authors: Sonja Surjanovic, Simon Fraser University
#'#          Derek Bingham, Simon Fraser University
#'# Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
#'#
#'# Copyright 2013. Derek Bingham, Simon Fraser University.
#'#
#'# THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
#'# FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
#'# derivative works, such modified software should be clearly marked.
#'# Additionally, this program is free software; you can redistribute it
#'# and/or modify it under the terms of the GNU General Public License as
#'# published by the Free Software Foundation; version 2.0 of the License.
#'# Accordingly, this program is distributed in the hope that it will be
#'# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#'# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#'# General Public License for more details.
#'#
#'# For function details and reference information, see:
#'# http://www.sfu.ca/~ssurjano/
#'#
#'##########################################################################
#'#
#'# INPUT:
#'#
#'# xx = c(x1, x2)
#'#
#'##########################################################################
#'
#'
#' @param xx Vector with two elements xx = c(x1, x2)
#'
#' @return y an single double
#' @export
#'
#' @examples
egg <- function(xx) {
  x1 <- xx[1]
  x2 <- xx[2]

  term1 <- -(x2+47) * sin(sqrt(abs(x2+x1/2+47)))
  term2 <- -x1 * sin(sqrt(abs(x1-(x2+47))))

  y <- term1 + term2
  return(y)
}


#' modifiedEgg
#'
#' The Eggholder Function applied to a matrix
#'
#' @param x matrix with the inital points
#'
#' @return a matrix
#' @export
#'
#' @examples
modifiedEgg <- function(x) {
  matrix(apply(x, # matrix
               1, # margin (apply over rows)
               function(x){
                 egg(x)
               }),
         , 1) # number of columns
}


#' schwef
#'
#'###########################################################################
#'#
#'# SCHWEFEL FUNCTION
#'#
#'# Authors: Sonja Surjanovic, Simon Fraser University
#'#          Derek Bingham, Simon Fraser University
#'# Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
#'#
#'# Copyright 2013. Derek Bingham, Simon Fraser University.
#'#
#'# THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
#'# FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
#'# derivative works, such modified software should be clearly marked.
#'# Additionally, this program is free software; you can redistribute it
#'# and/or modify it under the terms of the GNU General Public License as
#'# published by the Free Software Foundation; version 2.0 of the License.
#'# Accordingly, this program is distributed in the hope that it will be
#'# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#'# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#'# General Public License for more details.
#'#
#'# For function details and reference information, see:
#'# http://www.sfu.ca/~ssurjano/
#'#
#'##########################################################################
#'#
#'# INPUT:
#'#
#'# xx = c(x1, x2, ..., xd)
#'#
#'##########################################################################
#'
#' @param xx vector
#'
#' @return a double
#' @export
#'
#' @examples
schwef <- function(xx) {
  d <- length(xx)

  sum <- sum(xx*sin(sqrt(abs(xx))))

  y <- 418.9829*d - sum
  return(y)
}


#' modifiedschwef
#'
#' The schwefel Function applied to a matrix
#'
#' @param x matrix with the initial points
#'
#' @return a matrix
#' @export
#'
#' @examples
modifiedschwef <- function(x) {
  matrix(apply(x, # matrix
               1, # margin (apply over rows)
               function(x){
                 schwef(x)
               }),
         , 1) # number of columns
}


#' levy
#'
#'##########################################################################
#'#
#'# LEVY FUNCTION
#'#
#'# Authors: Sonja Surjanovic, Simon Fraser University
#'#          Derek Bingham, Simon Fraser University
#'# Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
#'#
#'# Copyright 2013. Derek Bingham, Simon Fraser University.
#'#
#'# THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
#'# FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
#'# derivative works, such modified software should be clearly marked.
#'# Additionally, this program is free software; you can redistribute it
#'# and/or modify it under the terms of the GNU General Public License as
#'# published by the Free Software Foundation; version 2.0 of the License.
#'# Accordingly, this program is distributed in the hope that it will be
#'# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#'# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#'# General Public License for more details.
#'#
#'# For function details and reference information, see:
#'# http://www.sfu.ca/~ssurjano/
#'#
#'##########################################################################
#'#
#'# INPUT:
#'#
#'# xx = c(x1, x2, ..., xd)
#'#
#'##########################################################################
#'
#' @param xx vector
#'
#' @return a dounle
#' @export
#'
#' @examples
levy <- function(xx) {
  d <- length(xx)
  w <- 1 + (xx - 1)/4

  term1 <- (sin(pi*w[1]))^2
  term3 <- (w[d]-1)^2 * (1+1*(sin(2*pi*w[d]))^2)

  wi <- w[1:(d-1)]
  sum <- sum((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))

  y <- term1 + sum + term3
  return(y)
}


#' modifiedlevy
#'
#' The Levy Funciton applied to a matrix
#'
#' @param x vector
#'
#' @return a matrix
#' @export
#'
#' @examples
modifiedLevy <- function(x) {
  matrix(apply(x, # matrix
               1, # margin (apply over rows)
               function(x){
                 levy(x)
               }),
         , 1) # number of columns
}


#' branin
#'
#'##########################################################################
#'# BRANIN FUNCTION
#'#
#'# Authors: Sonja Surjanovic, Simon Fraser University
#'#          Derek Bingham, Simon Fraser University
#'# Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
#'#
#'# Copyright 2013. Derek Bingham, Simon Fraser University.
#'#
#'# THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
#'# FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
#'# derivative works, such modified software should be clearly marked.
#'# Additionally, this program is free software; you can redistribute it
#'# and/or modify it under the terms of the GNU General Public License as
#'# published by the Free Software Foundation; version 2.0 of the License.
#'# Accordingly, this program is distributed in the hope that it will be
#'# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#'# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#'# General Public License for more details.
#'#
#'# For function details and reference information, see:
#'# http://www.sfu.ca/~ssurjano/
#'#
#'##########################################################################
#'#
#'# INPUTS:
#'#
#'# xx = c(x1, x2)
#'# a = constant (optional), with default value 1
#'# b = constant (optional), with default value 5.1/(4*pi^2)
#'# c = constant (optional), with default value 5/pi
#'# r = constant (optional), with default value 6
#'# s = constant (optional), with default value 10
#'# t = constant (optional), with default value 1/(8*pi)
#'#
#'##########################################################################
#'
#'
#' @param xx vector
#' @param a double
#' @param b double
#' @param c double
#' @param r double
#' @param s double
#' @param t double
#'
#' @return a one single double
#' @export
#'
#' @examples
branin <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi)) {

  x1 <- xx[1]
  x2 <- xx[2]

  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)

  y <- term1 + term2 + s
  return(y)
}

#' modifiedBranin
#'
#' The Branin Function applied to a matrix
#'
#' @param x matrix with the initial points
#'
#' @return m a matrix
#' @export
#'
#' @examples
modifiedBranin <- function(x) {
  matrix(apply(x, # matrix
               1, # margin (apply over rows)
               function(x){
                 branin(x)
               }),
         , 1) # number of columns
}
