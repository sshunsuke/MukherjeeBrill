% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MukherjeeBrill_facade.R
\name{util_MB_Colebrook}
\alias{util_MB_Colebrook}
\title{Colebrook correlation for the Darcy friction factor}
\usage{
util_MB_Colebrook(Re, roughness, D, tol = 1e-08, itMax = 10, warn = TRUE)
}
\arguments{
\item{Re}{Reynold number}

\item{roughness}{Pipe roughness}

\item{D}{Pipe diameter}

\item{tol}{Tolerance in Newton-Raphson method (optional)}

\item{itMax}{Maximum number of iteration  (optional)}

\item{warn}{If FALSE, not show warnings when Re <= 4000  (optional)}
}
\value{
Darcy friction factor
}
\description{
Calculate the Darcy friction factor with the Colebrook correlation.
As the correlation cannot be resolved explicitly, Newton-Raphson method is used.
}
