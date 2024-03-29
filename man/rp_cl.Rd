% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rp.R
\name{rp_cl}
\alias{rp_cl}
\alias{crqa_cl}
\title{Fast (C)RQA (command line crp)}
\usage{
rp_cl(
  y1,
  y2 = NULL,
  emDim = 1,
  emLag = 1,
  emRad = NA,
  DLmin = 2,
  VLmin = 2,
  theiler = 0,
  win = min(length(y1), ifelse(is.null(y2), (length(y1) + 1), length(y2)), na.rm = TRUE),
  step = win,
  JRP = FALSE,
  distNorm = c("EUCLIDEAN", "MAX", "MIN", "OP")[[1]],
  standardise = c("none", "mean.sd", "median.mad")[1],
  returnMeasures = TRUE,
  returnRPvector = FALSE,
  returnLineDist = FALSE,
  doPlot = c("noplot", "rp", "distmat")[[1]],
  path_to_rp = getOption("casnet.path_to_rp"),
  saveOut = FALSE,
  path_out = NULL,
  file_ID = NULL,
  silent = TRUE,
  surrogateTest = FALSE,
  targetValue = 0.05,
  useParallel = FALSE,
  ...
)
}
\arguments{
\item{y1}{Time series 1}

\item{y2}{Time series 2 for Cross Recurrence Analysis (default = \code{NULL})}

\item{emDim}{Embedding dimensions (default = \code{1})}

\item{emLag}{Embedding lag (default = \code{1})}

\item{emRad}{Radius on distance matrix (default = \code{1})}

\item{DLmin}{Minimum length of diagonal structure to be considered a line (default = \code{2})}

\item{VLmin}{Minimum length of vertical structure to be considered a line (default = \code{2})}

\item{theiler}{Theiler window (default = \code{0})}

\item{win}{Window to calculate the (C)RQA (default = minimum of length of \code{y1} or \code{y2})}

\item{step}{Stepsize for sliding windows (default = size of \code{win}, so no sliding window)}

\item{JRP}{Wether to calculate a Joint Recurrence Plot (default = \code{FALSE})}

\item{distNorm}{One of "EUCLIDEAN" (default), \verb{"MAX", "MIN"}, or \code{"OP"} for an Order Pattern recurrence matrix}

\item{standardise}{Standardise data: \code{"none"} (default), \code{"mean.sd"}, or \code{"median.mad"}}

\item{returnMeasures}{Return the (C)RQA measures? (default = \code{TRUE})}

\item{returnRPvector}{Return the recurrent points in a dataframe? (default = \code{FALSE})}

\item{returnLineDist}{Return the distribution of diagonal and horizontal line length distances (default = \code{FALSE})}

\item{doPlot}{Produce a plot of the recurrence matrix by calling \code{\link[=rp_plot]{rp_plot()}}, values can be \code{"rp"} (the thresholded recurrence matrix),\code{"distmat"} (the unthresholded recurrence matrix) or \code{"noplot"} (default = \code{"noplot"})}

\item{path_to_rp}{Path to the command line executable (default = path set during installation, use \code{getOption("casnet.path_to_rp")} to see)}

\item{saveOut}{Save the output to files? If \code{TRUE} and \code{path_out = NA}, the current working directory will be used (default = \code{FALSE})}

\item{path_out}{Path to save output if \code{saveOut = TRUE} (default = \code{NULL})}

\item{file_ID}{A file ID which will be a prefix to to the filename if \code{saveOut = TRUE} (default = \code{NULL}, an integer will be added tot the file name to ensure unique files)}

\item{silent}{Do not display any messages (default = \code{TRUE})}

\item{surrogateTest}{Perform surrogate tests. If \code{TRUE}, will run surrogate tests using default settings for a two-sided test of \eqn{H_0: The data generating process is a rescaled linear Gaussian process} at \eqn{\alpha = .05} (arguments \verb{ns = 39, fft = TRUE, amplitude = TRUE})}

\item{targetValue}{A value passed to \code{est_radius(...,type="fixed", targetMeasure="RR")} if \code{is.na(emRad)==TRUE}. This is useful for windowed analysis, it will estimate a new radius for each window.}

\item{useParallel}{Speed up calculations by using the parallel processing options provided by \code{parallel} to assign a seperate process/core for each window in windowed (C)RQA analysis and \code{\link[parallel:detectCores]{parallel::detectCores()}} with\code{logical = TRUE} to decide on the available cores (default = \code{FALSE})}

\item{...}{Additional parameters (currently not used)}
}
\value{
A list object containing 1-3 elements, depending on arguments requesting output.
\itemize{
\item \code{rqa_measures} - A list of the (C)RQA measures returned if \code{returnMeasures = TRUE}:
\enumerate{
\item RR = Recurrence rate
\item DET = Determinism
\item DET_RR = Ratio DET/RR
\item LAM = Laminarity
\item LAM_DET = Ratio LAM/DET
\item L_max = maximal diagonal line length
\item L_mean = mean diagonal line length
\item L_entr = Entropy of diagonal line length distribution
\item DIV =  Divergence (1/L_max)
\item V_max = maximal vertical line length
\item TT = Trapping time
\item V_entr = Entropy of vertical line length distribution
\item T1 = Recurrence times 1st type
\item T2 = Recurrence times 2nd type
\item W_max = Max interval length
\item W_mean = Mean of interval lengths
\item W_entr = Entropy of interval length distribution
\item W_prob = Probability of interval
\item F_min = F min
}
\item \code{rqa_rpvector} - The radius thresholded distance matrix (recurrence matrix), which can be visualised as a recurrence plot by calling \code{\link[=rp_plot]{rp_plot()}}. If a sliding window analysis is conducted this will be a list of matrices and could potentially grow too large to handle. It is recommended you save the output to disk by setting \code{saveOut = TRUE}.
\item \code{rqa_diagdist} - The distribution of diagonal line lengths
}
}
\description{
This function will run the \href{http://tocsy.pik-potsdam.de/commandline-rp.php}{commandline Recurrence Plots} executable provided by Norbert Marwan.
}
\details{
The \code{rp} executable is installed when the function is called for the first time and is renamed to \code{rp}, from a platform specific filename downloaded from http://tocsy.pik-potsdam.de/commandline-rp.php or extracted from an archive located in the directory:
\verb{...\\\\casnet\\\\commandline_rp\\\\}
The file is copied to the directory: \verb{...\\\\casnet\\\\exec\\\\}
The latter location is stored as an option and can be read by calling \code{getOption("casnet.path_to_rp")}.
}
\note{
The platform specific \code{rp} command line executables were created by Norbert Marwan and obtained under a Creative Commons License from the website of the Potsdam Institute for Climate Impact Research at \url{http://tocsy.pik-potsdam.de/}.

The full copyright statement on the website is as follows:

(C) 2004-2017 SOME RIGHTS RESERVED

University of Potsdam, Interdisciplinary Center for Dynamics of Complex Systems, Germany

Potsdam Institute for Climate Impact Research, Transdisciplinary Concepts and Methods, Germany

This work is licensed under a \href{https://creativecommons.org/licenses/by-nc-nd/2.0/de/}{Creative Commons Attribution-NonCommercial-NoDerivs 2.0 Germany License}.

More information about recurrence analysis can be found on the \href{http://www.recurrence-plot.tk}{Recurrence Plot} website.
}
\section{Troubleshooting}{

Some notes on resolving errors with \code{rp}.The script will first try to download the correct executable, if that fails, the copy will have... failed. It should be relatively easy to get \code{rp_cl()} working though, by using some custom settings:
\itemize{
\item \emph{Copy failed} - Every time the function \code{rp_cl()} is called it will check whether a log file \code{rp_instal_log.txt} is present in the \verb{...\\\\casnet\\\\exec\\\\} directory. If you delete the \code{rp_instal_log.txt} file, and call the function, another attempt will be made to download a copy of the executable.
\item \emph{Copy still fails and/or no permission to copy} - If you cannot acces the directory \verb{...\\\\casnet\\\\commandline_rp\\\\}, download the appropriate executable from the \href{http://tocsy.pik-potsdam.de/commandline-rp.php}{Commandline Recurrence Plots} page and copy to a directory you do have the rights to: \emph{execute} commands, \emph{write} and \emph{read} files. Make sure you rename the file to \code{rp} (\code{rp.exe} on Windows OS). Then, either pass the path to \code{rp} as the argument \code{path_to_rp} in the \code{rp_cl(.., path_to_rp = "YOUR_PATH")} function call, or, as a more permanent solution, set the \code{path_to_rp} option by calling \code{options(casnet.path_to_rp="YOUR_PATH")}.
\item \emph{Error in execution of \code{rp}} - This can have a variety of causes, the \code{rp} executable is called using \code{\link[=system2]{system2()}} and makes use of the \code{\link[=normalizePath]{normalizePath()}} function with argument \code{mustWork = FALSE}. Problems caused by specific OS, machine, or, locale problems (e.g. the \code{winslash} can be reported as an \href{https://github.com/FredHasselman/casnet/issues}{issue on Github}). One execution error that occurs when the OS is not recognised properly can be resolved by chekcing \code{getOption("casnet.rp_prefix")}. On Windows OS this should return an empty character vector, on Linux or macOS it should return \code{"./"}. You can manually set the correct prefix by calling \code{options(casnet.rp_prefix="CORRECT OS PREFIX")} and fill in the prefix that is correct for your OS
}
}

\seealso{
Other Recurrence Quantification Analysis: 
\code{\link{rp_measures}()},
\code{\link{rp_measures_main}()}
}
\concept{Recurrence Quantification Analysis}
