%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{PairedPSCBS}
\docType{class}
\alias{PairedPSCBS}


 \title{The PairedPSCBS class}

 \description{
  Package:  PSCBS \cr
\bold{Class PairedPSCBS}\cr

\code{list}\cr
\code{~~|}\cr
\code{~~+--}\code{PairedPSCBS}\cr

\bold{Directly known subclasses:}\cr
\cr

public class \bold{PairedPSCBS}\cr
extends list\cr



  A PairedPSCBS is an object containing the results from the
  Paired PSCBS method.
 }

 \usage{PairedPSCBS(fit=list(), ...)}

 \arguments{
   \item{fit}{A \code{\link[base]{list}} structure containing the Paired PSCBS results.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.data.frame} \tab  -\cr
 \tab \code{\link[PSCBS:callAB.PairedPSCBS]{callAB}} \tab Calls segments that are in allelic balance.\cr
 \tab \code{\link[PSCBS:callLOH.PairedPSCBS]{callLOH}} \tab Calls segments that are in LOH.\cr
 \tab \code{callSegments} \tab  -\cr
 \tab \code{\link[PSCBS:estimateDeltaAB.PairedPSCBS]{estimateDeltaAB}} \tab Estimate a threshold for calling allelic balance from DH.\cr
 \tab \code{\link[PSCBS:estimateDeltaLOH.PairedPSCBS]{estimateDeltaLOH}} \tab Estimate a threshold for calling LOH from DH.\cr
 \tab \code{\link[PSCBS:estimateKappa.PairedPSCBS]{estimateKappa}} \tab Estimate global background in segmented copy numbers.\cr
 \tab \code{extractC1C2} \tab  -\cr
 \tab \code{extractTCNAndDHs} \tab  -\cr
 \tab \code{\link[PSCBS:plotTracks.PairedPSCBS]{plotTracks}} \tab Plots parental specific copy numbers along the genome.\cr
 \tab \code{plotTracks2} \tab  -\cr
 \tab \code{plotTracksManyChromosomes} \tab  -\cr
}


 \bold{Methods inherited from list}:\cr
all.equal, as.data.frame, attachLocally, callHooks, relist, within


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   The \code{\link{segmentByPairedPSCBS}}() method returns an object of this class.
 }
\keyword{classes}