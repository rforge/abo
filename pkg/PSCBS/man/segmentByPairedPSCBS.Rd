%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentByPairedPSCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{segmentByPairedPSCBS}
\alias{segmentByPairedPSCBS.default}
\alias{segmentByPairedPSCBS}
 \alias{segmentByPairedPSCBS}

 \title{Segment total copy numbers and allele B fractions using the Paired PSCBS method}

 \description{
  Segment total copy numbers and allele B fractions using the Paired PSCBS method [1].
  This method requires matched normals.
  This is a low-level segmentation method.
  It is intended to be applied to one tumor-normal sample at the time.
 }

 \usage{\method{segmentByPairedPSCBS}{default}(CT, betaT, betaN, muN=NULL, chromosome=0, x=NULL, alphaTCN=0.009, alphaDH=0.001, undoTCN=Inf, undoDH=Inf, ..., flavor=c("tcn&dh", "tcn,dh", "sqrt(tcn),dh", "sqrt(tcn)&dh"), tbn=TRUE, joinSegments=TRUE, knownCPs=NULL, seed=NULL, verbose=FALSE)}

 \arguments{
   \item{CT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor total copy numbers to be segmented.}
   \item{betaT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor allele B fractions (BAFs).}
   \item{betaN}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J matched normal BAFs.}
   \item{muN}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genotype calls in
        \{0,1/2,1\} for AA, AB, and BB, respectively. If not given,
        they are estimated from the normal BAFs using
        \code{\link[aroma.light]{callNaiveGenotypes}} as described in [2].}
   \item{chromosome}{(Optional) An \code{\link[base]{integer}} scalar
       (or a \code{\link[base]{vector}} of length J contain a unique value).
       Only used for annotation purposes.}
   \item{x}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic locations.
            If \code{\link[base]{NULL}}, index locations \code{1:J} are used.}
   \item{alphaTCN, alphaDH}{The significance levels for segmenting total
        copy numbers (TCNs) and decrease-in-heterozygosity signals (DHs),
        respectively.}
   \item{undoTCN, undoDH}{Non-negative \code{\link[base]{numeric}}s.  If less than +\code{\link[base:is.finite]{Inf}},
        then a cleanup of segmentions post segmentation is done.
        See argument \code{undo} of \code{\link{segmentByCBS}}() for more
        details.}
   \item{...}{Additional arguments passed to \code{\link{segmentByCBS}}().}
   \item{flavor}{A \code{\link[base]{character}} specifying what type of segmentation and
     calling algorithm to be used.}
   \item{tbn}{If \code{\link[base:logical]{TRUE}}, \code{betaT} is normalized before segmentation
     using the TumorBoost method [2], otherwise not.}
   \item{joinSegments}{If \code{\link[base:logical]{TRUE}}, there are no gaps between neighboring
     segments.
     If \code{\link[base:logical]{FALSE}}, the boundaries of a segment are defined by the support
     that the loci in the segments provides, i.e. there exist a locus
     at each end point of each segment.  This also means that there
     is a gap between any neighboring segments, unless the change point
     is in the middle of multiple loci with the same position.
     The latter is what \code{DNAcopy::segment()} returns.
   }
   \item{knownCPs}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of known
     change point locations.}
   \item{seed}{An (optional) \code{\link[base]{integer}} specifying the random seed to be
     set before calling the segmentation method.  The random seed is
     set to its original state when exiting.  If \code{\link[base]{NULL}}, it is not set.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns the segmentation results as a \code{\link{PairedPSCBS}} object.
 }

 \details{
   Internally \code{\link{segmentByCBS}}() is used for segmentation.

   This segmentation method does \emph{not} support weights.

   The "DNAcopy::segment" implementation of CBS uses approximation
   through random sampling for some estimates.  Because of this,
   repeated calls using the same signals may result in slightly
   different results, unless the random seed is set/fixed.
 }

 \section{Missing and non-finite values}{
   The total copy number signals as well as any optional positions
   must not contain missing values, i.e. \code{\link[base]{NA}}s or \code{\link[base:is.finite]{NaN}}s.
   If there are any, an informative error is thrown.
   Allele B fractions may contain missing values, because such are
   interpreted as representing non-polymorphic loci.

   None of the input signals may have infinite values, i.e. -\code{\link[base:is.finite]{Inf}} or \code{\link[base:is.finite]{Inf}}.
   If so, an informative error is thrown.
 }

 \examples{

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# (note to package developers: this example data set may
#  be replaced in a future release of the package)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS")
data <- R.utils::loadObject(pathname)

# Order
o <- order(data$position)
data <- data[o,]
str(data)
R.oo::attachLocally(data)
x <- position


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
CTs <- dropSegmentationOutliers(CT, chromosome=1, x=x, verbose=-10)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x,
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
plotTracks(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == 12)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bootstrap segment level estimates
# (used by the AB caller, which, if skipped here,
#  will do it automatically)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- bootstrapTCNandDHByRegion(fit, verbose=-10)
print(fit)
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in allelic balance (AB)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explicitly estimate the threshold in DH for calling AB
# (which be done by default by the caller, if skipped here)
deltaAB <- estimateDeltaAB(fit, flavor="qq(DH)", verbose=-10)
print(deltaAB)
## [1] 0.1657131

fit <- callAB(fit, delta=deltaAB, verbose=-10)
print(fit)
plotTracks(fit)

# Even if not explicitly specified, the estimated
# threshold parameter is returned by the caller
stopifnot(fit$params$deltaAB == deltaAB)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in loss-of-heterozygosity (LOH)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explicitly estimate the threshold in C1 for calling LOH
# (which be done by default by the caller, if skipped here)
deltaLOH <- estimateDeltaLOH(fit, flavor="minC1|nonAB", verbose=-10)
print(deltaLOH)
## [1] 0.625175

fit <- callLOH(fit, delta=deltaLOH, verbose=-10)
print(fit)
plotTracks(fit)

# Even if not explicitly specified, the estimated
# threshold parameter is returned by the caller
stopifnot(fit$params$deltaLOH == deltaLOH)
}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   Internally, \code{\link[aroma.light]{callNaiveGenotypes}} is used to
   call naive genotypes, \code{\link[aroma.light]{normalizeTumorBoost}} is
   used for TumorBoost normalization, and \code{\link{segmentByCBS}}() is used
   to segment TCN and DH separately.
 }

 \references{
  [1] A.B. Olshen, H. Bengtsson, P. Neuvial, P. Spellman, R.A. Olshen, V.E. Seshan, \emph{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation}, Bioinformatics 2011 (in press)\cr
 \cr
   [2] H. Bengtsson, P. Neuvial and T.P. Speed, \emph{TumorBoost: Normalization of allele-specific tumor copy numbers from a single pair of tumor-normal genotyping microarrays}, BMC Bioinformatics, 2010.
 }


\keyword{IO}