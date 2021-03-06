\name{karyotype}
\alias{karyotype}


\title{
karyotype plotting function
}

\description{
A general overview of the whole genome of the sample.
}

\usage{
karyotype(chr,start,end,int,ai,mchr,mpos,mval,schr,spos,sval,name='',xlim=c(0,2.5),ylim=c(0,1),MAPD=" ",MHOF=" ")
}

\arguments{
  	\item{chr}{					%segs$chr
  	Chromosomes name of segs object. Final version of segs created from patchwork.Medians_n_AI().
	}
  	\item{start}{				%segs$start
	Start position of segment.
	}
	\item{end}{					%segs$end
  	End position of segment.	
	}
	\item{int}{					%segs$mean
  	Normalized coverage.
	}
	\item{ai}{					%segs$ai
  	Allelic imbalance.
	}
	\item{mchr}{
	Chromosomes name of kbsegs object. Final version of kbsegs created from patchwork.smoothing().
	}
	\item{mpos}{
  	Middle position of kbsegs segments.
	}
	\item{mval}{
  	Ratio between kbsegs coverage and refcoverage.
	}
	\item{schr}{
  	Chromosome names of alf object. Final version of alf created from patchwork.alleledata().
	}
	\item{spos}{
  	Positions in alf.
	}
	\item{sval}{
	Inverse ratio between allele maximum and allele minimum.
	}
	\item{name}{
  	Default is the name of the input BamFile from patchwork.plot().
  	Plot generated as "name_overview.jpg" in your working directory.
	}
  	\item{xlim}{
  	Default is c(0,2.5). Limit of x axis.
	}
  	\item{ylim}{
  	Default is 0:1. Limit of y axis.
	}
  	\item{MAPD}{
  	Placeholder for later sample measures.
	}
  	\item{MHOF}{
  	Placeholder for later sample measures.
	}
}

\details{
The plot is a overview, for a closer look see the plots generated by karyotype_chroms().

Walkthrough of the plot:

\bold{TOP} \cr
Vertical axis: Allelic Imbalance \cr
Horizontal axis: Normalized Coverage \cr

The individual chromsomes in color plotted against the complete genome in grey.

\bold{TOP LOWER} \cr
Vertical axis: Normalized Coverage \cr
Horizontal axis: Coordinate & Chromosome\cr

Normalized coverage over coordinate for whole genome, chromosomally. See bottom for chromosome axis.

\bold{MIDDLE LOWER} \cr

Cytoband information for whole genome, chromosomally.

\bold{BOTTOM LOWER} \cr
Vertical axis: Allelic Imbalance \cr
Horizontal axis: Coordinate & Chromosome \cr

Allelic imbalance over coordinate for whole genome, chromosomally.

}

%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Markus Mayrhofer, \email{markus.mayrhofer@medsci.uu.se}
Sebastian DiLorenzo, \email{sebastian.dilorenzo@medsci.uu.se}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{patchwork.plot}} \cr
\code{\link{karyotype_chroms}}
}

%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%Run examples of karyotype stuff. or not, maybe just point back to patchwork.plot.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line