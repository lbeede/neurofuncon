\name{neurofun3d}
\alias{neurofun3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
neurofun3d
}
\description{
 This function estimates and plots a Gaussian Graphical Model (GGM via graphical lasso) in 3D describing the functional connectivity inferred from electro-physiological data (for example, calcium imaging data) collected from neurons on multiple slices of a 3D space.
}
\usage{
neurofun3d(traces, position, all.stimuli.ids, nlayers, ggm.stimuli.ids, num.neurons = 25, num.edges = 25, rho = seq(5,0.001,-0.001), Plot = TRUE, node.size = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traces}{
    n x d matrix where n is the number of samples and d is the total number of neurons.
}
  \item{position}{
     A d x 3 matrix containing the 3D coordinates (x,y,z) of the d neurons.
}
  \item{all.stimuli.ids}{
     A vector containing the stimuli ids (numeric) related to the n samples.
}
  \item{nlayers}{
    The layers (z-axis) desired to create the GGM.
}
  \item{ggm.stimuli.ids}{
     The stimuli ids desired to create the GGM.
}
  \item{num.neurons}{
     Number of neurons with the highest activation means desired to create the GGM.
}
  \item{num.edges}{
     Number of edges desired to create the GGM.
}
  \item{rho.seq}{
      Sequence of penalty parameters of the graphical lasso. % Don't forget ebic implementation
}
  \item{Plot}{
      Boolean value to plot the graph.
  }
    \item{node.size}{
      Maximum size of the highest activation means nodes
  }
}
\details{
The graph produced by neurofun3d will have nodes with sizes proportional to neurons’ average activities and colors depending on their variances. The edges indicate positive (blue) or negative (red) covariance between neurons in response to the analyzed stimulus.
}
\value{
%%  ~Describe the value returned
\item{nlayers}{The layer id of the GGM generated}

\item{location.positive}{The starting and end points (x,y,z) of the edges with positive neuronal functional correlation}

\item{location.negative}{The starting and end points (x,y,z) of the edges with negative neuronal functional correlation}

\item{layers}{All information about the layer's neurons. V1, V2, V3 are the x, y, z location of the neurons. 'averages' is the calculated averages of the neurons' activations. 'variances' is the calculated variances of the neurons' activations. 'size' is the plotted size of the nodes. If the corresponding neuron was NOT one of the top num.neurons highest activation neurons, size will be 0.5. Else, size will be equivalent to the average. 'alpha' is the plotted opacity of the nodes. If the corresponding neuron was NOT one of the top num.neurons highest activation neurons, alpha will be 0.4. Else, alpha will be 1.}

\item{rho}{The calculated rho value calculated using graphical lasso}

\item{true.num.edges}{Due to the nature of the code, sometimes the user inputed num.edges will not be equivalent to the true number of edges plotted.}

\item{true.num.neurons}{Due to the nature of the code, sometimes the user inputed num.neurons will not be equivalent to the true number of neurons plotted.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lauren Miako Beede <lbeede2@nd.edu> and Giuseppe Vinci <gvinci@nd.edu>

Department of Applied and Computational Mathematics and Statistics, University of Notre Dame
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(4) # to be reproducible

# Creating artificial data for this example
d <- 100 # number of neurons
n <- 1000 # number of time bins
Sigma = matrix(rnorm(d^2),nc=d)
Sigma = .1*cov2cor(Sigma\%*\%t(Sigma))

TRACES <- exp(MASS::mvrnorm(n,mu=rep(0,d),Sigma=Sigma))
POSITIONS <- cbind(matrix(round(runif(d*2)*1000),nc=2), sample(1:10,size=d,replace=TRUE)*100)
STIMULI <- sample(1:5,size=n,replace=TRUE)

# Example plot with 15 neuron-nodes on layer 100, 300, 500 and 20 edges, about stimuli 2 and 4
neurofun3d(traces=TRACES, position=POSITIONS, all.stimuli.ids=STIMULI, nlayer=c(100,300,500), ggm.stimuli.ids=c(2,4), num.neurons = 15, num.edges = 20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
