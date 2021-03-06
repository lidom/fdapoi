% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_Functions.R
\docType{data}
\name{emotion_rating}
\alias{emotion_rating}
\title{Video Rating Data}
\format{A data frame with 65 rows and 168 variables}
\usage{
data(emotion_rating)
}
\description{
A dataset containing the emotion ratings reported from n=65 participants 
while watching an affective online video on the persecution of African albinos. 
A version of the video can be found online at YouTube (www.youtube.com). 
Link to the video: https://youtu.be/9F6UpuJIFaY 
The video clip used in the experiment corresponds to the first 112 sec. 
The first six data points (<1 sec.) are removed from the trajectories as they 
contain some obviously erratic components.
}
\details{
\itemize{
  \item Y: Final overall rating (Y=0 denotes 'I feel negative' and Y=1 denotes 'I feel not negative')
  \item X_1 to X_167: Continuous emotion ratings (discretization points 1 to 167) 
}
}
\keyword{datasets}
