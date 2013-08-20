\name{cookie}
          \docType{data}
          \alias{cookie}
          \title{Near-Infrared (NIR) Spectroscopy of Biscuit Doughs}
          \description{
            This data set contains measurements from quantitative NIR spectroscopy. The example studied arises from an experiment done to test the feasibility of NIR spectroscopy to measure the composition of biscuit dough pieces (formed but unbaked biscuits). Two similar sample sets were 
made up, with the standard recipe varied to provide a large 
range for each of the four constituents under investigation: 
fat, sucrose, dry flour, and water. The calculated percentages 
of these four ingredients represent the 4 responses. There 
are 40 samples in the calibration or training set (with 
sample 23 being an outlier) and 
a further 32 samples in the separate prediction or validation set (with example 21 considered as an outlier). 

An NIR reflectance spectrum is available for each dough 
piece. The spectral data consist of 700 points measured 
from 1100 to 2498 nanometers (nm) in steps of 2 nm. 
          }
          \usage{data(cookie)}
          \format{A data frame of dimension 72 x 704. The first 700 columns correspond to the NIR reflectance spectrum, the last four columns correspond to the four constituents fat, sucrose, dry flour, and water. The first 40 rows correspond to the calibration data, the last 32 rows correspond to the prediction data.}
          \references{Please cite the following papers if you use this data set.

          P.J. Brown, T. Fearn, and M. Vannucci (2001) \emph{Bayesian Wavelet Regression on Curves with Applications to a Spectroscopic Calibration Problem}. Journal of the American Statistical Association, 96, pp. 398-408. 

     B.G. Osborne, T. Fearn, A.R. Miller, and S. Douglas (1984) \emph{Application of Near-Infrared Reflectance Spectroscopy to Compositional Analysis of Biscuits and Biscuit Dough}. Journal of the Science of Food and Agriculture, 35, pp. 99 - 105.
          }
         \examples{
    data(cookie) # load data
    X<-as.matrix(cookie[,1:700]) # extract NIR spectra
    Y<-as.matrix(cookie[,701:704]) # extract constituents
    Xtrain<-X[1:40,] # extract training data
    Ytrain<-Y[1:40,] # extract training data
    Xtest<-X[41:72,] # extract test data
    Ytest<-Y[41:72,] # extract test data
    }
          \keyword{datasets}
