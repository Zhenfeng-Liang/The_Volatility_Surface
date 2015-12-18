library(stinepack)

# Function to compute total variance given k and t
# texp is a vector of times to expiration
sviW <- function(sviMatrix,texp,k,t){

# Vector of SVI variance for a given strike k
sviWk <- function(k){
m <- dim(sviMatrix)[1];
sapply(1:m,function(i){svi(sviMatrix[i,],k)}); 
}

# Function to compute interpolated variance for any strike and expiration
wInterp <- function(k,t){stinterp(texp,sviWk(k),t)$y};

# Vectorized function that returns implied total variance for vectors k and t
return(sapply(k,function(k1){wInterp(k1,t)}));

} # End of sviW

#######################################################################################################
# SVI local variance function
sviLocalVar <- function(sviMatrix,texp,k,t){

sviMatrix_columns <- dim(sviMatrix)[2];

# Check to see if number of columns is correct
if (sviMatrix_columns != 5) {
print("Error: Wrong number of columns to be SVI data");
return(NA);
}

deltak <- 0.01;
deltat <- 0.01*(t>0.01)+t/2*(t<=0.01); # Make sure we don't end up with negative times!


# Estimate derivatives
wkt <- sviW(sviMatrix,texp,k,t); 
wktm <- sviW(sviMatrix,texp,k,t-deltat); 
wktp <- sviW(sviMatrix,texp,k,t+deltat);
dwdt <- (wktp-wktm)/(2*deltat);# Central difference!
wkmt <- sviW(sviMatrix,texp,k-deltak,t);
wkpt <- sviW(sviMatrix,texp,k+deltak,t);
dwdk <- (wkpt-wkmt)/(2*deltak);
d2wdk2 <- (wkpt+wkmt-2*wkt)/(deltak^2);

tmp <- dwdt/(1 - k/wkt*dwdk + 1/4*(-1/4 - 1/ wkt + k^2/wkt^2)*(dwdk)^2 + 1/2*d2wdk2);

return(tmp);

}
