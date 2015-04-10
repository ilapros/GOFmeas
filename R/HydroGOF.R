
# library(lmom)    ## needed to fit the L-moments, generate the Monte Carlo sample from the Kappa distribution and draw the L-moment plot
# library(ellipse) ## needed to draw the ellipse



HWt4<-function(x,Distribution)
{
  # Return H&W population value of t4 based on t3 value (x) and Distribution
  # Coefficients from Hosking & Wallis (1997),
  if(Distribution=="GLO") {pdist<-c(0.16667,0,0.83333,0,0,0,0,0,0)}
  if(Distribution=="GEV") {pdist<-c(0.10701,0.11090,0.84838,-0.06669,0.00567,-0.04208,0.03763,0,0)}
  if(Distribution=="GNO") {pdist<-c(0.12282,0,0.77518,0,0.12279,0,-0.13638,0,0.11368)}
  if(Distribution=="PE3") {pdist<-c(0.12240,0,0.30115,0,0.95812,0,-0.57488,0,0.19383)}
  
  px<-x^(0:8)
  yLmom<-t(pdist)%*%px
  HWt4<-yLmom
  HWt4
}


Reg.sim<-function(rmom, N.years, Nsim=500)
{
  #
  # Monte Carlo simulation of independent AMAX events from a homogeneous region
  # and output regional average L-moment ratios.
  #
  # The regional averages are formed as record-weighted averages (Hosking & Wallis, 1997)
  #
  #
  # Input:
  # ======
  # rmom:       vector containing the first four regional L-moment ratios (l1,l2,t3,t4) from the observed region
  # N.years:    vector containing the record-length at each of the i=1,..,M.station sites
  # N.sim:      number of replica regions to be generated; default = 500
  #
  #
  # Output:
  # =======
  # mcmom:    (Nsim*4) matrix of L-moments containing regional average L-moments for each of the Nsim generated regions
  #
  #
  #==================================================================================
  
  # M.stations: number of sites in the homogeneous region
  M.stations<-length(N.years)
  
  # Fit kappa of GLO distribution (if t4R is too large, see H&W 1997, section 4.3.5
  if((rmom[4]>=((5*rmom[3]^2+1)/6) || rmom[4]<=((5*rmom[3]^2-1)/4))) {
    rpar<-pelglo(rmom)
  } else {
    rpar<-pelkap(rmom)}
  
  # Simulate Nsim kappa distributed regions
  mcmom<-matrix(nrow=Nsim, ncol=4)
  for(i in 1:Nsim)
  {
    t.region<-matrix(nrow=M.stations,ncol=4)
    u<-runif(sum(N.years))
    if(rmom[4]>=((5*rmom[3]^2+1)/6) | rmom[4]<=((5*rmom[3]^2-1)/4))
      amax<-quaglo(u,rpar)
    else
      amax<-quakap(u,rpar)
  
    for(m in 1:M.stations)
    {
      flag<-1
      if(m==1) flag<-0
      t.region[m,]<-samlmu(amax[(1+flag*sum(N.years[1:(m-1)])):sum(N.years[1:m])])
    }
    mcmom[i,]<-apply(t.region,2,mean) # regional average L-SK and L-kur for m'th region
  }
  
  mcmom
}



zdist<-function(mcmom,rmom,Nsim)
{
  #
  # Hosking & Wallis (1993) GOF measure. Requires output from Reg.sim
  #
  # Input:
  # ======
  # mcmom:  (Nsim*4) Matrix containing regional average L-moments (l1,l2,t3,t4) for the Nsim MC generated homogeneous regions
  # rmom:   regional average L-moment ratios
  # Nsim:   number of simulations used in the creation of mcmom
  #
  # Output:
  # =======
  # zdist:    HW goodness-of-fit measure Z evaluated using each of the four distributions
  #           GLO, GEV, GNO, PE3
  #
  
  B4 <- sum((mcmom[,4]-rmom[4]))/Nsim
  s4 <- sqrt((sum((mcmom[,4]-rmom[4])^2)-Nsim*B4^2)/(Nsim-1))
  
  t4dist<-0
  t4dist[1]<-lmrglo(para=pelglo(rmom),nmom=4)[4] #GLO
  t4dist[2]<-lmrgev(para=pelgev(rmom),nmom=4)[4] #GEV
  t4dist[3]<-lmrgno(para=pelgno(rmom),nmom=4)[4] #LN3 
  t4dist[4]<-lmrpe3(para=pelpe3(rmom),nmom=4)[4] #PE3
  
  zdist<-(t4dist-rmom[4]+B4)/s4
  zdist
}



kp_gof<-function(mcmom, rmom, Nsim, plot, conf.lev, dist.type = "mahalanobis", ...)
{
  #
  # Kjeldsen & Prosdocimi GOF measure. Requires output from Reg.sim
  #
  # Input:
  # ======
  # mcmom:  (Nsim*4) Matrix containing regional average L-moments (l1,l2,t3,t4) for the Nsim MC generated homogeneous regions
  # rmom:   regional average L-moment ratios
  # Nsim:   number of simulations used in the creation of mcmom
  # plot:   logical, should a plot of the measure be displayed
  # conf.level: confidence level (1-alpha). If for a given distribution the distance between the bias-corrected tR vector and the theoretical L-skewness and L-kurtosis line is bigger than the qchisq(1-alpha, 2) value the distribution is not considered as an acceptable distribution. 
  # dist.type: the type of distance to be used to determine the distance between the bias-corrected regional parameters and the theoretical L-skewness and L-kurtosis line, default set to "mahalanobis", the other option ("geometric") computes the simple geometric distance
  #
  # Output:
  # ======= 
  # afsDist: KP goodness-of-fit measure T evaluated using each of the four distributions
  #            GLO, GEV, GNO, PE3. If a distribution is not accepted at the conf.lev value, NA is returned.
  #  
  # Link:
  # =====
  # HWt4:   return theoretical value of t4 for a given value of t3 and a specified distribution (GLO, GEV, GNO, PE3)
  #
  #
  #
  #===============================================================================
  
  # First, calculate bias, variance and covariance of regional averages of (t3,t4)
  B3<-sum((mcmom[,3]-rmom[3]))/Nsim
  B4<-sum((mcmom[,4]-rmom[4]))/Nsim
  
  v3 <-  (sum((mcmom[,3]-rmom[3])^2)-Nsim*B3^2)/(Nsim-1)
  v4 <-  (sum((mcmom[,4]-rmom[4])^2)-Nsim*B4^2)/(Nsim-1)
  v34 <- (sum((mcmom[,3]-rmom[3])*(mcmom[,4]-rmom[4]))-Nsim*B3*B4)/(Nsim-1)
  
  StR<-matrix(c(v3,v34,v34,v4),nrow=2,ncol=2)
  center=c(rmom[3],rmom[4])
  
  ## computes the mahalabnobis distances from a point on the distribution line to the bias-corrected estimates
  mahalanobisDistr <- function(xval, Dist, ...) mahalanobis(x = c(xval,HWt4(xval,Dist)), ...)
  outM <- rep(-9999, 4)
  names(outM) <- c("GLO","GEV","GNO","PE3")
  # find minimum mahalabnobis distance for each distribution
  for(j in c("GLO","GEV","GNO","PE3")) 
    outM[j] <- nlm(mahalanobisDistr,0.1,cov=StR, Dist = j, center = center - c(B3,B4))$minimum   
  ## is the minimum distance smaller than the conf.lev-th quantile of the chisq?
  outM[outM > qchisq(conf.lev, 2)] <- NA  

  if(plot==TRUE){
  	# little sneaky function that computes the t4 and checks if the mahalanobis difference is smaller than qchisq(conf.level,2) 
		# needed to determine which parts of the lines should be bold
		HWt4Sel <- function(x,Dist, center, StR, conf.lev) 
		  ifelse(mahalanobisDistr(x,Dist=Dist, center, StR) < qchisq(conf.lev,2), HWt4(x,Dist), NA)
	    ### a more elegant way would be to find the yLmom for which the mahalanobis distance exceeds qchisq - 
		  
        # Plot ellipse as well as bias-corrected regional values of (t3,t4).
        points(ellipse(StR,centre=(center-c(B3,B4)),level=conf.lev),type="l", lwd=2, ...)
        points(center[1]-B3,center[2]-B4,pch="+",cex=2, ...)
        x.int <- seq(-0.5,1, by = 0.002) ## very fine x grid 
        # Plot thick lines for part of distribution located within ellipse
        cc <- c(4,3,1,5); names(cc) <- c("GLO","GEV","GNO","PE3")
        for(j in c("GLO","GEV","GNO","PE3"))
        lines(x.int, 
            apply(as.matrix(x.int),1, HWt4Sel, Dist = j, center = center-c(B3,B4), StR=StR, conf.lev = conf.lev) , lwd = 2, col =cc[j])
     }
  afsDist <- c(-5,-5,-5,-5)
  if(dist.type == "geometric") {
    geomDistr <- function(x, Dist, center) dist(rbind(center, c(x,HWt4(x,Dist))) )  
    for(j in c("GLO","GEV","GNO","PE3")) 
       afsDist[j] <- nlm(geomDistr,0.1, Dist = j, center = center - c(B3,B4))$minimum
  }
  if(dist.type == "mahalanobis") afsDist <- outM
  afsDist
  
} # end of function

#===============================================================================



#' @name GOFmeasures
#' @title Goodness of fit measures
#' @description The main function to compute the Hosking and Wallis (1997) and Kjeldsen and Prosdocimi (2014) goodness of fit measures.  
#' @details The function can take information on the annual maxima information either as a list of stations in \code{stations} or as the summarised lmom matrix (\code{lmom}) and information of the length of each series (\code{n.amax}).
#' @references Hosking, J. R. M., and J. R. Wallis (1997), Regional frequency analysis: an approach based on L-moments, Cambridge University Press.  
#' @references Kjeldsen T. R. and Prosdocimi I. (2015) A bivariate extension of the Hosking and Wallis goodness-of-fit measure for regional distributions. Water Resources Research. DOI: 10.1002/2014WR015912
#' @author Thomas Kjeldsen and Ilaria Prosdocimi
#' 
#' @param stations a list of amax, with one list components per station (see \code{Details})
#' @param lmom the matrix of sample lmoments: each line corresponds to the lmoments of for each station (see \code{Details})
#' @param n.amax a vector with the length of years recorded in each station (see \code{Details})
#' @param mcmom  an (Nsim*4) Matrix containing regional average L-moments (l1,l2,t3,t4) for the Nsim MC generated homogeneous regions. By default this is left \code{NULL} and the mcmom matrix is computed. 
#' @param Nsim the number of simulations used in the creation of mcmom, if \code{mcmom} is left null.
#' @param type indicates which GOF measure should be calculated, default is both "HW" and "KP"
#' @param plot logical vector indicating if the lmoment diagram with the ellipses intersecting the theoretical lines should be displayed
#' @param conf.lev confidence level (1-alpha). If for a given distribution the distance between the bias-corrected tR vector and the theoretical L-skewness and L-kurtosis line is bigger than the qchisq(1-alpha, 2) value the distribution is not considered as an acceptable distribution. 
#' @param dist.type  the type of distance to be used in the Kjeldsen and Prosdocimi measure to determine the distance between the bias-corrected regional parameters and the theoretical L-skewness and L-kurtosis line, default set to "mahalanobis", the other option ("geometric") computes the simple geometric distance
#' 
#' @return An object of the GOFmeas class. 
#' \item{GOFtable}{a table with the Goodness of fit results}
#' \item{selectedDist}{the selected distribution which minimises the goodness of fit measures calculated}
#' \item{mcmom}{the (Nsim*4) mcmom matrix of sample lmoments from the MCgenerated homogeneous region. The matrix can be reused if needed}
#' \item{lmom}{the matrix of sample lmoments for the region under study}
#' \item{rmom}{the weighted average regional lmoment. Weights are given according to the length of sample at each station}
#' \item{n.amax}{a vector with the number of years recorded in each station}
#' @seealso \code{\link{print.GOFmeas}} \code{\link{summary.GOFmeas}} 
#' @examples set.seed(546)
#' # generate a region with (mostly) GEV distribution
#' library(lmom)
#' stats <- list(quagev(runif(30,0,1),c(289, 45 ,-0.22)),
#'               quagev(runif(45,0,1),c(189, 40 ,-0.15)),
#'               quagev(runif(25,0,1),c(122, 10 ,-0.24)),
#'               quagev(runif(43,0,1),c( 59,  8 ,-0.18)),
#'               quagev(runif(32,0,1),c( 62, 10 ,-0.21)),
#'               quagev(runif(28,0,1),c( 91,  9 ,-0.25)),
#'               quaglo(runif(27,0,1),c(202, 25 ,-0.17)))
#' tt <- GOFmeasures(stats)
#' tt
#' summary(tt)
#' @export 
#' @importFrom ellipse ellipse
#' @import lmom
GOFmeasures <- function(stations=NULL,lmom=NULL,n.amax=NULL,Nsim=500,mcmom=NULL,type=c("HW_GOF","KP_GOF"),
                   plot=FALSE,conf.lev=0.9, dist.type = "mahalanobis", ...){
  # a unified function which computes the HW distance measure and the Kjeldsen & Prosdocimi GOF measure.
  #
  # Input:
  # ======
  # amax information given either as a list of stations or as the summarised lmom matrix and information of the length of each series in n.amax
  # stations: a list of amax, with one list components per station -
  # alternatively to the list of stations one can input
  # lmom:  the matrix of sample lmoments for each station
  # n.amax: a vector with the length of years recorded in each station
  # for the GOF measures to be calculated the following is needed
  # mcmom:  (Nsim*4) Matrix containing regional average L-moments (l1,l2,t3,t4) for the Nsim MC generated homogeneous regions
  #         if mcmom=NULL computed via Reg.sim.
  #         Otherwise it should be a matrix of lmoments from a simulated region with characteristics similar to the one under study
  # Nsim:   number of simulations used in the creation of mcmom
  # type: indicates which GOF measure should be calculated, default is both "HW" and "KP"
  # plot: indicates if the lmoment diagram with the ellipses intersecting the theoretical lines should be plotted
  # conf.level: confidence level (1-alpha). If for a given distribution the distance between the bias-corrected tR vector and the theoretical L-skewness and L-kurtosis line is bigger than the qchisq(1-alpha, 2) value the distribution is not considered as an acceptable distribution. 
  # dist.type: the type of distance to be used in the Kjeldsen and Prosdocimi measure to determine the distance between the bias-corrected regional parameters and the theoretical L-skewness and L-kurtosis line, default set to "mahalanobis", the other option ("geometric") computes the simple geometric distance
  
  # Output:
  # =======
  # GOFtable: the table with the Goodness of fit results
  # selectedDist: the selected distribution which minimises the goodness of fit measure
  # mcmom: the mcmom matrix of sample lmoments from the MC generated homogeneous region. The matrix can be reused if needed
  # lmom: the matrix of sample lmoments for the region under study
  # rmom: the weighted average regional lmoment
  # n.amax: a vector with the length of years recorded in each station
  #
  if(!any(type %in% c("KP_GOF","HW_GOF"))) {warning("measure type not recognised, doing KP_GOF and HW_GOF"); type <- c("KP_GOF","HW_GOF")}
  if(is.list(stations) & !any(is.null(lmom) & is.null(n.amax))) warning("when a list of stations is given, lmoments and n.amax are directly calculated")
  if(!is.list(stations) & (is.null(lmom) | is.null(n.amax))) stop("please provide either lmom and n.amax or a list of station records")
  if(is.list(stations)){
    n.amax <- unsplit(lapply(stations,function(x) length(x[!is.na(x)])),seq(1,length(stations)))
    if(any(n.amax < 4)) warning("Some stations have less than 4 years of data: these will be excluded from the analysis")
    lmom <- do.call(rbind,lapply(stations,samlmu))
    lmom <- lmom[n.amax>3,]; n.amax <- n.amax[n.amax>3]
  }
  if(nrow(lmom)!=length(n.amax)) stop("the number of years does not match with the lmom matrix dimension")
  M.stations <- length(n.amax) ## total number of stations
  
  # Calculate record-length weighted regional average L-moment ratios (rmom)
  sum.n<-sum(n.amax)
  w<-n.amax/sum.n  # weights for each site
  rmom<-apply(lmom,2,weighted.mean,w) ## regional moments
  
  if(is.null(mcmom)) mcmom <- Reg.sim(rmom=rmom,N.years=n.amax,Nsim=Nsim)
  
  out<-matrix(-5,ncol=4,nrow=0)
  colnames(out) <- c("GLO","GEV","GNO","PE3")
  # Hosking & Wallis measure
  if("HW_GOF" %in% type) out<-rbind(out,zdist(mcmom=mcmom,rmom=rmom,Nsim=Nsim))
  
  if(plot){
    lmrd(distributions = "GLO GEV GNO PE3", 
         xlim = range(c(0,0.6,quantile(mcmom[,3],c(0.005,0.995)))), 
         ylim = range(c(0,0.4,quantile(mcmom[,4],c(0.005,0.995)))))
    points(lmom[,3],lmom[,4],pch="+",cex=1.0,col="darkgrey") 
  }
  # KP measure
  if("KP_GOF" %in% type) out<-rbind(out,kp_gof(mcmom=mcmom,rmom=rmom,Nsim=Nsim,plot=plot,conf.lev=conf.lev,dist.type=dist.type, ...))
  rownames(out) <- type
  
  ## a little trick in case the KP output consists of all NA
  out <- cbind(out,rep(5500,NROW(out)))
  colnames(out) <- c("GLO","GEV","GNO","PE3","none")
  zz<-t(apply(abs(out), 1, function(x) {c(min(x,na.rm=TRUE),names(x)[which.min(x)])}))
  zz<-as.data.frame(zz)
  names(zz)<-c("minvalues","Distribution")
  zz$minvalues<-as.numeric(as.character(zz$minvalues))
  zz$Distribution<-as.character(zz$Distribution)
  zz[,1][round(zz[,1]) == 5500] <- -999
  out <- out[,1:4]
  res<-list(GOFtable=out,selectedDist=zz,mcmom=mcmom,lmom=lmom,rmom=rmom,n.amax=n.amax)
  class(res) <- "GOFmeas"
  res
}


# methods for class GOFmeas

#' @name GOFmeas
#' @export
print.GOFmeas <- function(object) print(round(object$GOFtable,4))

#' @name GOFmeas
#' @title The GOFmeas class 
#' @description Measures of goodness of fit 
#' @param object an object of the GOF class, typically the output of a GOFmeasures call
#' @seealso \code{\link{GOFmeasures}}
#' @export
summary.GOFmeas <- function(object) data.frame(minvalues = round(object$selectedDist[,1],4), Distribution = object$selectedDist[,2], row.names = rownames(object$selectedDist))



