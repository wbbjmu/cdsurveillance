## function to conduct space-time permutation using Poisson scan statistics

library(reliaR)
library(stringr)
library(scanstatistics)
library("Rglpk")
library(reshape2)

## function
sp_perm <- function(postcode=NULL, cases=NULL, date=NULL, lat=NULL, long=NULL, max_timeinc=NULL, max_distance=0.5, n=NULL, 
                    bootstraps=5,p_method=NULL,significance_level=NULL, data){
  set.seed(123)
  ## function to define gumbel_pvalue
  gumbel_pvalue <- function(observed, replicates, method = "ML") {
    # Fit Gumbel distribution to Monte Carlo replicates
    gumbel_mu <- NA
    gumbel_sigma <- NA
    if (method == "ML") {
      gum_fit <- gum.fit(replicates, show = FALSE)
      gumbel_mu <- gum_fit$mle[1]
      gumbel_sigma <- gum_fit$mle[2]
    } else {
      gumbel_sigma <- sqrt(6 * var(replicates) / pi^2)
      gumbel_mu <- mean(replicates) + digamma(1) * gumbel_sigma
    }
    
    pvalue <- pgumbel(observed, gumbel_mu, gumbel_sigma, lower.tail = FALSE)
    
    return(list(pvalue = pvalue, 
                gumbel_mu = gumbel_mu, 
                gumbel_sigma = gumbel_sigma))
  }
  
  # For greater precision, upped digits to 10.
  options(digits=10)
  combined_sample_result <- as.data.frame(matrix())
  combined_total <- matrix(,1,1)
  
  # Format dataset.
  
  # Limit the dataset to five variables in the following order: ID, case status, date (in YYYY-MM-DD format), latitude and longitude.
  # Based on what the user told R, ClustR extracts these variables and puts them in the order it wants.
  overall_data <- data[,c(postcode,cases,date,lat,long)]
  
  # Rename the columns for clearer use in later code.
  colnames(overall_data) <- c("postcode","cases","date","lat","long")
  
  # Checking if options are correctly specified
  if (is.null(postcode)) stop("You must specify an identification variable.")
  if (is.null(cases)) stop("You must specify an disease status variable.")
  if (is.null(date)) stop("You must specify a date variable.")
  if (is.null(lat)) stop("You must specify a latitude variable.")
  if (is.null(long)) stop("You must specify an longitude variable.")
  if (is.null(data)) stop("You must specify what dataframe to use.")
  
  
  # Tell R that the birthdate values are dates.
  overall_data$D <- as.Date(overall_data$date)
  
  # Determine the minimum date in the dataset and create a new column repeating the minumum date at every row.
  min_date <- min(overall_data$D)
  overall_data$min_date <- min_date
  
  # Calculate the number of days that have passed for each location between the earliest date in the dataset and each location's date and put this in a new column called *time*.
  overall_data$time <- as.numeric(difftime(overall_data$D,overall_data$min_date), units = "days")
  
  # Drop any person with missing data.
  overall_data <- overall_data[complete.cases(overall_data),]
  
  # Evaluate User Options
  
  # Get information about what the starting and ending dates for the data are.
  date_range = c(as.character(overall_data$date[which.min(overall_data$time)]), as.character(overall_data$date[which.max(overall_data$time)]))
  mintime <- min(overall_data$time)
  maxtime <- max(overall_data$time)
  time_range_vector <- mintime:maxtime
  
  # Drop the variables for birthdate and minimum date. (Now that we're done with them)
  data_true <- overall_data[,c("postcode","cases","time","lat","long")]
  
  # Limit the time ranges the user can look at.
  if (is.null(max_timeinc)){
    max_timeinc=ceiling(maxtime/2)
  } else {
    try(if(max_timeinc>maxtime) stop(paste("The time inc must be less than the time range in your data:",maxtime,"days")))
    max_timeinc <- ceiling(max_timeinc/2)
  }
  
  # Limit the distance the user can look at.
  try(if(max_distance<0) stop("The max_distance must be proportion from 0 to 1"))
  try(if(max_distance>1) stop("The max_distance must be proportion from 0 to 1"))
  
  ## calculate distances for each pair of locations
  dist_x <- data_true$long*pi/180
  dist_y <- data_true$lat*pi/180
  
  cos_x <- matrix(rep(dist_x,length(dist_x)),nrow = length(dist_x),byrow = F)-matrix(rep(dist_x,length(dist_x)),nrow = length(dist_x),byrow = T)
  cos_x <- cos(cos_x)
  
  sin_y <- sin(matrix(rep(dist_y,length(dist_y)),nrow = length(dist_y),byrow = F))*sin(matrix(rep(dist_y,length(dist_y)),nrow = length(dist_y),byrow = T))
  cos_y <- cos(matrix(rep(dist_y,length(dist_y)),nrow = length(dist_y),byrow = F))*cos(matrix(rep(dist_y,length(dist_y)),nrow = length(dist_y),byrow = T))
  
  dist_location <- 6371*acos(pmin(cos_y*cos_x+sin_y,1))  ## matrix of distances
  max_radius <- max(dist_location)*max_distance  ## max distances
  
  
  
  #The size of samples to bootstrap.
  
  # If the user did not specify the size of the samples to take: if sample size is small, boostrap as large of samples as possible; if sample size is large, boostrap samples of 1000.
  # If the user did specify the size of the samples to take: use that number.
  if (is.null(n)){
    if(nrow(data)< 1001) {
      samples <- nrow(data)
    } else {
      samples <- 1000
    }
  } else {
    samples <- n
  }
  
  # Permutated Data
  # Create a 1x1 matrix called *total* to hold the null distribution.
  
  total <- matrix(,1,1)
  
  ## Repeat the following process [bootstraps] times in order to get a stable distribution:
  glr_max <- 1
  
  # Create a new dataset called *data_perm* with case/control status scrambled randomly across participants (using the same proportion as found in the original dataset). This is your permuted dataset.
  # shuffle the time for each location
  # Create a new dataset called *data_perm* with case/control status scrambled randomly across participants (using the same proportion as found in the original dataset). This is your permuted dataset.
  data_perm <- data_true
  
  ## 生成所有地区与时间的可能组合
  cases_time <- data_perm$time[!duplicated(data_perm$time)]
  cases_location <- data_perm$postcode[!duplicated(data_perm$postcode)]
  
  data_perm <- data_perm[order(data_perm$postcode,data_perm$time),]
  
  
  simu_potential <- dcast(data = data_perm[,c("postcode","time","cases")],
                          postcode~time,value.var = cases)
  
  
  ############################################################
  
  #NA变为0
  simu_potential[is.na(simu_potential)] <- 0
  #丢弃postcode变量
  rownames(simu_potential) <- simu_potential$postcode
  dropvar <- names(simu_potential) %in% c("postcode")
  simu_potential <- simu_potential[!dropvar]
  
  row_sum <- apply(simu_potential, 1, sum)
  col_sum <- apply(simu_potential, 2, sum)
  
  
  
  A <- matrix(c(rep(c(rep(1,ncol(simu_potential)),rep(0,nrow(simu_potential)*ncol(simu_potential))),nrow(simu_potential)-1),rep(1,ncol(simu_potential))),
              nrow = nrow(simu_potential),ncol = ncol(simu_potential)*nrow(simu_potential), byrow = T)
  
  
  B <- matrix(c(rep(c(rep(c(1,rep(0,ncol(simu_potential)-1)),nrow(simu_potential)),0),ncol(simu_potential)-1),c(rep(c(1,rep(0,ncol(simu_potential)-1)),nrow(simu_potential)-1),1)),
              nrow = ncol(simu_potential),ncol = ncol(simu_potential)*nrow(simu_potential), byrow = T)
  
  C <- rbind(A,B)
  
  b <- c(row_sum,col_sum)
  
  solve_result <- Rglpk_solve_LP(obj=rep(0,ncol(simu_potential)*nrow(simu_potential)),C,dir = rep("==",nrow(C)),rhs = b,types = "I")
  
  raw_new <- matrix(solve_result$solution,nrow = nrow(simu_potential),ncol = ncol(simu_potential),byrow = T)
  
  raw_new <- as.data.frame(raw_new)
  names(raw_new) <- names(simu_potential)
  raw_new$postcode <- row.names(simu_potential)
  
  perm_newdata <- melt(raw_new,ID="postcode")
  names(perm_newdata)[c(2,3)] <- c("time","cases")
  perm_newdata <- perm_newdata[perm_newdata$cases!=0,]
  
  xy_info <- data_perm[!duplicated(data_perm$postcode),c("postcode","lat","long")]
  perm_newdata <- merge(perm_newdata,xy_info,by="postcode")
  
  perm_newdata$time <- as.numeric(perm_newdata$time)
  
  # Rename *data_perm* as *data*.
  data <- perm_newdata
  # Number each row of the permuted dataset starting at 1 and ending at [number of rows in the dataset] in a new column called *num_id*.
  r <- nrow(data)
  data$num_id <- seq(1,r)
  
  ## calculate total cases in dataset
  total_cases <- sum(data$cases)
  ## calculate the expected cases in each location with its time
  for (cxk in 1:nrow(data)) {
    data$expec_cases[cxk] <- (sum(data$cases[data$postcode==data$postcode[cxk]])*sum(data$cases[data$time==data$time[cxk]]))/total_cases
  }
  
  ## define function to calculate 
  dist_cal <- function(long1,long2,lat1,lat2){
    x1 <- long1*pi/180
    x2 <- long2*pi/180
    y1 <- lat1*pi/180
    y2 <- lat2*pi/180
    
    dist_lyf <- 6371*acos(pmin(cos(y1)*cos(y2)*cos(x1-x2)+sin(y1)*sin(y2),1))
    return(dist_lyf)
  }
  
  
  for (xz in 1:samples) {
    zyx <- sample(r,1)
    
    lat <- data[zyx,"lat"]
    long <- data[zyx,"long"]
    time <- data[zyx,"time"]
    
    cylinder_data <- data[(dist_cal(long1 = long,long2 = data$long,lat1 = lat,lat2 = data$lat)<max_radius)&(abs(data$time-time)<max_timeinc),]
    cylinder_data$distance <- dist_cal(long1 = long,long2 = cylinder_data$long,lat1 = lat,lat2 = cylinder_data$lat)
    cylinder_data$time_diff <- abs(cylinder_data$time-time)
    
    cylinder_data_order <- cylinder_data[order(cylinder_data$time_diff,cylinder_data$distance),]
    rm(cylinder_data)
    
    glr <- 0
    for (lh in 1:nrow(cylinder_data_order)) {
      glr_data <- cylinder_data_order[1:lh,]
      ob_cases <- sum(glr_data$cases)
      exp_cases <- sum(glr_data$expec_cases)
      glr[lh] <- ob_cases*log(ob_cases/exp_cases)+(total_cases-ob_cases)*log((total_cases-ob_cases)/(total_cases-exp_cases))
    }
    
    glr_max[xz] <- max(glr)
    
  }
  
  rm(data_perm, data, cylinder_data_order,glr_data)
  glr_max <- glr_max[order(glr_max)]
  log(glr_max)
  
  # (End of For Loop for repeating this process [bootstaps] times.)
  # Calculate and save the cut offs for the 97.5th, 99th, and 100th percentile of your chi-square values based on your randomly sampled points from your null distribution.
  
  # Quantiles are based on the distribution of glr over every bootstrap performed.
  quantile(glr_max, c(.975,.99, 1), na.rm = TRUE)
  nsfpct <- quantile(glr_max, .975, na.rm = TRUE)
  nnpct <- quantile(glr_max, .99, na.rm  = TRUE)
  max <- quantile(glr_max, 1, na.rm = TRUE)
  ## create the cumulative distribution  function for glr
  percentile <- ecdf(glr_max)
  
  # Real Data
  data <- data_true
  
  r <- nrow(data)
  data$num_id <- seq(1,r)
  
  ## calculate total cases in dataset
  total_cases <- sum(data$cases)
  
  ## calculate the expected cases in each location with its time
  for (cxk in 1:nrow(data)) {
    data$expec_cases[cxk] <- (sum(data$cases[data$postcode==data$postcode[cxk]])*sum(data$cases[data$time==data$time[cxk]]))/total_cases
  }
  
  data$p_value <- NA
  data$p_gumbel <- NA
  data$p_mc <- NA
  for (xz in 1:r) {
    lat <- data[xz,"lat"]
    long <- data[xz,"long"]
    time <- data[xz,"time"]
    
    cylinder_data <- data[(dist_cal(long1 = long,long2 = data$long,lat1 = lat,lat2 = data$lat)<max_radius)&(abs(data$time-time)<max_timeinc),]
    cylinder_data$distance <- dist_cal(long1 = long,long2 = cylinder_data$long,lat1 = lat,lat2 = cylinder_data$lat)
    cylinder_data$time_diff <- abs(cylinder_data$time-time)
    
    cylinder_data_order <- cylinder_data[order(cylinder_data$time_diff,cylinder_data$distance),]
    rm(cylinder_data)
    
    glr <- 0
    for (lh in 1:nrow(cylinder_data_order)) {
      glr_data <- cylinder_data_order[1:lh,]
      ob_cases <- sum(glr_data$cases)
      exp_cases <- sum(glr_data$expec_cases)
      glr[lh] <- ob_cases*log((ob_cases/exp_cases))+(total_cases-ob_cases)*log((total_cases-ob_cases)/(total_cases-exp_cases))
    }
    
    data$glr_max[xz] <- max(glr)
    
    glr_mc <- c(glr_max,max(glr))
    glr_mc <- glr_mc[order(glr_mc)]
    
    data$p_mc[xz] <- sprintf("%.10f",1-(which(glr_mc==max(glr))[1]/length(glr_mc)))
    
    data$p_value[xz] <- sprintf("%.10f",1-percentile(max(glr)))
    data$glr_index[xz] <- which.max(glr)
    
    ## when glr is infinite, gumbel_pvalue cannot be calculated, so reset the gumbel_pvalue to 0 in this occasion
    
    tryCatch(
      (if (is.infinite(max(glr))){
        data$p_gumbel[xz] <- 0.000000000000000000
      } else {
        tryCatch(
          (data$p_gumbel[xz] <- as.numeric(gumbel_pvalue(max(glr),replicates = glr_max,method="ML")[1])), error = function(e){
            data$p_gumbel[xz] <- as.numeric(gumbel_pvalue(max(glr),replicates = glr_max,method="MoM")[1]) 
          }
        )
      }), error= function(e){
        data$p_gumbel[xz] <- NA
      }
      
    ) 
  }
  data$p_fdr <- NA
  data$p_fdr <- sprintf("%.10f",p.adjust(data$p_value, method = "fdr", n = nrow(data)))
  
  ## check significance level
  if (is.null(significance_level)){
    significance_level <- 0.05
  } 
  
  try(if(significance_level>1) stop("The significance level should be less than 1"))
  try(if(significance_level<0) stop("The significance level should be more than 0"))
  
  ## check methods user specified to calculate the p value
  if (is.null(p_method)){
    p_method= "MC"
    index=4
  } else if(p_method== "crude"){
    index=1
  } else if(p_method== "FDR"){
    index=2
  } else if(p_method== "Gumbel"){
    index=3
  } else if(p_method== "MC"){
    index=4
  }
  
  try(if(!(p_method %in% c("crude", "FDR", "Gumbel","MC"))) stop("The method to calculate p_value should be one of them:'crude','FDR','Gumbel','MC'."))
  p_names <- c("p_value","p_fdr","p_gumbel","p_mc")
  
  ## calculate number of P < 0.05
  
  p_qualify <- length(which(data[,p_names[index]]<significance_level))
  if (p_qualify==0){
    return(paste("No cluster was found in your data at the significance_level:",significance_level,sep = ""))
  } else {
    
    ## find cases with P < 0.05
    report_cases <- data[data[,p_names[index]]<significance_level,]
    report_cases$cluster <- seq(1:nrow(report_cases))
    
    ## report the total number of cluster
    cluster_number <- nrow(report_cases)
    hzt <- 1
    cluster_detail <- list()
    for (hzt in 1:nrow(report_cases)) {
      index <- report_cases$num_id[hzt]
      glr_index <- report_cases$glr_index[hzt]
      
      lat <- data[data$num_id==index,"lat"]
      long <- data[data$num_id==index,"long"]
      time <- data[data$num_id==index,"time"]
      
      cylinder_data <- data[(dist_cal(long1 = long,long2 = data$long,lat1 = lat,lat2 = data$lat)<max_radius)&(abs(data$time-time)<max_timeinc),]
      cylinder_data$distance <- dist_cal(long1 = long,long2 = cylinder_data$long,lat1 = lat,lat2 = cylinder_data$lat)
      cylinder_data$time_diff <- abs(cylinder_data$time-time)
      
      cylinder_data_order <- cylinder_data[order(cylinder_data$time_diff,cylinder_data$distance),]
      rm(cylinder_data)
      
      cylinder_data_order <- cylinder_data_order[1:glr_index,]
      cylinder_data_order$date <- cylinder_data_order$time+min_date
      cylinder_data_order <- cylinder_data_order[,c("postcode","cases","date","lat","long","expec_cases","distance","time_diff")]
      cluster_detail[[hzt]] <- cylinder_data_order
      location_order <- cylinder_data_order$postcode[!duplicated(cylinder_data_order$postcode)]
      location_order <- location_order[order(location_order)]
      cluster_location_range <- str_c(location_order,collapse = ",")
      cluster_date_range <- str_c(range(cylinder_data_order$date),collapse = " to ")
      
      report_cases$locationIDs_included[hzt] <- cluster_location_range
      report_cases$time_frame[hzt] <- cluster_date_range
    }
    
    
    ## Return the results
    cluster_info <- list(report_cases=report_cases,cluster_detail=cluster_detail)
    return(cluster_info)
  }
}




