#' Simulates a fishery using the Framework for Integrated Stock and Habitat Evaluation to guide their management decisions 
#' 
#' Simulate a fishery using an underlying Pella Tomlinson Surplus Production model, 
#'  add influences of climate change on biological growth rates, and make management decisions
#'  about how to change the fishing pressure based on target and limit  
#' @param b initial biomass of the fishery
#' @param r intrinsic growth rate of species 
#' @param r_s pace and magnitude of climate change as a percent decline in productivity per year (refer to references/climate_scenarios.md for examples)
#' @param r_p_s a proxy estimate for the pace and magnitude of climate change as a percent decline in productivity per year
#' @param error the amount of sampling error as a decimal (i.e. 0.1 for 10% error)
#' @param p the shape parameter of the surplus production model; default is 0.2 for the Pella Tomlinson model 
#' @param k the estimated carrying capacity relative to the initial biomass; default is 10000
#' @param years the length of the simulation; default is 100 years
#' @param hcr the harvest control rule used if the fishery falls between the target and limit expressed as the decimal 1-hcr (i.e. if you will reduce fishing pressure by 10% the hcr is 0.9)
#' @return a dataframe with 11 columns and rows corresponding the years parameter: 
#' \describe{
#' \item{b}{biomass in year (y)}
#' \item{c}{catch in year (y)}
#' \item{year}{year of the simulation}
#' \item{r}{growth rate in year (y)}
#' \item{r_p}{perceived growth rate in year (y) based on proxy estimate from the input r_p_s}
#' \item{f}{fishing mortality rate in year (y)}
#' \item{f_msy}{the fishing mortality rate that gives the biomass for maximum sustainable yield}
#' \item{f_msy_p}{the perceived fishing mortality rate that gives biomass for maximum sustainable yield, this is based on if and how managers are accounting for the influence of climate change on the underlying growth of the species}
#' \item{f_ratio}{the ratio of fishing mortality in year (y) to the fishing mortality that gives the biomass for maximum sustainable yield (i.e. f/f_msy)}
#' \item{f_ratio_p}{the perceived f/f_msy ratio based on the f_msy_p value instead of the true f_msy}
#' \item{f_ratio_err}{an f/f_msy ratio used to make a fishery management decision drawn from a lognormal distribution based on the mean of f_ratio_p and the sampling error defined by the error parameter}
#' }


sim_fishery <- function(b,
                        r,
                        r_s,
                        r_p_s,
                        error,
                        p = 0.2,
                        k = 10000,
                        years = 100,
                        hcr) {
  
  # Setup the results dataframe 
  results <- data.frame(
    b = rep(NA, years), c = rep(NA, years), 
    year = 1:years, r = rep(NA, years), r_p = rep(NA, years),
    f = rep(NA, years), f_msy = rep(NA, years),
    f_msy_p = rep(NA, years), f_ratio = rep(NA, years),
    f_ratio_p = rep(NA, years), f_ratio_err = rep(NA, years)
  ) 
  
  # Set the initial result for the outputs in year 1 using the input values
  results$b[1]   <- b
  results$r[1]   <- r
  results$r_p[1] <- r
  
  # Initial f assuming catch = surplus
  f_int <- (results$r[1] / p) * (1 - ((results$b[1] / k) ^ p)) 
  

  # Function to calculate fmsy
  fmsy <-function(r,p) {
    r * (1 / (1 + p))
  } 
  
  # f_msy and f_ratio calculations: 
  r_calc1   <- results$r[1]
  r_calc_p1 <- results$r_p[1]
  
  results$f_msy[1]     <- fmsy(r=r_calc1, p=p) 
  results$f_msy_p[1]   <- fmsy(r=r_calc_p1, p=p)  
  results$f_ratio[1]   <- f_int/results$f_msy[1] 
  results$f_ratio_p[1] <- f_int/results$f_msy_p[1]
  
  # f_ratio_err calculation:
  # Log transform to avoid pullin negative values
  mu_1 <- log(results$f_ratio_p[1]) 
  cv   <- error
  sd_1 <- sqrt(log(cv^2 + 1))

  results$f_ratio_err[1] <- rlnorm(1, meanlog = mu_1, sdlog = sd_1)
  
  
  # Decide how to change f based on the f_ratio estimate with error
  if (results$f_ratio_err[1] >= 2) {
    results$f[1] <- hcr * f_int 
  } 
  
  if (results$f_ratio_err[1] > 1.1 & results$f_ratio_err[1] < 2) {
    results$f[1] <- hcr * f_int 
  }  

  if (results$f_ratio_err[1] > 1 & results$f_ratio_err[1] < 1.1) {
    results$f[1] <- f_int  
  }

  if (results$f_ratio_err[1] < 1) {
    results$f[1] <- 1.05 * f_int 
  } 
  
  # Adjust catch if the f_ratio was over the limit to simulate "near closure"
  if (results$f_ratio_err[1] >= 2) {
    results$c[1] <- results$f[1] * results$b[1] * 0.05 
  } 

  if (results$f_ratio_err[1] < 2) {
    results$c[1] <- results$f[1] * results$b[1] 
  }
  
  
# Loop the model over the specified number of years
# Management decisions (i.e. changing fishing pressure (f)) are only made during assessment years
  
  # Write a "not contained in" function" for non-assessment years
  `%not_in%` <- purrr::negate(`%in%`)
  
  for (t in 2:years) {
 
# Assessment year results: 
      if (results$year[t] %in% assess_int) {
      
      # Growth rates
      results$r[t]   <- results$r[t-1] + (r_s * results$r[t-1])
      results$r_p[t] <- results$r_p[t-1] + (r_p_s * results$r_p[t-1])
      
      # f and f_msy  
      r_calc2   <- results$r[t] 
      r_calc_p2 <- results$r_p[t] 
      
      results$f_msy[t]     <- fmsy(r=r_calc2, p=p)
      results$f_msy_p[t]   <- fmsy(r=r_calc_p2, p=p)
      results$f_ratio[t]   <- results$f[t-1]/results$f_msy[t-1] 
      results$f_ratio_p[t] <- results$f[t-1]/results$f_msy_p[t-1] 
      
      # f_ratio_err calculation:
      if (results$f_ratio_p[t] != 0) {
        
        mu_2 <- log(results$f_ratio_p[t])
      
        results$f_ratio_err[t] <- rlnorm(1, meanlog = mu_2, sdlog = sd_1)
      }
      # Use if the fishery is closed
      if (results$f_ratio_p[t] == 0) {
        results$f_ratio_err[t] <- 0
      }
      
      # Decide how to change f based on the f_ratio estimate with error
      if (results$f_ratio_err[t] >= 2) {
        results$f[t] <- hcr * results$f[t-1] 
      } 
      
      if (results$f_ratio_err[t] > 1.1 & results$f_ratio_err[t] < 2) {
        results$f[t] <- hcr * results$f[t-1] 
      }  
       
      if (results$f_ratio_err[t] > 1 & results$f_ratio_err[t] < 1.1) {
        results$f[t] <- results$f[t-1]  
      }
     
      if (results$f_ratio_err[t] < 1) {
        results$f[t] <- 1.05 * results$f[t-1] 
      } 
      
      
      #Calculate the biomass
      results$b[t] <- results$b[t-1] + (results$r[t-1] / p) * 
        results$b[t-1] * (1 - ((results$b[t-1] / k) ^ p)) - results$c[t-1]
      
      # Adjust catch if the f_ratio was over the limit to simulate "near closure"
      if (results$f_ratio_err[t] >= 2) {
        results$c[t] <- results$f[t] * results$b[t] * 0.05 #only keep 5% of the catch 
      } 

      if (results$f_ratio_err[t] < 2) {
        results$c[t] <- results$f[t] * results$b[t]
      }
      
    } 
    
# Non-assessment year results: 
    if (results$year[t] %not_in% assess_int) {
      
      # Growth rates
      results$r[t]   <- results$r[t-1] + (r_s * results$r[t-1])
      results$r_p[t] <- results$r_p[t-1] + (r_p_s * results$r_p[t-1])
      
      # Calculate biomass 
      results$b[t] <- results$b[t-1] + (results$r[t-1] / p) * 
        results$b[t-1] * (1 - ((results$b[t-1] / k) ^ p)) - results$c[t-1]
      
      # Fishing pressure stays the same as prior year
      results$f[t] <- results$f[t-1]
      
      # Calculate catch 
      results$c[t] <- results$f[t] * results$b[t]
      
      # f_msy results need to continue to update annually to capture changes in r and r_p
      r_calc3   <- results$r[t]
      r_calc_p3 <- results$r_p[t]
      
      results$f_msy[t]   <- fmsy(r=r_calc3, p=p) 
      results$f_msy_p[t] <- fmsy(r=r_calc_p3, p=p)
      
      # f_ratio estimates are only updated during assessment years
      results$f_ratio_err[t] <- results$f_ratio_err[t-1]
      results$f_ratio[t]     <- results$f[t-1]/results$f_msy[t-1]
      results$f_ratio_p[t]   <- results$f_ratio_p[t-1]
      
    }
  }
  return(results)
}
