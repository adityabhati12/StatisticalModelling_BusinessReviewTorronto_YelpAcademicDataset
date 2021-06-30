#Installing Packages
library(reshape)
library(tidyr)
library(jsonlite)
library(dplyr)
library(Amelia)
library(mice)
library(lattice)
library(ggplot2)
library(MCMCpack)
library(mcmc)
library(tidyverse)
########################################################################################

YelpMain <- readRDS("Statistical Modeling/review_business_merge.rds")

YelpBusiness <- readRDS("Statistical Modeling/abcd.rds")

###########################################################################################################

# Merging datasets

#business dataset inspection
summary(YelpMain)
str(YelpMain)
names(business)
names(review)
##########################################################################################################

#CHECKING FOR NA VLAUES
sum(is.na(YelpBusiness$is_open))
missmap(YelpBusiness,main="Yelp Data - Missings Map", col=c("black", "red"), legend=TRUE) 
#insignificant Na values in the data so we ignore it
#tidying category feature in business_torronto
###########################################################################################################

# business_torronto%>%select(-starts_with("hours"), -starts_with("attribute")) %>% unnest(categories) %>%
#   select(name, categories)%>%group_by(categories)%>%summarise(n=n())%>%arrange(desc(n))%>%head(20)

###########################################################################################################

#splitting the category
# class(business_torronto$categories)
# cat <- cbind(business_torronto[1, ]$business_id, business_torronto[1, ]$categories[[1]])
# for(i in 2:nrow(business_torronto)) cat <- rbind(cat, cbind(business_torronto[i, ]$business_id, business_torronto[i, ]$categories[[1]]))
# 
# cat_names <- names(sort(table(cat[, 2]), decreasing = TRUE))[1:10]
# 
# cat_bus_ind_mat <- t(sapply(tapply(cat[, 2], cat[, 1], function(y) cat_names %in% y), function(x) x*1))
# 
# colnames(cat_bus_ind_mat) <- cat_names
# df_cat <- data.frame(ind = rownames(cat_bus_ind_mat), cat_bus_ind_mat)
# business_merge <- merge(business_torronto, df_cat, by.x = "business_id", by.y = "ind")
# main_review <- merge(review, business_merge , by.x = "business_id", by.y = "business_id")

############################################################################################################

#CHECKING FOR NA VALUES AND IMPUTING IT
colnames(main_review)[colSums(is.na(main_review)) > 0]
# #IMPUTING attributes.RestaurantsPriceRange2 WITH 1792 NA VALUES
# main_review12 <- main_review
# main_review12 <- matrix(main_review12$attributes.RestaurantsPriceRange2)
# tempData <- mice(main_review$attributes.RestaurantsPriceRange2,m=5,maxit=50,meth='pmm',seed=500)


#############################################################################################################

#converting to factor and numeric
YelpBusiness$neighborhoodID <- as.numeric(YelpBusiness$neighborhood)
YelpBusiness$neighborhoodID <- as.factor(YelpBusiness$neighborhoodID)

business.agg <- YelpBusiness %>% group_by(YelpBusiness$neighborhood, YelpBusiness$neighborhoodID) %>% summarise(Total = n())

YelpBusiness$business.stars <- as.numeric(YelpBusiness$business.stars)

#YelpBusiness$neighborhood <- as.numeric(YelpBusiness$neighborhood)



#Comparing Differnt neighbourhoods without sampling

ggplot(YelpBusiness) + geom_boxplot(aes(YelpBusiness$neighborhood, YelpBusiness$business.stars, fill = YelpBusiness$neighborhood),show.legend=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

sort(tapply(as.numeric(YelpBusiness$business.stars), YelpBusiness$neighborhood, mean), decreasing = TRUE)

sort(tapply(as.numeric(YelpBusiness$business.stars), YelpBusiness$neighborhood, median), decreasing = TRUE)

sort(tapply(as.numeric(YelpBusiness$business.stars), YelpBusiness$neighborhood, sd), decreasing = TRUE)

nbr <- sort(tapply(as.numeric(YelpBusiness$business.stars), YelpBusiness$neighborhood, mean), decreasing = TRUE)[1:2]

YelpBusiness.nbr <- subset(YelpBusiness, YelpBusiness$neighborhood %in% names(nbr))

ggplot(YelpBusiness, aes(business.stars)) + stat_bin()

t.test(YelpBusiness.nbr$business.stars~ YelpBusiness.nbr$neighborhood , data = YelpBusiness.nbr, var.equal =TRUE)


###################################################################################################################
#Comparing different neighbours before Gibbs sampling

  ggplot(YelpBusiness) + geom_boxplot(aes(x = reorder(neighborhood, business.stars, median), business.stars, fill = reorder(neighborhood, business.stars, median)), show.legend=FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


  ggplot(YelpBusiness, aes(business.stars))+ stat_bin()

ggplot(data.frame(size = tapply(YelpBusiness$business.stars, YelpBusiness$neighborhood, length), mean_score = tapply(YelpBusiness$business.stars, YelpBusiness$neighborhood, mean)), aes(size, mean_score)) + geom_point()

####################################################################################################################
#gibbs Sampling
compare_m_gibbs <- function(y, ind, maxiter = 5000)
  {
  
  ### weakly informative priors
  a0 <- 1/2 ; b0 <- 2 ## tau_w hyperparameters
  eta0 <-1/2 ; t0 <- 2 ## tau_b hyperparameters
  mu0<-3 ; gamma0 <- 1
  ### 
  m <- nlevels(ind)
  ybar <- theta <- tapply(y, ind, mean)
  tau_w <- mean(1 / tapply(y, ind, var ),na.rm=TRUE) ##within group precision
  mu <- mean(theta)
  tau_b <-var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  an <- a0 + sum(n_m)/2
  ###
  
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  
  ### MCMC algorithm
  for(s in 1:maxiter) 
  {
    
    # sample new values of the thetas
    for(j in 1:m) 
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    bn <- b0 + ss/2
    tau_w <- rgamma(1, an, bn)
    
    #sample a new value of mu
    gammam <- m * tau_b + gamma0
    mum <- (mean(theta) * m * tau_b + mu0 * gamma0) / gammam
    mu <- rnorm(1, mum, 1/ sqrt(gammam)) 
    
    # sample a new value of tau_b
    etam <- eta0 + m/2
    tm <- t0 + sum((theta-mu)^2)/2
    tau_b <- rgamma(1, etam, tm)
    
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}

fit <- compare_m_gibbs(YelpBusiness$business.stars, YelpBusiness$neighborhoodID)


apply(fit$params, 2, mean)
apply(fit$params, 2, sd)

raftery.diag(as.mcmc(fit))mean(1/sqrt(fit$params[,3]))
plot(as.mcmc(fit$params))


theta_hat <- apply(fit$theta, 2, mean)
a <- ggplot(data.frame(size = tapply(YelpBusiness$business.stars, YelpBusiness$neighborhood, length), theta_hat = theta_hat), aes(size, theta_hat)) + geom_point()
plot(a)

##############################################################################################################################################################
#plot of the the clusters
datasim=as.data.frame(fit$theta)
datasim


colnames(datasim) <- business.agg$`YelpBusiness$neighborhood`
datasim$ID <- seq.int(nrow(datasim))


md <- melt(datasim,id="ID")

ggplot(md) + 
  geom_boxplot(aes(x = reorder(variable, value, mean), value, fill = reorder(variable, value, mean)), show.legend=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

###################################################################################################################
#Comparing Two neighbourhoods 
compare_2_gibbs <- function(y, ind, mu0 = 50, tau0 = 1/625, del0 = 0, gamma0 = 1/625, a0 = 0.5, b0 = 50, maxiter = 5000)
{
  y1 <- y[ind == 1]
  y2 <- y[ind == 2]
  
  n1 <- length(y1) 
  n2 <- length(y2)
  
  ##### starting values
  mu <- (mean(y1) + mean(y2)) / 2
  del <- (mean(y1) - mean(y2)) / 2
  
  mat_store <- matrix(0, nrow = maxiter, ncol = 3)
  #####
  
  ##### Gibbs sampler
  an <- a0 + (n1 + n2)/2
  
  for(s in 1 : maxiter) 
  {
    
    ##update tau
    bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
    tau <- rgamma(1, an, bn)
    ##
    
    ##update mu
    taun <-  tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    mu <- rnorm(1, mun, sqrt(1/taun))
    ##
    
    ##update del
    gamman <-  gamma0 + tau*(n1 + n2)
    deln <- ( del0 * gamma0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  return(mat_store)
}

comp_n <- filter(md, md$variable == "City Place" | md$variable == "Ossington Strip")
fit_new <- compare_2_gibbs(comp_n$value, comp_n$variable)

plot(as.mcmc(fit_new))

raftery.diag(as.mcmc(fit_new))

y1_sim <- rnorm(5000, fit_new[, 1] + fit_new[, 2], sd = 1/sqrt(fit_new[, 3]))
y2_sim <- rnorm(5000, fit_new[, 1] - fit_new[, 2], sd = 1/sqrt(fit_new[, 3]))

ggplot(data.frame(y_sim_diff = y1_sim - y2_sim)) + stat_bin(aes(y_sim_diff))

mean(y1_sim > y2_sim)


ggplot(comp_n) + geom_boxplot(aes(variable, value, fill = variable)) + geom_jitter(aes(variable, value, shape = comp_n$variable))

################################################################################################################################
QUESTION - 2
YelpMain <- read_rds("D:/aditya/TRINITY/Statistical Modeling/reviews_business_merge.rds")

data1 <- YelpMain





fit_mcmc <- MCMCregress(data1$stars.x ~data1$useful+data1$funny+data1$cool+data1$review_count+data1$Food+data1$Nightlife+data1$Bars+
                          data1$Sandwiches+data1$Breakfast...Brunch+data1$Chinese+data1$Canadian..New.+data1$Cafes+data1$Coffee...Tea + data1$neighborhood + data1$stars.y)

fit_mcmc1 <-MCMCregress(data1$stars.x ~data1$useful+data1$funny+data1$cool+data1$review_count+data1$Food+data1$Nightlife+data1$Bars+
                            data1$Sandwiches+data1$Breakfast...Brunch+data1$Canadian..New.+data1$Cafes+ data1$stars.y,B0 = 1e-2, marginal.likelihood = "Chib95")
attr(fit_mcmc1, "logmarglike")
#left
plot(fit_mcmc1)
##


fit_mcmc2 <-  MCMCregress(data1$stars.x ~data1$useful+data1$funny+data1$cool+data1$review_count+data1$Food+data1$Nightlife+data1$Bars+
                            data1$Sandwiches+data1$Breakfast...Brunch+data1$Canadian..New.+data1$Cafes, B0 = 1e-5, marginal.likelihood = "Chib95")
attr(fit_mcmc2, "logmarglike") 

dummy <- data1
dummy <- dummy[ -c(1,2:5) ]
dummy
drops=c("review_count" ,
        "useful" ,
        "cool" ,
        "PriceRange" ,
        "Canadian..New." ,
        "Chinese" ,
        "Coffee...Tea" ,
        "Sandwiches" ,
        "Nightlife" ,
        "Bars" ,
        "Cafes" ,
        "Breakfast...Brunch" ,
        "Food" ,
        "funny",
        "stars.y" )
dummy=dummy[,(names(dummy) %in% drops)]

fit_mcmc <- MCMCregress(data1$business.x ~data1$useful+data1$funny+data1$cool+data1$review_count+data1$Food+data1$Nightlife+data1$Bars+data1$Sandwiches+data1$Breakfast...Brunch+data1$Chinese+data1$Canadian..New.+data1$Cafes+data1$Coffee...Tea + data1$neighborhood + data1$stars.y)

fit_mcmc
par(mar=c(2,2,2,2))
summary(fit_mcmc)
plot(fit_mcmc)

beta_mean <- apply(fit_mcmc1, 2, mean)
coeff <-as.data.frame(beta_mean)
coeff[-1] <- coeff[order(coeff$beta_mean),-1]
data1$attributes.RestaurantsPriceRange2 <- as.numeric(data1$attributes.RestaurantsPriceRange2)

df_dummy <- data1[, -4] 
  # remove response variable
#df_dummy$sex <- as.numeric(df_dummy$sex) # for convenience, change class of sex variable - trickier if variable has several levels

factor_columns <- c("useful", "funny", "cool", "review_count", "Food", "Nightlife", "Bars", "Sandwiches", "Breakfast...Brunch", "Chinese", "Canadian..New.", "Cafes", "Coffee...Tea")
strdf_dummy[factor_columns] <- lapply(strdf_dummy[factor_columns], numeric)




strdf_dummy <- cbind(1, as.matrix(df_dummy)) # add dummy variable for intercept, convert to matrix class
foo <-beta_mean[-15]
pred_fit <- as.factor(strdf_dummy) %*% as.matrix(foo) # make prediction, ignore variance parameter
p <-ggplot(data1$stars.y,pred_fit ) 
data1$predicted <-pred_fit
p + geom_jitter()
plot_external <- ggplot(data1)+
  aes(x = data1$stars.y)+
  # fitted values
  geom_jitter(aes(y=predicted), color = "purple")+
  geom_point(aes(y=stars.y), color = "red")
plot_external




#####################################################################################################
#linear Regression

 fit_lm <- lm(data1$stars.x ~data1$useful+data1$funny+data1$cool+data1$review_count+data1$Food+data1$Nightlife+data1$Bars+data1$Sandwiches+data1$Breakfast...Brunch+data1$Chinese+data1$Canadian..New.+data1$Cafes+data1$Coffee...Tea+ data1$neighborhood)
 library(caret)
 varImp(fit_lm, scale = F)
 
 step_aic <- step(fit_lm)

 step_BIC <- step(fit_lm, k=log(nrow(data1))) 
 
 
 #after removie 
 fit_lm1 <- lm(data1$stars.x ~data1$useful+data1$funny+data1$cool+data1$review_count+data1$Food+data1$Nightlife+data1$Bars+data1$Sandwiches+data1$Breakfast...Brunch+data1$Chinese+data1$Cafes+data1$Coffee...Tea+ data1$neighborhood)
 varImp(fit_lm1, scale = F)
 summry_lm <- summary(fit_lm1)
 rmse <- function(summry_lm)mean(summry_lm$residuals^2)
 
  plot(predict(fit_lm1),data1$stars.x,
      xlab="predicted",ylab="actual")
 abline(a=0,b=1)
 
 ####################################################################################################################
 #question 3
 
 bus_cat <- readRDS("D:/aditya/TRINITY/Statistical Modeling/bus_cat.rds")
 
 
 ad <- bus_cat[,2:30]
 
 md <- melt(bus_cat)
 
 ggplot(md) + 
   geom_point(aes(x = reorder(md$neighborhood), value, fill = reorder(variable, value)), show.legend=FALSE)+
   theme(axis.text.x = element_text(angle = 90, hjust = 1))
 
 
 
 ggplot(md, aes(x=md$neighborhood))
 
 
 anova <-  aov(YelpBusiness$business.stars~YelpBusiness$neighborhood)
 summary(anova)
 TukeyHSD(anova)
 
 
 
 
