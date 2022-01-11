#####################################################################################
#     Nonlinear regression of the Baranyi and Roberts model (1994) as adapted by    #
#     Fischer et al. (2014)                                                         #
#####################################################################################
rm(model.fit.LAG)
rm(model.fit.NOLAG)
setwd("C:/path/to/file")
#myPaths <- .libPaths("C:/path/to/set") # set path to files
install.packages("ggplot2")
library(ggplot2)
install.packages("nlstools")
library(nlstools)

# #Function to be fitted is a logistic growth model with a lag-phase
# lag.logistic <- function(n0,#population size at time t = 0
#                          r, #intrinsic growth rate 
#                          K, #carrying capacity
#                          lambda, #length of the lag phase
#                          t #time 
# ){
#    K * n0 * (exp(r * lambda) + exp(r * t) - 1) / (K * exp(r * lambda) + n0 * (exp(r * t) -1 ))
# }

#Function to be fitted is a logistic growth model withOUT a lag-phase
g.logistic <- function(n0,#population size at time t = 0
                         r, #intrinsic growth rate 
                         K, #carrying capacity
                         t #time 
){
  K * n0 *  exp(r * t) / (K  + n0 * (exp(r * t) - 1 ))
}

# #test whether function actually works
# test.times<- seq(0,30, by = 2) 
# test.data <- sapply(X = test.times, FUN =lag.logistic, n0 = 0.05,r = 2.17, K = 0.17,lambda = 0.87)
# 
# ggplot(data = data.frame(t = test.times, n = test.data), aes(x = t, y = n)) + 
#   geom_path()
# 
# #get some test values by adding random noise to data
# ver.data <- rep(test.data , 3) + sapply(0.1*sqrt(rep(test.data,3)), rnorm, n = 1, mean = 0)

# load in my data
rm(Data)
Data = read.csv("file_name.csv", header = T)
Data



#fit model to data- test data
# rm(model.fit) # remove previous model fits from the environment
# model.fit <- nls(n ~ lag.logistic(n0, r, K, lambda, t) , #formula
#                  data = list(t = rep(test.times,3), n = ver.data), #data to be fitted
#                  start = list(n0 = .05, r =2 , K = 5, lambda = 0.5))
# 
# model.fit.LAG <- nls(n ~ lag.logistic(n0, r, K, lambda, t ), #formula
#                  data = list(t = rep(test.times,3), n = ver.data), #data to be fitted
#                  start = list(n0 = .05, r =2 , K = 5, lambda = 0.5)) #start values for the parameters (try different values with an educated guess)
# 
# model.fit.NOLAG <- nls(n ~ g.logistic(n0, r, K,  t), #formula
#                      data = list(t = rep(test.times,3), n = ver.data), #data to be fitted
#                      start = list(n0 = .05, r =2 , K = 5)) #start values for the parameters (try different values with an educated guess)
# 
# 
# fit models to my data
# model.fit.LAG <- nls(n ~ lag.logistic(n0, r, K, lambda, t), #formula
#                       data = Data, #data to be fitted
#                       start = list(n0 = 0.001, r =2 , K = 0.1, lambda = 1.0)) #start values for the parameters (try different values with an educated guess)


model.fit.NOLAG <- nls(n ~ g.logistic(n0, r, K,  t), #formula
                       data = Data, #data to be fitted
                       start = list(n0 = 0.001, r =2 , K = 0.60)) #start values for the parameters (try different values with an educated guess)


ggplot() + 
  geom_point(aes(x = t, y = n), data = Data) + 
  geom_point(aes(t,pn), color = "blue", data = data.frame(t = Data$t, pn = predict(model.fit.NOLAG)[1:length(Data$t)]))+ #prediction of the fitted model
  geom_line(aes(t,pn), color = "blue", data = data.frame(t = Data$t, pn = predict(model.fit.NOLAG)[1:length(Data$t)]))+#prediction of the fitted model
  scale_x_continuous(name = "Time (h)", limits = c(0,9))+
  scale_y_continuous(name = "OD 600nm", limits = c(0, 0.400))+
  ggtitle("Strain Replicate")+
  theme(plot.title = element_text(size=16))+
  theme(axis.text=element_text(size=16, face="bold"))+
  theme(axis.title=element_text(size=16,face="bold"))


#check if the fit is okay
#summary(model.fit.LAG)
summary(model.fit.NOLAG)

#tab_model(model.fit.NOLAG)
# s <- summary(model.fit.NOLAG)
# capture.output(s, file = "myfile.csv")

#Plot test data
# ggplot() + 
#   geom_path(aes(x = t, y = n),data = data.frame(t = test.times,  n= log10(test.data)))+ #actual model values
#   geom_point(aes(x = t, y = nvar), data = data.frame(t = rep(test.times,3),  nvar = log10(ver.data))) + #values with some random noise
#   #geom_point(aes(t,pn), color = "red", data = data.frame(t = test.times, pn = predict(model.fit.LAG)[1:length(test.times)]))+ #prediction of the fitted model
#   geom_point(aes(t,pn), color = "blue", data = data.frame(t = test.times, pn = predict(model.fit.NOLAG)[1:length(test.times)])) #prediction of the fitted model

#check residuals
ggplot() + 
  #geom_point(aes(t,pn), color = "red", data = data.frame(t = rep(test.times,3), pn = residuals(model.fit.LAG)))+ #prediction of the fitted model- test data
  geom_point(aes(t,pn), color = "blue", data = data.frame(t = Data$t, pn = residuals(model.fit.NOLAG))) #prediction of the fitted model- my data


#plot(y = residuals(model.fit.LAG), x = fitted(model.fit.LAG))
plot(y = residuals(model.fit.NOLAG), x = fitted(model.fit.NOLAG))

# qqnorm( residuals(model.fit.LAG))
# qqline(residuals(model.fit.LAG))
qqnorm( residuals(model.fit.NOLAG))
qqline(residuals(model.fit.NOLAG))

#determine Goodness-of-Fit by correlation of the predicted values to the data
#cor(predict(model.fit.LAG), Data$n) #0.993422 (test data)
cor(predict(model.fit.NOLAG), Data$n) #0.9920472 (test data)
#both have  a very high  correlation > 0.99  thus very good fit

##compare models based on AIC (test data)
#AIC(model.fit.LAG)
#AIC(model.fit.NOLAG)
# > AIC(model.fit.LAG)
# [1] 5.387139
# > AIC(model.fit.NOLAG)
# [1] 12.46236
# in this case the model with a lag-phase fits better

