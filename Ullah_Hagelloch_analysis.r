

# Get the Data And Run MCMC -----------------------------------------------
library(epinet)
cov_matrix = HagellochDyadCov
times_matrix = HagellochTimes
colnames(times_matrix)
head(times_matrix)
head(cov_matrix)
mcmcinput <- MCMCcontrol(nsamp = 20000000, thinning = 40, 
                          burnin = 100000, seed = 1,
                         etapropsd =c(rep(0.05, times = 3), 0.005))

priors <- priorcontrol(etaprior = c(1, 3, 1, 3, 1, 3, 1, 3),
                        bprior = c(1, 3), tiprior = c(1, 3), teprior = c(1, 3),
                        keprior = c(1, 10), kiprior = c(1, 10), betapriordist="inversegamma",
                       thetaipriordist="inversegamma", thetaepriordist="inversegamma",
                       kipriordist="inversegamma",kepriordist="inversegamma",
                       etapriordist="gamma")

c1 = cov_matrix[,4]
c2 = cov_matrix[,5]
c3= cov_matrix[,6]


out_inversegamma <- epinet(~c1+c2+c3, 
              epidata = times_matrix, dyadiccovmat = cov_matrix, mcmcinput = mcmcinput,
              priors = priors)
out_inversegamma

# Trace_Plots of Network Parameters ---------------------------------------
coeff = out_inversegamma$eta
par(mfrow=c(2,2))
plot(coeff[,1],col='orchid')
title(main='Trace Plots of Intercept')
plot(coeff[,2],col='orchid')
title(main='Trace Plots of Classroom 1')
plot(coeff[,3],col='orchid')
title(main='Trace Plots of Classroom 2')
plot(coeff[,4],col='orchid')
title(main='Trace Plots of House Distance')


# Summary of the Network Posterior Parmaters --------------------------------------
matrix_statistics = matrix(NA,nrow = 4,ncol = 4)
m1 = mean(coeff[,1][800000:1000000])
sd1 = sd(coeff[,1][800000:1000000])
upper1 = m1+1.96 * (sd1/sqrt(200000))
lower1 = m1-1.96 * (sd1/sqrt(200000))
matrix_statistics[1,] = c(m1,sd1,upper1,lower1)
m2 = mean(coeff[,2][800000:1000000])
sd2 = sd(coeff[,2][800000:1000000])
upper2 = m2+1.96 * (sd2/sqrt(200000))
lower2 = m2-1.96 * (sd2/sqrt(200000))
matrix_statistics[2,] = c(m2,sd2,upper2,lower2)
m3 = mean(coeff[,3])
sd3 = sd(coeff[,3])
upper3 = m3+1.96 * (sd3/sqrt(1000000))
lower3 = m3-1.96 * (sd3/sqrt(1000000))
matrix_statistics[3,] = c(m3,sd3,upper3,lower3)
m4 = mean(coeff[,4][800000:1000000])
sd4 = sd(coeff[,4][800000:1000000])
upper4 = m4+1.96 * (sd4/sqrt(200000))
lower4 = m4-1.96 * (sd4/sqrt(200000))
matrix_statistics[4,] = c(m4,sd4,upper4,lower4)
df = data.frame(matrix_statistics)
`colnames<-`(df,c('Mean','Standard_Deviation','95 % CI Upper Limit','95 % CI lower Limit'))


# Trace Plots for the other Bayesian Parameters ---------------------------

hist(coeff[,4],col='orchid',probability = TRUE)
lines(density(coeff[,4]),col='red',lwd=4)

par(mfrow=c(1,1))
post_thetai = out_inversegamma$thetai
hist(post_thetai,col='orchid',probability = TRUE)
lines(density(post_thetai),col='red',lwd=4)


plot(post_thetai,col='orchid')
title('Posterior_theta_infectiour')

par(mfrow=c(1,1))
post_thetae = out_inversegamma$thetae 
hist(post_thetae,col='orchid',probability = TRUE)
lines(density(post_thetae),col='red',lwd=4)

plot(post_thetae,col='orchid')
title('Posterior_theta_exposure')

par(mfrow=c(1,1))
post_beta = out_inversegamma$beta
hist(post_beta,col='orchid',probability = TRUE)
lines(density(post_beta),col='red',lwd=4)

plot(post_beta,col='orchid')
title('Posterior_beta')


par(mfrow=c(1,1))
post_ki = out_inversegamma$ki
hist(post_ki,col='orchid',probability = TRUE)
lines(density(post_ki),col='red',lwd=4)

plot(post_ki,col='orchid')
title('Posterior_ki')


par(mfrow=c(1,1))
post_ke = out_inversegamma$ke 
hist(post_ke,col='orchid',probability = TRUE)
lines(density(post_ke),col='red',lwd=4)

plot(post_ke,col='orchid')
title('Posterior_ke')


post_iexp =out_inversegamma$initexp
post_iexp[1]



post_exptime=out_inversegamma$initexptime
hist(post_exptime,col='orchid',probability = TRUE)
lines(density(post_exptime),col='red',lwd=4)


plot(post_exptime,col='orchid')
title('Posterior_exposure_time')



# Summary of Bayesian Pramters --------------------------------------------

m1=mean(post_thetai)
sd1 = sd(post_thetai)

c('Mean','SD')
df1 = data.frame(m1,sd1)
`colnames<-`(df1, c('Mean','SD'))
cat('Summary of Posterior theta_infectioous')


m1=mean(post_thetae)
sd1 = sd(post_thetae)

c('Mean','SD')
df1 = data.frame(m1,sd1)
`colnames<-`(df1, c('Mean','SD'))
cat('Summary of Posterior theta_exposure')


m1=mean(post_ke)
sd1 = sd(post_ke)

c('Mean','SD')
df1 = data.frame(m1,sd1)
`colnames<-`(df1, c('Mean','SD'))
cat('Summary of Posterior ke')


m1=mean(post_ki)
sd1 = sd(post_ki)

c('Mean','SD')
df1 = data.frame(m1,sd1)
`colnames<-`(df1, c('Mean','SD'))
cat('Summary of Posterior ki')

m1=mean(post_exptime)
sd1 = sd(post_exptime)

c('Mean','SD')
df1 = data.frame(m1,sd1)
`colnames<-`(df1, c('Mean','SD'))
cat('Summary of Posterior exposure time')


m1=mean(coeff[,4])
sd1 = sd(coeff[,4])

c('Mean','SD')
df1 = data.frame(m1,sd1)
`colnames<-`(df1, c('Mean','SD'))
cat('Summary of House Distance')


# Write the Results in PDF ------------------------------------------------

pdf(file="SampleHagellochTree", height = 20, width = 7)
plot(out, index = 100)
dev.off()

write.epinet(out, "HagellochOutput")

# Auto Correlation for Each Markov Chain ----------------------------------
par(mfrow=c(1,1))
acf(post_beta)


# Estimator for Each Bayesian Parameter Along With Standard Error ---------
library(mcmcse)
cat('Estimator and SE for Beta')
mcse(post_beta)

cat('Estimator and SE for theta_exposure')
mcse(post_thetae)

cat('Estimator and SE for theta_infectious')
c(mcse(post_thetai)[1],mcse(post_thetai)[2])


cat('Estimator and SE for k_infectious')
c(mcse(post_ki)[1],mcse(post_ki)[2])


cat('Estimator and SE for k_exposure')
c(mcse(post_ke)[1],mcse(post_ke)[2])


cat('Estimator and SE for exposure_time')
c(mcse(post_exptime)[1],mcse(post_exptime)[2])

cat('Epidemic Most Probably started by Person Number : ')
post_iexp[1]

cat('Estimate for Intercept alongside SE')
c(mcse(coeff[,1])[1],mcse(coeff[,1])[2])


cat('Estimate for Classroom 1 alongside SE')
c(mcse(coeff[,2])[1],mcse(coeff[,2])[2])

cat('Estimate for Classroom2 alongside SE')
c(mcse(coeff[,3])[1],mcse(coeff[,3])[2])

cat('Estimate for House Distance alongside SE')
c(mcse(coeff[,4])[1],mcse(coeff[,4])[2])





# Bayesian Mean Estimates vs Sample Size ---------------------------------------
estvssamp(post_beta,col='red',lwd=4,main = 'Mean Estimates vs Sample Size for Beta')
estvssamp(post_thetae,col='red',lwd=4,
          main = 'Mean Estimates vs Sample Size for theta_exposure')

estvssamp(post_thetai,col='red',lwd=4,
          main = 'Mean Estimates vs Sample Size for theta_infectious')
estvssamp(post_exptime,col='red',lwd=4,
          main = 'Mean Estimates vs Sample Size for exposure_time')

estvssamp(post_ke,col='red',lwd=4,
          main = 'Mean Estimates vs Sample Size for k_exposure')
estvssamp(post_ki,col='red',lwd=4,
          main = 'Mean Estimates vs Sample Size for k_infectious')

estvssamp(coeff[,1],col='red',lwd=4,
          main = 'Mean Estimates vs Sample Size for Intercept')

estvssamp(coeff[,2],col='red',lwd=4,
          main = 'Mean Estimates vs Sample Size for Classroom1')

estvssamp(coeff[,3],col='red',lwd=4,
          main = 'Mean Estimates vs Sample Size for Classroom2')
estvssamp(coeff[,4],col='red',lwd=4,
          main = 'Mean Estimates vs Sample Size for House Distance')



# Effective Sample Sizes for Each Markov Chain ----------------------------
effective_size = rep(NA,10)
effective_size[1]=ess(as.matrix(post_beta))
effective_size[2]=ess(as.matrix(post_thetae))
effective_size[3]=ess(as.matrix(post_thetai))
effective_size[4]=ess(as.matrix(post_ke))
effective_size[5]=ess(as.matrix(post_ki))
effective_size[6]=ess(as.matrix(post_exptime))
effective_size[7]=ess(as.matrix(coeff[,1]))
effective_size[8]=ess(as.matrix(coeff[,2]))
effective_size[9]=ess(as.matrix(coeff[,3][400000:500000]))
effective_size[10]=ess(as.matrix(coeff[,4]))
plot(effective_size,type = 'b',col='red',lwd=3)







# MCMC Diagnostics --------------------------------------------------------
library(coda)
autocorr(mcmc(post_beta),lags = c(50,100,500))
autocorr(mcmc(post_thetae),lags = c(50,100,500))
autocorr(mcmc(post_thetai),lags = c(50,100,500))
autocorr(mcmc(post_exptime),lags = c(50,100,500))
autocorr(mcmc(post_ki),lags = c(50,100,500))
autocorr(mcmc(post_ke),lags = c(50,100,500))
autocorr(mcmc(coeff[,1]),lags = c(50,100,500))
autocorr(mcmc(coeff[,2]),lags = c(50,100,500))
autocorr(mcmc(coeff[,3]),lags = c(50,100,500))
autocorr(mcmc(coeff[,4]),lags = c(50,100,500))
par(mfrow=c(1,1))

#epidemic parameters diagnostics
heidel.diag(mcmc(post_thetai))
heidel.diag(mcmc(post_ki))
heidel.diag(mcmc(post_thetae))
heidel.diag(mcmc(post_ke))
heidel.diag(mcmc(post_beta)) 
heidel.diag(mcmc(post_exptime))
#Network parameters diagnostics
heidel.diag(mcmc(coeff[,1])) 
heidel.diag(mcmc(coeff[,2]))
heidel.diag(mcmc(coeff[,3])[450000:500000])
heidel.diag(mcmc(coeff[,4])) #Household

geweke.diag(mcmc(post_thetai))
geweke.diag(mcmc(post_ki))
geweke.diag(mcmc(post_thetae))
geweke.diag(mcmc(post_ke))
geweke.diag(mcmc(post_beta)) 
geweke.diag(mcmc(post_exptime))



# HPD CI -------------------------------------------------------------
HPDinterval(mcmc(post_beta))
HPDinterval(mcmc(post_thetae))
HPDinterval(mcmc(post_thetai))
HPDinterval(mcmc(post_ki))
HPDinterval(mcmc(post_ke))
HPDinterval(mcmc(post_exptime))
HPDinterval(mcmc(coeff[,4]))

# Posterior Mean of Infectious and exposed process ------------------------
hist(post_ke*post_thetae,col='orchid',probability = TRUE)
lines(density(post_ke*post_thetae),col='red',lwd=4)

hist(post_ki*post_thetai,col='orchid',probability = TRUE)
lines(density(post_ki*post_thetai),col='red',lwd=4)

# average_number_of_contacts ----------------------------------------------
r0=187 * 0.046 * (1-((1/(1+(post_beta*post_thetai)))^post_ki))
HPDinterval(mcmc(r0))
mean(r0)
sd(r0)
mcse(r0)
hist(r0,col='orchid',probability = TRUE)
lines(density(r0),col='red',lwd=4)
