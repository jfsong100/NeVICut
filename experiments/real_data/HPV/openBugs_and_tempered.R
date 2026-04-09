### R code from vignette source 

rm(list = ls()) # clear environment
current_path <- rstudioapi::getSourceEditorContext()$path
base_dir <- dirname(current_path)
output.path <- file.path(base_dir, "downstream")

###################################################
### code chunk number 1: Load Packages
###################################################
library(BRugs)
library(coda)
library(RColorBrewer)


###################################################
### code chunk number 2: Data
###################################################
nhpv <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,
               143, 229, 696, 93)
ncases <- c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066,
          26751, 75815, 150302, 354993, 3683043, 507218)


###################################################
### code chunk number 3: Full Model
###################################################
fullmodel <- c(
"model {",
"",
"   ## Disease model",
"   for (j in 1:13) {",
"      ncases[j] ~ dpois(mean[j])",
"      log(mean[j]) <- theta1 + phi[j] * theta2 + log(Npop[j] * 1.0E-3) ",
"   }",
"   ",
"   theta1 ~ dnorm(0, 1.0E-3)",
"   theta2 ~ dnorm(0, 1.0E-3) I(0,)",
"",
"   ## Measurement model",
"   for (j in 1:13) {",
"      nhpv[j] ~ dbin(phi[j], Npart[j])",
"   }",
"",
"   ## Exposure model",
"   for (j in 1:13) {",
"      phi[j] ~ dunif(0, 1)",
"   }",
"",
"}")
writeLines(fullmodel, "FullModel.bug")


###################################################
### code chunk number 4: Cut Model
###################################################
cutmodel <- c(
"model {",
"",
"   ## Disease model",
"   for (j in 1:13) {",
"      ncases[j] ~ dpois(mean[j])",
"      log(mean[j]) <- theta1 + p.cut[j] * theta2 + log(Npop[j] * 1.0E-3) ",
"   }",
"",
"   theta1 ~ dnorm(0, 1.0E-3)",
"   theta2 ~ dnorm(0, 1.0E-3)",
"",
"   ## Cuts",
"   for (j in 1:13) {",
"      p.cut[j] <- cut(phi[j])",
"   }",
"",
"   ## Measurement model - below the cut",
"   for (j in 1:13) {",
"      nhpv[j] ~ dbin(phi[j], Npart[j])",
"   }",
"",
"   ## Exposure model - below the cut",
"   for (j in 1:13) {",
"      phi[j] ~ dunif(0, 1)",
"   }",
"",
"}")
writeLines(cutmodel, "CutModel.bug")


###################################################
### code chunk number 5: fitting function
###################################################
fit <- function(model, N.iterations, thin)
{
    data <- c("nhpv","Npart","ncases","Npop")
    modelCheck(model)
    modelData(bugsData(data))
    modelCompile(1)
    modelGenInits()
    modelUpdate(N.iterations, thin=thin)
    samplesSet(c("theta1", "theta2"))  # monitor both theta1 and theta2
    modelUpdate(N.iterations, thin=thin)
    #buildMCMC(c("theta1", "theta2"))   # build chains for both
}


###################################################
### code chunk number 6: Disable Samplers
###################################################
modelDisable('logit/log-linear block glm')
modelDisable('logit/log-linear rejection')
modelDisable('adaptive metropolis 1D')
modelDisable('adaptive acceptance rate')
modelDisable('state dependent scale')
modelDisable('slice')


###################################################
### code chunk number 7: Draw samples from cut model
###################################################
time1 = Sys.time()
updaters <- c('logit/log-linear block glm',
              'logit/log-linear rejection',
              'adaptive metropolis 1D')
cut.samples <- vector("list", 3)
for (i in 1:3) {
    modelEnable(updaters[i])
    modelEnable(updaters[i])
    fit("CutModel.bug", N.iterations=100000, thin=20)
    
    # Extract samples of theta1 and theta2
    theta1_samp <- samplesSample("theta1")
    theta2_samp <- samplesSample("theta2")
    
    # Store as a list of data.frames for easy access
    cut.samples[[i]] <- list(theta1 = theta1_samp, theta2 = theta2_samp)
    #cut.samples[[i]] <- fit("CutModel.bug", N.iterations=1000, thin=20)
    modelDisable(updaters[i])
}
time2 = Sys.time()
naive_cut_time = difftime(time2, time1, units = 'secs')

###################################################
### code chunk number 8: Plot Function
###################################################
plot.samples <- function(s, n_chains=1)
{
    ymax <- sapply(s, function(x) { max(density(unlist(x))$y) })
    xrange <- sapply(s, function(x) { quantile(unlist(x), c(0.002, 0.998)) })

    plot(NA, xlim=range(xrange), ylim=c(0, max(ymax)), ylab="density", 
         xlab=expression(theta[2]))
    pal <- brewer.pal(3, name="Dark2")
    for (i in 1:3) {
        for (j in 1:n_chains) {
            lines(density(s[[i]][[j]]), col=pal[i])
        }
    }

    legend("topright", col=pal, lty=1, lwd=2,
           legend=c("log-linear block glm",
           "log-linear rejection",
           "adaptive metropolis 1D"))
}


###################################################
### code chunk number 9: figA1
###################################################
extract_param <- function(cut_samples, param_name, n_chains = 1) {
    lapply(cut_samples, function(updater) {
        param_vec <- updater[[param_name]]
        n_samples_per_chain <- length(param_vec) / n_chains
        # Split into list of chains
        split(param_vec, rep(1:n_chains, each = n_samples_per_chain))
    })
}




theta1_list <- extract_param(cut.samples, "theta1")
names(theta1_list) <- updaters
plot.samples(theta1_list)
theta2_list <- extract_param(cut.samples, "theta2")
names(theta2_list) <- updaters
plot.samples(theta2_list)

theta_samples_naive_cut = list(theta21_samples = theta1_list,
                             theta22_samples = theta2_list,
                             time = naive_cut_time)
saveRDS(theta_samples_naive_cut, file = file.path(output.path, "theta2_samples_naive_cut.rds"))




###################################################
### JS: get theta1 and theta2 samples
###################################################

theta1_samples <- samplesSample("theta1")
theta2_samples <- samplesSample("theta2")

###################################################
### code chunk number 10: Conjugate updater for p
###################################################
update.p <- function(nhpv, Npart) {

    shape1 <- 1 + nhpv
    shape2 <- 1 + (Npart - nhpv)

    rbeta(length(shape1), shape1, shape2)
}


###################################################
### code chunk number 11: Poisson log likelihood
###################################################
poisson.loglik <- function(Npop, ncases, theta1, theta2, p)
{
    mean <- Npop * 1.0E-3 * exp(theta1 + p * theta2)
    sum(dpois(ncases, mean, log=TRUE))
}


###################################################
### code chunk number 12: Updater for theta1
###################################################
update.theta1 <- function(Npop, ncases, theta1, theta2, p, STEP)
{
    logdensity <- function(a) {
        poisson.loglik(Npop, ncases, a, theta2, p)
    }

    logdensity0 <- logdensity(theta1)
    theta1.new <- theta1 + rnorm(1, mean=0, sd=STEP)
    logdensity1 <- logdensity(theta1.new)

    R <- exp(logdensity1 - logdensity0)
    if (runif(1) < R) theta1.new else theta1
}


###################################################
### code chunk number 13: Updater for theta2
###################################################
update.theta2 <- function(Npop, ncases, theta1, theta2, p, STEP)
{
    logdensity <- function(b) {
        poisson.loglik(Npop, ncases, theta1, b, p)
    }

    logdensity0 <- logdensity(theta2)
    theta2.new <- theta2 + rnorm(1, mean=0, sd=STEP)
    logdensity1 <- logdensity(theta2.new)

    R <- exp(logdensity1 - logdensity0)
    if (runif(1) < R) theta2.new else theta2
}


###################################################
### code chunk number 14: Tempered
###################################################
sample.temper <- function(niter, step.theta1=0.1, step.theta2=1, ntemper=1)
{
    theta1 <- rnorm(1, mean=-1.66, sd=0.344)
    theta2 <- rnorm(1, mean=13.05, sd=3.39)
    p <- rep(13, 0.5)

    theta1.samples <- theta2.samples <- numeric(niter)

    for (i in 1:niter) {

        p.old <- p
        p.new <- update.p(nhpv, Npart)

        for (t in 1:ntemper) {
            pt <- ((ntemper - t) * p.old + t * p.new) / ntemper
            theta1 <- update.theta1(Npop, ncases, theta1, theta2, pt, step.theta1)
            theta2 <- update.theta2(Npop, ncases, theta1, theta2, pt, step.theta2)
        }
        p <- p.new

        theta1.samples[i] <- theta1
        theta2.samples[i] <- theta2
    }

    list("theta1" = theta1.samples, "theta2" = theta2.samples)
}


###################################################
### code chunk number 15: Function to run tempered sampler
###################################################
run.temper <- function(nsim=1000, nburnin=1000, seed=121876) {

    set.seed(seed)

    Nstep <- 8

    #Based on pilot runs, we set thinning interval to get approximately
    #independent samples from each chain (and hence the same effective
    #sample size).
    thin <- c(50, 30, 12, 6, 3, 2, 2, 1)

    theta1.samples <- theta2.samples <- vector("list", Nstep)
    for (t in 1:Nstep) {
        ntemper <- 2^(t-1)
        new.samp <- sample.temper(thin[t] * (nsim + nburnin), ntemper=ntemper)

        theta1 <- as.mcmc(new.samp$theta1)
        theta1 <- window(theta1, start=nburnin * thin[t] + 1, thin=thin[t])
        theta1.samples[[t]] <- theta1
        
        theta2 <- as.mcmc(new.samp$theta2)
        theta2 <- window(theta2, start=nburnin * thin[t] + 1, thin=thin[t])
        theta2.samples[[t]] <- theta2
    }

    list(theta1=theta1.samples, theta2=theta2.samples)
}


###################################################
### code chunk number 16: Plot function for tempered sampler
###################################################
plot.temper <- function(s) {

    plot(density(s[[1]], bw="nrd"), type = "n", main="", 
         xlab=expression(theta[2]))
    Nstep <- length(s)
    for (t in 1:Nstep) {
        lines(density(s[[t]], bw="nrd"), col=grey(1 - t/Nstep), lwd=2)
    }
    legend("topright", col=grey(1 - (1:Nstep)/Nstep), lty=1, lwd=2,
           legend=formatC(2^(0:(Nstep-1)), width=2), title="Steps")
}


###################################################
### code chunk number 17: Generate tempered samples
###################################################
tempered.samples <- run.temper()
plot.temper(tempered.samples$theta2)


###################################################
### code chunk number 18: figA2
###################################################
plot.temper(tempered.samples$theta2)


###################################################
### JS: tempered theta1 vs theta2 samples
###################################################
time1 = Sys.time()
# Run a single chain of the tempered sampler (for example, using ntemper=1)
set.seed(123)
result <- sample.temper(niter = 100000, step.theta1 = 0.1, step.theta2 = 1, ntemper = 128)
time2 = Sys.time()
plummer_time = difftime(time2, time1, units = 'secs')

# Extract theta1 and theta2 samples
theta1_samples <- result$theta1
theta2_samples <- result$theta2

# Create a scatter plot of theta1 versus theta2
plot(theta1_samples, theta2_samples,
     xlab = expression(theta[1]),
     ylab = expression(theta[2]),
     main = "Scatter Plot of theta1 vs theta2",
     pch = 16, col = rgb(0, 0, 1, 0.5))  # pch and color for better visualization

# Plot density for theta1 and theta2 in separate panels
par(mfrow = c(1, 2))
plot(density(theta1_samples), main = expression("Density of " * theta[1] * " (ntemper = 128)"),
     xlab = expression(theta[1]), col = "blue", lwd = 2)
plot(density(theta2_samples), main = expression("Density of " * theta[2] * " (ntemper = 128)"),
     xlab = expression(theta[2]), col = "red", lwd = 2)
par(mfrow = c(1, 1))

theta_samples_plummer = list(theta21_samples = theta1_samples,
                             theta22_samples = theta2_samples,
                             time = plummer_time)
saveRDS(theta_samples_plummer, file = file.path(output.path, "theta2_samples_plummer.rds")
)





