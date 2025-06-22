# Example R-INLA script for zero-inflated coregionalized modelling
# and extreme value analysis following Opitz et al. (2018)

# Libraries
library(INLA)

# ------------------------------------------------------------------
# 1. Data
# Set 'simulate_example' to TRUE to create a synthetic data set for
# illustration. If set to FALSE, provide the path to your own CSV file
# containing the columns listed below.
simulate_example <- TRUE

if (simulate_example) {
    set.seed(123)
    n <- 200
    coords <- cbind(runif(n, 0, 100), runif(n, 0, 100))
    Year <- sample(1995:2024, n, replace = TRUE)
    Month <- sample(1:12, n, replace = TRUE)
    X1 <- 20 + rnorm(n, 0, 3)
    X2 <- 2 + rlnorm(n, 0, 0.3)
    X3 <- runif(n, 0.1, 2)
    field <- rnorm(n)
    mu1 <- exp(1 + 0.1 * X1 - 0.05 * X2 + 0.03 * X3 + 0.1 * field)
    mu2 <- exp(0.5 + 0.05 * X1 - 0.03 * X2 + 0.02 * X3 + 0.1 * field)
    Y1 <- rnbinom(n, mu = mu1, size = 2)
    Y2 <- rnbinom(n, mu = mu2, size = 1)
    Y1[runif(n) < 0.4] <- 0
    Y2[runif(n) < 0.3] <- 0
    mydata <- data.frame(
        Y1 = Y1,
        Y2 = Y2,
        X = coords[, 1],
        Y = coords[, 2],
        X1 = X1,
        X2 = X2,
        X3 = X3,
        Month = Month,
        Year = Year
    )
} else {
    # Data frame 'mydata' should contain the following columns:
    #   Y1  - disease counts
    #   Y2  - mortality counts
    #   X, Y - UTM coordinates
    #   X1  - sea surface temperature (SST)
    #   X2  - maximum wave height
    #   X3  - water turbidity
    #   Month - month of the year
    #   Year  - year from 1995 to 2024
    mydata <- read.csv("path/to/data.csv")
}

coords <- as.matrix(mydata[, c("X", "Y")])
n <- nrow(mydata)

# ------------------------------------------------------------------
# 2. Spatial mesh and SPDE specification
mesh <- inla.mesh.2d(loc = coords, cutoff = 1, max.edge = c(10, 40))
spde <- inla.spde2.pcmatern(mesh,
    prior.range = c(50, 0.5),
    prior.sigma = c(1, 0.01))

# Replicate index for the two responses
replicate <- rep(1:2, each = n)
indexs <- inla.spde.make.index("field", spde$n.spde, n.group = 2)
# Design matrix for spatial field (counts)
coords_rep <- rbind(coords, coords)
A <- inla.spde.make.A(mesh, loc = coords_rep, group = replicate)

# Separate index for the extreme-value models (single group)
index_ext <- inla.spde.make.index("field_ext", spde$n.spde)
A_ext <- inla.spde.make.A(mesh, loc = coords)

# ------------------------------------------------------------------
# 3. Zero-inflated count model (disease and mortality jointly)
response <- c(mydata$Y1, mydata$Y2)

covariates <- data.frame(
    Intercept = 1,
    SST = c(mydata$X1, mydata$X1),
    WaveHeight = c(mydata$X2, mydata$X2),
    Turbidity = c(mydata$X3, mydata$X3),
    Month = c(mydata$Month, mydata$Month),
    Year = c(mydata$Year, mydata$Year))

stack_counts <- inla.stack(
    data = list(y = response),
    A = list(A, 1),
    effects = list(
        c(indexs, list()),
        covariates
    ),
    tag = "counts")

formula_counts <- y ~ -1 + Intercept + SST + WaveHeight + Turbidity +
    f(Month, model = "rw1") +
    f(Year, model = "rw1") +
    f(field, model = spde, group = field.group,
      control.group = list(model = "exchangeable"))

res_counts <- inla(
    formula_counts,
    family = "zeroinflatednbinomial1",
    data = inla.stack.data(stack_counts),
    control.predictor = list(A = inla.stack.A(stack_counts), compute = TRUE))

# ------------------------------------------------------------------
# 4. Extreme value modelling for SST following Opitz et al. (2018)
# Step 1: exceedance probability over threshold q
q_sst <- quantile(mydata$X1, 0.90)
exceed <- as.integer(mydata$X1 > q_sst)

stack_ext1 <- inla.stack(
    data = list(y = exceed),
    A = list(A_ext, 1),
    effects = list(
        c(index_ext, list()),
        data.frame(Month = mydata$Month, Year = mydata$Year)
    ),
    tag = "ext1")

formula_ext1 <- y ~ 1 + f(Month, model = "rw1") + f(Year, model = "rw1") +
    f(field_ext, model = spde)

res_ext1 <- inla(
    formula_ext1, family = "binomial",
    data = inla.stack.data(stack_ext1),
    control.predictor = list(A = inla.stack.A(stack_ext1), compute = TRUE))

# Step 2: magnitude of exceedances (log-Gamma model)
sel <- mydata$X1 > q_sst
excess <- log(mydata$X1[sel] - q_sst)
A_excess <- inla.spde.make.A(mesh, loc = coords[sel, ])

stack_ext2 <- inla.stack(
    data = list(y = excess),
    A = list(A_excess, 1),
    effects = list(
        c(index_ext, list()),
        data.frame(Month = mydata$Month[sel], Year = mydata$Year[sel])
    ),
    tag = "ext2")

formula_ext2 <- y ~ 1 + f(Month, model = "rw1") + f(Year, model = "rw1") +
    f(field_ext, model = spde)

res_ext2 <- inla(
    formula_ext2, family = "gamma",
    data = inla.stack.data(stack_ext2),
    control.predictor = list(A = inla.stack.A(stack_ext2), compute = TRUE))

# Step 3: estimation of high quantile
prob_exc <- res_ext1$summary.fitted.values[, "mean"]
mean_exc <- res_ext2$summary.fitted.values[, "mean"]
shape <- res_ext2$summary.hyperpar["size for y", "mean"]
q95_sst <- q_sst + exp(mean_exc) * qgamma(0.95, shape = shape, rate = 1)

mydata$SST_q95 <- q95_sst[1:n]

# Repeat the extreme value analysis for X2 and X3 if required

# The updated covariate 'SST_q95' can be included in the main model
covariates$SST_q95 <- c(mydata$SST_q95, mydata$SST_q95)

# End of script
