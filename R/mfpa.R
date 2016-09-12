#'Multi-frequential periodogram analysis
#'
#'This function performs multi-frequential periodogram analysis for univariate 
#'temporal or spatial data series collected with equal intervals. Compared with 
#'the traditional periodogram used in spectral analysis, this method can detect 
#'overlapping signals with fractional frequencies. Fitting a joint 
#'polynomial-trigonometric model is achieved by Ordinary Least Squares (OLS) 
#'regression. The function also performs autocorrelation analysis of OLS 
#'residuals up to a number of lags determined by the user.
#'
#'@aliases mfpa plot.mfpa print.mfpa
#'  
#'@param y	Vector of \emph{n} observations (vector of integer or real numbers, or 
#'  one-column matrix).
#'@param MaxNFreq	Maximum number of frequencies to be estimated in the stepwise 
#'  procedure (e.g. 2).
#'@param MinFreq	Minimum value for frequency estimates (e.g. 3.0).
#'@param MaxFreq	Maximum value for frequency estimates (e.g. 10). Must be larger
#'  than MinFreq and smaller than half of the number of observations in the 
#'  series. If unspecified by the user \code{(MaxFreq=NA)}, \code{MaxFreq} is 
#'  set to \emph{n}/4 by the function.
#'@param ntrend	Number (0 to 3) of orthogonal polynomial components estimating 
#'  the broad-scale trend, to be included in the joint polynomial-trigonometric 
#'  model. Use 0 to estimate no trend component, only an intercept.
#'@param nlags	Number of lags to be used for autocorrelation analysis of OLS 
#'  residuals. Use 0 to bypass this analysis.
#'@param alpha	Significance threshold for including frequencies.
#'@param x An object of class \code{mfpa}
#'@param xlab,ylab Labels for x and y axes
#'@param \dots Further arguments passed to or from other methods
#'  
#'@details
#'
#'The fitting of a joint polynomial-trigonometric model is limited to ordinary 
#'least squares (OLS), with autocorrelation analysis of OLS residuals up to a 
#'certain lag. Orthogonal polynomials are used to model broad-scale trends, 
#'whereas cosines and sines model the periodic structures at intermediate 
#'scales. See Dutilleul (2011, section 6.5) and Legendre & Legendre (2012, 
#'section 12.4.4) for details. OLS regression could be replaced by an \emph{estimated 
#'generalized least squares} (EGLS) procedure, as described in Dutilleul (2011). 
#'
#'In spectral analysis in general and in mfpa in particular, the cosines and 
#'sines are considered jointly in the search for the dominant frequency 
#'components since they are both required to fully account for a frequency 
#'component in a linear model. So, when either the cosine or the sine is 
#'significant, this is sufficient indication that a significant frequency 
#'component has been found. But see the first paragraph of the ‘Recommendations 
#'to users’ below.
#'
#'The periodic phenomenon corresponding to each identified 
#'frequency is modelled by a cosine and a sine. The first pair (‘cos 1’, ‘sin 
#'1’) corresponds to the first frequency, the second pair to the second 
#'frequency, and so on. An intercept is also computed, as well as a polynomial 
#'broad-scale trend if argument ntrend > 0. The coefficients shown for each 
#'periodic component (‘cos’ and ‘sin’) are the OLS regression coefficients. The 
#'tests of significance producing the p-values (called ‘prob’ in the output 
#'file) are 2-tailed parametric t-tests, as in standard OLS regression.
#'
#'A global
#'R-square statistic for the periodogram is computed as the variance of the 
#'fitted values divided by the variance of the data series. An R-squared 
#'corresponding to each frequency is also returned.
#'
#'In the Dutilleul 
#'periodogram, the time unit is the length of the data series (in time units: 
#'seconds, hours, days, etc.). Hence, the \emph{frequency} identified by a Dutilleul 
#'periodogram is the number of cycles of the periodic signal (how many full or 
#'partial cycles) along the time series. That number is an integer when the 
#'series contains an integer number of cycles; it may also be a real number when
#'the number of cycles is fractional. The periodogram can identify several 
#'periodic phenomena with different frequencies. The estimated frequencies could
#'be divided by an appropriate constant to produce numbers of cycles per second 
#'or day, or per meter or km, depending on the study.
#'
#'To find the \emph{period} (number
#'of days, hours, etc.) of the process generating a periodic signal in the data,
#'divide the length of the series (in days, hours, etc.) by the frequency 
#'identified by Dutilleul’s periodogram. Recommendations to users The mfpa code 
#'estimates the periodic frequencies to be included in the model through a 
#'combination of a stepwise procedure and non-linear optimisation. Following 
#'that, the contributions of the ‘cos’ and ‘sin’ components of all frequencies 
#'in the model are estimated by multiple linear regression in the presence of 
#'the intercept and trends (if any). Because the mfpa method estimates 
#'fractional frequencies, the cos-sin combinations are not orthogonal among the 
#'identified frequencies, and unnecessary frequencies may be selected as 
#'‘significant’. 
#'
#'1. It is important that users of this periodogram have 
#'hypotheses in mind about the frequencies of the processes that may be 
#'operating on the system under study and the number of periodic components they
#'are expecting to find. If one asks for more components than the number of 
#'periodic phenomena at work on the system, the ‘real’ frequency usually has a 
#'strong or fairly strong R-squared and it is followed by other components with 
#'very small R-squared. Selection of frequencies of interest should thus be 
#'based more upon examination of the R-squares of the components rather than on 
#'the p-values. For short series in particular, the adjusted R-squared is an 
#'unbiased estimate of the variance of the data explained by the model. Even 
#'series of random numbers can produce ‘significant’ frequencies for periodic 
#'components; the associated (adjusted) R-squares will, however, be very small. 
#'
#'2. Function mfpa cannot detect frequencies < 1 (smaller than one cycle in the 
#'series) or larger than (\emph{n}–1) where \emph{n} is the number of observations in the 
#'series, the latter case corresponding to periods smaller than the interval 
#'between successive observations. When a periodic component with such a period 
#'is present in the data, Dutilleul’s periodogram can detect harmonics of that 
#'frequency. Recommendation: when a frequency is detected that does not seem to 
#'correspond to a hypothesized process, one could check, using simulated data, 
#'if it could be produced by a process operating at a temporal scale (period) 
#'smaller than the interval between successive observations.  An example is 
#'shown in Example 2. 
#'
#'3. When analysing a time series with unknown periodic 
#'structure, it is recommended to try first with more than one frequency, say 2 
#'or 3, and also with a trend. Eliminate the non-significant components, step by
#'step, in successive runs, starting with the trend(s), then eliminate the 
#'weakly significant periodic components, until there are only highly 
#'significant components left in the model.
#'
#'@return A list containing the following elements:
#'  
#'  \itemize{ 
#'  \item \code{frequencies}: Vector of estimated frequencies of the model periodic components
#'  and associated R-squared. The frequencies are numbers of cycles in the whole
#'  (temporal or spatial) series under study. 
#'  \item \code{coefficients}: Data frame containing
#'  OLS slope estimates, starting with the intercept, then the orthogonal 
#'  polynomials modelling trend in increasing order, followed by the cosine and 
#'  sine coefficients (alternating) in the order of the estimated frequencies. 
#'  Columns: (1) \code{coefficient}: the OLS intercept or slope estimates; (2) 
#'  \code{prob}: the associated probabilities. 
#'  \item \code{predicted}: A vector (length \emph{n}) of 
#'  predicted response values (fitted values), including the trend if any. The 
#'  data and predicted values can be plotted together using function \code{plot.mfpa}; 
#'  type plot(name.of.output.object). The data values are represented by red 
#'  circles and the fitted values by a black line. 
#'  \item \code{auto_coeff}: If nlags > 0: data
#'  frame containing the following columns. (1) \code{lag}: lags at which 
#'  autocorrelation analysis of the OLS residuals is performed; (2) 
#'  \code{auto_r}: vector of sample autocorrelation coefficients calculated from
#'  OLS residuals for each lag; (3) \code{prob}: vector of probabilities 
#'  associated with the tests of significance of the sample autocorrelation 
#'  coefficients. 
#'  \item \code{y}: the original data series (one-column matrix). 
#'  \item \code{X}: the matrix of explanatory variables; it contains a column of 
#'  “1” to estimate the intercept, a column for each of the trend components (if
#'  any), and two columns for each frequency component, each frequency being 
#'  represented by a cosine and a sine. 
#'  \item \code{r.squared.global}: The global R-squared of
#'  the model and the adjusted R-squared.
#'  }
#'  
#'@references
#'
#'Dutilleul, P. 1990. Apport en analyse spectrale d’un périodogramme modifié et 
#'modélisation des séries chronologiques avec répétitions en vue de leur 
#'comparaison en fréquence. Doctoral Dissertation, Université Catholique de 
#'Louvain, Louvain-la-Neuve, Belgium. Dutilleul, P. 1998. Incorporating scale in
#'ecological experiments: data analysis. Pp. 387-425 in: D. L. Peterson & V. T. 
#'Parker [eds.] Ecological scale – Theory and applications. Columbia University 
#'Press, New York. Dutilleul, P. 2001. Multi-frequential periodogram analysis 
#'and the detection of periodic components in time series. Commun. Stat. - 
#'Theory Methods 30, 1063–1098. Dutilleul, P. R. L. 2011. Spatio-temporal 
#'heterogeneity — Concepts and analyses. Cambridge University Press, Cambridge. 
#'Dutilleul, P. and C. Till. 1992. Evidence of periodicities related to climate 
#'and planetary behaviors in ring-width chronologies of Atlas cedar (Cedrus 
#'atlantica) in Morocco. Can. J. For. Res. 22: 1469-1482. Legendre, P. and P. 
#'Dutilleul. 1992. Introduction to the analysis of periodic phenomena. 11-25 in:
#'M. A. Ali [ed.] Rhythms in fishes. NATO ASI Series, Vol. A-236. Plenum, New 
#'York. Legendre, P. and L. Legendre. 2012. Numerical Ecology. 3rd English 
#'edition. Elsevier, Amsterdam.
#'
#'@author  Guillaume Larocque <glaroc@@gmail.com> and Pierre Legendre.
#'  
#'@examples
#'
#' ### Example 1
#'
#' # Simulate data with frequencies 2.3 and 6.1 and a random component, n = 100. 
#' # No trend, no autocorrelated residuals.
#'
#' y <- as.matrix(0.4*(sin(2.3*2*pi*(1:100)/100)) +
#' 0.4*(sin(6.1*2*pi*(1:100)/100)) + 0.2*rnorm(100))
#'
#' res <- mfpa(y, MaxNFreq = 2, MinFreq = 2, ntrend = 0, nlags = 0)
#'
#' # Compute the periods associated with the two periodic components. Each
#' # frequency in element $frequencies is a number of cycles in the whole series.
#' # The periods are expressed in numbers of time intervals of the data series. In
#' # this example, if the data are measured every min, the periods are in min.
#'
#' periods <- 100/res$frequencies$frequency 
#'
#' # Draw the data series and the fitted (or predicted) values
#'
#' plot(res)
#'
#' ### Example 2
#'
#' # Generate hourly periodic data with tide signal (tide period T = 12.42 h)
#' # during 1 year, hence 24*365 = 8760 hourly data. See
#' # https://en.wikipedia.org/wiki/Tide.
#'
#' # In this simulation, constant (c = 0) puts the maximum value of the cosine at
#' # midnight on the first day of the series.
#'
#' periodic.component <- function(x, T, c) cos((2*pi/T)*(x+c))
#'
#' tide.h <- periodic.component(1:8760, 12.42, 0)
#'
#' # The number of tides in the series is: 8760/12.42 = 705.314 tidal cycles
#' # during one year.
#'
#' # Sample the hourly data series once a day at 12:00 noon every day. The
#' # periodic signal to be detected has a period smaller then the interval between
#' # consecutive observations and its frequency is larger than (n–1). The sequence
#' # of sampling hours for the tide.h data is:
#'
#' h.noon <- seq(12, 8760, 24) 
#' tide.data <- tide.h[h.noon]
#' length(tide.data)   
#'
#' # The series contains 365 sampling units
#'
#' # Compute Dutilleul’s multi-frequential periodogram
#'
#' res.noon <- mfpa(tide.data, MaxNFreq = 1, MinFreq = 2, ntrend = 1, nlags = 2)
#'
#' # Examine the frequency detected by the periodogram, element
#' # res.noon$frequencies. This is a harmonic of the tide signal in the original
#' # data series tide.h.
#'
#' # Compute the period of the signal in the data series sampled hourly:
#'
#' period <- 365/res.noon$frequencies$frequency
#'
#' # Draw the data series and the adjusted values
#'
#' plot(res.noon)
#'
#' # Repeat this analysis after addition of random noise to the tide data
#'
#' tide.noise <- tide.data + rnorm(365, 0, 0.25)
#'
#' res.noise <- mfpa(tide.noise, MaxNFreq = 1, MinFreq = 2, ntrend = 1, nlags = 2)
#'
#' plot(res.noise)
#'
#'@importFrom stats optim
#'@export mfpa
#'  

'mfpa' <-
    function(y,
        MaxNFreq = 2,
        MinFreq = 3,
        MaxFreq = NA,
        ntrend = 0,
        nlags = 0,
        alpha = 0.05)
    {
        #### Internal functions
        
        ## Optimization to fine-tune the frequency estimation
        mfpa_LSoptim <- function(x, y, mf, n, xt) {
            for (i in c(1:mf)) {
                xt <- cbind(xt, cos(2 * pi * x[i] * c(1:n) / n), sin(2 * pi * x[i] * c(1:n) /
                        n))
            }
            out <- -1 * c(t(y) %*% xt %*% solve(t(xt) %*% xt, t(xt) %*% y))
        }
        
        ## Brute force estimation of frequencies
        mfpa_LS <- function(y, X, freqs, mfreqs, ntrend) {
            NF <- length(freqs)
            n <- nrow(y)
            ncX <- ncol(X)
            IM <- matrix(0, NF, 1)
            y <- as.matrix(y)
            for (f1 in c(1:NF)) {
                X[, ncX - 1] <- cos(2 * pi * freqs[f1] * c(1:n) / n)
                X[, ncX] <- sin(2 * pi * freqs[f1] * c(1:n) / n)
                #lower <- max(freqs[mfreqind-1],0), upper = freqs[mfreqind+1])
                IM[f1] <-
                    tryCatch(
                        as.matrix(t(y) %*% X) %*% solve(t(X) %*% X, t(X) %*% y),
                        error = function(e)
                            e = 0
                    )
            }
            m <- max(IM)
            mfreqind <- which.max(IM)
            mf <- length(mfreqs) + 1
            xt <- X[, 1:c(ntrend + 1)]
            ctrl <- list(pgtol = 0, trace = 6)
            mfreq <-
                optim(
                    c(mfreqs, freqs[mfreqind]),
                    mfpa_LSoptim,
                    y = y,
                    mf = mf,
                    n = n,
                    xt = xt,
                    method = "BFGS"
                )
            out <- list(m = m, mfreq = mfreq$par)
        }
        
        ## Compute adjusted R-square from R-square
        adj.R2 <- function(R2, n, m)
            1 - (1 - R2) * (n - 1) / (n - m - 1)
        
        #### END Internal functions
        
        y <- as.matrix(y)
        n <- nrow(y)
        if (is.na(MaxFreq)) {
            # The maximum frequency value (cycles) is equal to one fourth of the
            # time series length.
            MaxF <- n / 4
        } else{
            if (MaxFreq < MinFreq | MaxFreq > (n / 2)) {
                stop(
                    'Invalid maximum frequency (MaxFreq). It must be larger than or equal to
                    MinFreq AND smaller than or equal to half the number of observations in the series'
                )
            } else{
                MaxF = MaxFreq
            }
        }
        res <- 0.5
        freqs <- seq(MinFreq, MaxF, by = res)
        IM <- matrix(0, 2, 1)
        TSS <- sum(y * y)
        alpha <- 0.1
        mfreqs <- NULL
        for (i in 1:MaxNFreq) {
            if (i == 1) {
                # Definition of trend components
                X <- matrix(0, n, 2 + 1 + ntrend)
                if (ntrend > -1)
                    X[, 1] <- matrix(1, n, 1)
                if (ntrend > 0) {
                    t1 <- c(1:n)
                    X[, 2] <- t1 - sum(t1) / n
                }
                if (ntrend > 1) {
                    t2 <- c(1:n) ^ 2
                    n2 <- sum(abs(X[, 2]) ^ 2) ^ (1 / 2)
                    X[, 3] <- t2 - sum(t2) / n - sum(X[, 2] * t2) * X[, 2] /
                        n2 ^ 2
                }
                if (ntrend > 2) {
                    t3 <- c(1:n) ^ 3
                    n3 <- sum(abs(X[, 3]) ^ 2) ^ (1 / 2)
                    X[, 4] = t3 - sum(t3) / n - sum(X[, 2] * t3) * X[, 2] /
                        n2 ^ 2 - sum(X[, 3] * t3) * X[, 3] / n3 ^ 2
                }
                if (ntrend > 3) {
                    stop('The number of orthogonal polynomials must be 3 or less.')
                }
                XP <- X
                yp <- 0
            } else{
                X <-
                    cbind(X, matrix(0, n, 2)) # Frequencies are added in a stepwise procedure.
            }
            LSRes <- mfpa_LS(y, X, freqs, mfreqs, ntrend)
            IM[2] <- LSRes$m
            mfreq <- LSRes$mfreq
            if (i > 1) {
                F <- (0.5 * (IM[2] - IM[1])) / ((TSS - IM[2]) / (n - 2 * i - 1 - ntrend))
                p <- pf(F, 2, (n - 2 * i - 1 - ntrend), lower.tail = FALSE)
                if (p > alpha) {
                    break
                } else{
                    IM[1] <- IM[2]
                    mfreqs <- mfreq
                }
            } else{
                IM[1] <- IM[2]
                mfreqs <- mfreq
            }
            X[, ncol(X) - 1] <-
                cos(2 * pi * round(tail(mfreq, 1) * 10) / 10 * c(1:n) / n)
            X[, ncol(X)] <-
                sin(2 * pi * round(tail(mfreq, 1) * 10) / 10 * c(1:n) / n)
        }
        X <- X[, !colSums(abs(X)) == 0]
        A <- solve(t(X) %*% X, t(X))
        b <- A %*% y # Regression slopes
        yp <- X %*% b # Predicted values
        f.r.squared = matrix(NA, length(mfreqs))
        for (i in 1:length(mfreqs)) {
            ind <- c(ntrend + (2 * i), ntrend + 1 + (2 * i))
            Af <- solve(t(X[, -ind]) %*% X[, -ind], t(X[, -ind]))
            bf <- Af %*% y # Regression slopes
            ypf <- X[, -ind] %*% bf # Predicted values
            f.r.squared[i] <-
                (crossprod(yp) - crossprod(ypf)) / crossprod(y)
        }
        yres <- y - yp # Residuals
        MSE <- sum(yres ^ 2) / (n - qr(X)$rank)
        VarB <- diag(A %*% t(A) * MSE)
        S <- sqrt(VarB)
        tstat <- b / S
        Probb <- (1 - pt(abs(tstat), n - ncol(X))) * 2
        r.squared.global <- vector(length = 2)
        r.squared.global[1] <- Rsq <- var(yp) / var(y)
        r.squared.global[2] <- adj.R2(Rsq, n, (ncol(X) - 1))
        
        if (nlags > 0) {
            # Definition of lags
            lags <- c(1:nlags)
            nlags <- length(lags)
            # Autocorrelation analysis of OLS residuals
            HH <-
                dist(cbind(c(1:n), matrix(1, n, 1)), upper = TRUE, diag = TRUE)
            yMat <- kronecker(matrix(1, 1, n), yres) - mean(yres)
            MSE <- sum(yres * yres) / (n - qr(X)$rank)
            if (nlags != 0) {
                r <- matrix(0, nlags - 1, 1)
                S2 <- matrix(0, nlags - 1, 1)
                yProd <- yMat * t(yMat)
                for (i in 1:nlags) {
                    islag <- as.vector(as.matrix(HH)) == lags[i]
                    tMat <- as.vector(yProd) * islag
                    r[i] <- sum(tMat) / (2 * (n - 1) * var(yres))
                    S2[i] <- (1 + 2 * sum(r[1:(i - 1)] ^ 2)) / n
                }
                Z <- r / sqrt(S2)
                Probr <- 2 * (1 - pnorm(abs(Z), mean = 0, sd = 1))
            } else{
                Z <- 0
            }
        } else{
            lags <- numeric()
            r <- numeric()
            Probr <- numeric()
        }
        
        cresults <- data.frame(coefficient = b, prob = round(Probb, 5))
        nsc <- (ncol(X) - ntrend - 1) / 2
        if (ntrend != 0) {
            rt <- paste(rep("trend order", ntrend), 1:ntrend)
        } else{
            rt <- {
            }
        }
        row.names(cresults) <-
            c('intercept', rt, paste(rep(c('cos', 'sin'), nsc), kronecker(1:(nsc), c(1, 1))))
        
        freqresults <-
            data.frame(frequency = mfreqs, r.squared = f.r.squared)
        
        auto_coeff <- data.frame(lag = lags,
            auto_r = r,
            prob = Probr)
        
        out <-
            list(
                frequencies = freqresults,
                coefficients = cresults,
                predicted = yp,
                auto_coeff = auto_coeff,
                y = y,
                X = X,
                r.squared.global = r.squared.global
            )
        class(out) <- c("mfpa", class(res))
        return(out)
        }

#' @rdname mfpa
#' @export
'plot.mfpa' <-
    function (x,
        xlab = "" ,
        ylab = "Values" ,
        ...) {
        plot(
            x$y,
            type = "p",
            col = "red",
            ylim = range(x$y) * 1.1,
            ylab = ylab,
            xlab = xlab
        )
        lines(x$predicted)
    }

#' @rdname mfpa
#' @export
'print.mfpa' <-
    function (x, ...) {
        cat("Multi-frequential periodogram analysis\n\n")
        print(x$frequencies)
        cat("\n")
        cat(paste0(
            "Global R-squared   : ",
            format(x$r.squared.global[1], digits = 5)
        ), "\n")
        cat(paste0(
            "Adjusted R-squared : ",
            format(x$r.squared.global[2], digits = 5)
        ), "\n")
        cat("\n")
        print(x$coefficients)
        if (length(x$lags) > 0) {
            cat("\n")
            print(x$auto_coeff)
        }
    }
