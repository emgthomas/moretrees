# ------------------------------------------------------------------------------------- #
# ---------------- Minor edit of logit.spike from BoomSpikeSlab package --------------- #
# ------------------------------------------------------------------------------------- #


logit.spike.edit <- function (formula, niter, data, subset, prior = NULL, na.action = options("na.action"), 
          contrasts = NULL, drop.unused.levels = TRUE, initial.value = NULL, 
          ping = niter/10, nthreads = 0, clt.threshold = 2, mh.chunk.size = 10, 
          proposal.df = 3, sampler.weights = c(DA = 0.333, RWM = 0.333, 
                                               TIM = 0.333), seed = NULL, ...) 
{
  stopifnot(is.numeric(sampler.weights), length(sampler.weights) == 
              3, "DA" %in% names(sampler.weights), "TIM" %in% names(sampler.weights), 
            "RWM" %in% names(sampler.weights), all(sampler.weights >= 
                                                     0), all(sampler.weights <= 1), abs(sum(sampler.weights) - 
                                                                                          1) < 0.01)
  has.data <- !missing(data)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), 
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- drop.unused.levels
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  response <- model.response(mf, "any")
  if (!is.null(dim(response)) && length(dim(response)) > 1) {
    stopifnot(length(dim(response)) == 2, ncol(response) == 
                2)
    ny <- response[, 1] + response[, 2]
    response <- response[, 1]
  }
  else {
    response <- response > 0
    ny <- rep(1, length(response))
  }
  design <- model.matrix(mt, mf, contrasts)
  if (is.null(prior)) {
    prior <- LogitZellnerPrior(design, response, ...)
  }
  stopifnot(inherits(prior, "SpikeSlabPriorBase"))
  stopifnot(inherits(prior, "LogitPrior"))
  if (!is.null(initial.value)) {
    if (inherits(initial.value, "logit.spike")) {
      stopifnot(colnames(initial.value$beta) == colnames(design))
      beta0 <- as.numeric(tail(initial.value$beta, 1))
    }
    else if (inherits(initial.value, "glm")) {
      stopifnot(colnames(initial.value$beta) == colnames(design))
      beta0 <- coef(initial.value)
    }
    else if (is.numeric(initial.value)) {
      stopifnot(length(initial.value) == ncol(design))
      beta0 <- initial.value
    }
    else {
      stop("initial.value must be a 'logit.spike' object, a 'glm' object,", 
           "or a numeric vector")
    }
  }
  else {
    beta0 <- prior$mu
  }
  stopifnot(is.matrix(design), nrow(design) == length(response), 
            length(prior$mu) == ncol(design), length(prior$prior.inclusion.probabilities) == 
              ncol(design), all(ny >= response), all(response >= 
                                                       0))
  if (is.null(prior$max.flips)) {
    prior$max.flips <- -1
  }
  if (!is.null(seed)) {
    seed <- as.integer(seed)
  }
  sampler.weights <- sampler.weights[c("DA", "RWM", "TIM")]
  ans <- .Call("logit_spike_slab_wrapper", as.matrix(design), 
               as.integer(response), as.integer(ny), prior, as.integer(niter), 
               as.integer(ping), as.integer(nthreads), beta0, as.integer(clt.threshold), 
               as.integer(mh.chunk.size), sampler.weights, seed)
  ans$prior <- prior
  class(ans) <- c("logit.spike", "glm.spike")
  ans$contrasts <- attr(design, "contrasts")
  ans$xlevels <- .getXlevels(mt, mf)
  ans$call <- cl
  ans$terms <- mt
  #### Editing out computation of fitted logits because they cause a huge memory spike
  #### (around three times the memory needed for the rest of the process)
  # fitted.logits <- design %*% t(ans$beta)
  # log.likelihood.contributions <- response * fitted.logits + 
  #   ny * plogis(fitted.logits, log.p = TRUE, lower.tail = FALSE)
  # ans$log.likelihood <- colSums(log.likelihood.contributions)
  # sign <- rep(1, length(response))
  # sign[response/ny < 0.5] <- -1
  # ans$deviance.residuals <- sign * sqrt(rowMeans(-2 * log.likelihood.contributions))
  # p.hat <- sum(response)/sum(ny)
  # ans$null.log.likelihood <- sum(response * log(p.hat) + (ny - 
  #                                                           response) * log(1 - p.hat))
  # fitted.probabilities <- plogis(fitted.logits)
  # ans$fitted.probabilities <- rowMeans(fitted.probabilities)
  # ans$fitted.logits <- rowMeans(fitted.logits)
  if (!is.null(initial.value) && inherits(initial.value, "logit.spike")) {
    ans$beta <- rbind(initial.value$beta, ans$beta)
  }
  ans$response <- response
  if (any(ny != 1)) {
    ans$trials <- ny
  }
  colnames(ans$beta) <- colnames(design)
  
  #### Editing out the below as it causes code to return design matrix 
  #### This is a problem for me because design matrix is huge
  # if (has.data) {
  #   ans$training.data <- data
  # }
  # else {
  #   ans$training.data <- mf
  # }
  class(ans) <- c("logit.spike", "lm.spike", "glm.spike")
  return(ans)
}
