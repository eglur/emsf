source("blackjack.R")
source("data.plot.R")
source("dp.R")
source("emsf.R")
source("util.R")


# counting (maximum likelihood)
# ne: # de jogos para avaliar a política
bj.ml <- function(num.episodes = 3e4, epsilon = 0.15, tc = 1000, df = 0.999, max.iter = 200, ne = 1e6)
{

  P <- array(0, c(203, 203, 2))
  r <- matrix(0, 203, 2)
  C <- array(0, c(203, 203, 2)) # counting
  rs <- matrix(0, 203, 2) # sums

  # structure for the specific problem
  P[201:203,,] <- 0
  P[201,201,] <- 1
  P[202,202,] <- 1
  P[203,203,] <- 1

  #   for (i in c(201, 202, 203)) P[i,i,1:2] <- 1
  pi <- sample(1:2, 203, TRUE)

  R <- array(0, num.episodes %/% tc)
  T <- array(0, num.episodes %/% tc)

  start.time <- proc.time()[1]
  time.eval <- 0
  for (ep in 1:num.episodes)
  {

    ph <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))
    dh <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))

    while (ph$sum.cards < 12) ph <- bj.draw.card(ph) # nothing to do otherwise

    s <- list(ph = ph, dh = dh, idx = 0)
    j <- bj.state.index(s)

    t <- NULL
    while (s$idx == 0)
    {

      i <- j
      if (runif(1) < epsilon) a <- sample(1:2, 1)
      else a <- pi[i]

      t <- bj.transition(s, a)

      s <- t$s
      j <- bj.state.index(s)

      C[i,j,a] <- C[i,j,a] + 1
      rs[i,a] <- rs[i,a] + t$r

    }

   if (ep %% tc == 0)
   {
      for (u in 1:2)
      {
         cnt <- apply(C[,,u], 1, sum)
         for (k in 1:length(cnt))
         {
            if (cnt[k] > 0)
            {
               P[k,,u] <- C[k,,u] / cnt[k]
               r[k,u] <- rs[k,u] / cnt[k]
            }
            else P[k,,u] <- 1 / ncol(P[,,u])
         }
      }

      pi <- policy.iteration(r, P, df = df, max.iter = max.iter)$pi

      tt <- proc.time()[1]
      R[ep / tc] <- bj.evaluate.policy.fast(pi, ne)

#       print(R[ep / tc])
#
      time.eval <- time.eval + (proc.time()[1] - tt)

      T[ep / tc] <- (proc.time()[1] - start.time) - time.eval
   }


  }

  list(R = R, T = T)
}



bj.emsf <- function(m, alpha = 1e-1, beta = NULL,
                    num.episodes = 5e3, epsilon = 0.15, tc = 10, df = 0.999, max.iter = 200, ne = 1e6)
{
  if (is.null(beta)) beta <- alpha * 0.5

  n <- 203

  D  <- array(runif(n * m * 2), c(n, m, 2))
#   D  <- array(1, c(n, m, 2))
  for (u in 1:2) D[,,u] <- D[,,u] / apply(D[,,u], 1, sum)

  K  <- matrix(runif((m-3) * n), (m-3), n)
#   K  <- matrix(1, (m-3), n)
  K <- K / apply(K, 1, sum)

  K <- rbind(K, matrix(0, 3, n))
  K[m-2, 201] <- 1
  K[m-1, 202] <- 1
  K[m, 203] <- 1

  D[201:203,,1:2] <- 0
  D[201, m-2,1:2] <- 1
  D[202, m-1,1:2] <- 1
  D[203, m, 1:2] <- 1

  rp <- array(0, n)
  rc <- array(0, n) # counting

  pi <- sample(1:2, n, TRUE)

  R <- array(0, num.episodes %/% tc)
  T <- array(0, num.episodes %/% tc)

  start.time <- proc.time()[1]
  time.eval <- 0
  for (ep in 1:num.episodes)
  {

    ph <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))
    dh <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))

    while (ph$sum.cards < 12) ph <- bj.draw.card(ph) # nothing to do otherwise

    s <- list(ph = ph, dh = dh, idx = 0)
    j <- bj.state.index(s)

    t <- NULL
    while (s$idx == 0)
    {

      i <- j
      if (runif(1) < epsilon) a <- sample(1:2, 1)
      else a <- pi[i]

      t <- bj.transition(s, a)

      s <- t$s
      j <- bj.state.index(s)

      w <- D[i,,a] * K[,j]
      w <- w / sum(w)

      D[i,,a] <- (1 - alpha) * D[i,,a] + alpha * w
      K[,j] <- (1 - beta) * K[,j] + beta * w

      rp[j] <- rp[j] + t$r
      rc[j] <- rc[j] + 1

    }

   if (ep %% tc == 0)
   {

      K <- K / apply(K, 1, sum)

      rt <- array(0, n)
      for (i in 1:length(rc)) if (rc[i] > 0) rt[i] <- rp[i] / rc[i]
      rb <- K %*% rt

      pi <- pisf(D, K, rb, df = df, max.iter = max.iter)$pi

      tt <- proc.time()[1]
      R[ep / tc] <- bj.evaluate.policy.fast(pi, ne)
      time.eval <- time.eval + (proc.time()[1] - tt)

      T[ep / tc] <- (proc.time()[1] - start.time) - time.eval
   }

  }

  list(R = R, T = T)
}


bj.emsf.full <- function(m = 10, alpha = 1, tcc = 100,
                    num.episodes = 3e4, epsilon = 0.15, tc = 1000, df = 0.999, max.iter = 200, ne = 1e6)
{

  if (is.null(tcc)) tcc <- tc

  n <- 203

  D  <- array(runif(n * m * 2), c(n, m, 2))
#   D  <- array(1, c(n, m, 2))
  for (u in 1:2) D[,,u] <- D[,,u] / apply(D[,,u], 1, sum)

  K  <- matrix(runif((m-3) * n), (m-3), n)
#   K  <- matrix(1, (m-3), n)
  K <- K / apply(K, 1, sum)

  K <- rbind(K, matrix(0, 3, n))
  K[m-2, 201] <- 1
  K[m-1, 202] <- 1
  K[m, 203] <- 1

  D[201:203,,1:2] <- 0
  D[201, m-2,1:2] <- 1
  D[202, m-1,1:2] <- 1
  D[203, m, 1:2] <- 1


  Dh <- D
  Kh <- K
  Dh[1:200, ,] <- 0
  Kh[1:(m-3),] <- 0

  rp <- array(0, n)
  rc <- array(0, n) # counting
  rb <- array(0, m)

  pi <- sample(1:2, n, TRUE)

  R <- array(0, num.episodes %/% tc)
  T <- array(0, num.episodes %/% tc)

  start.time <- proc.time()[1]
  time.eval <- 0
  for (ep in 1:num.episodes)
  {

    ph <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))
    dh <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))

    while (ph$sum.cards < 12) ph <- bj.draw.card(ph) # nothing to do otherwise

    s <- list(ph = ph, dh = dh, idx = 0)
    j <- bj.state.index(s)

    t <- NULL
    while (s$idx == 0)
    {

      i <- j
      if (runif(1) < epsilon) a <- sample(1:2, 1)
      else a <- pi[i]

      t <- bj.transition(s, a)

      s <- t$s
      j <- bj.state.index(s)

      w <- D[i,,a] * K[,j]
      g <- sum(w)
      if (g > 0)
      {
        w <- w / g
        Dh[i,,a] <- Dh[i,,a] + w
        Kh[,j] <- Kh[,j] + w
      }

      rp[j] <- rp[j] + t$r
      rc[j] <- rc[j] + 1

    }

   if (ep %% tcc == 0)
    {

      for (u in 1:2)
      {
        for (i in 1:nrow(Dh[,,u]))
        {
          ss <- sum(Dh[i,,u])
          if (ss > 0) D[i,,u] <- (1 - alpha) *  D[i,,u] + alpha * Dh[i,,u] / ss
        }
      }

      K <- (1 - alpha) * K + alpha * Kh / apply(Kh, 1, sum)


      rt <- array(0, n)
      for (i in 1:length(rc)) if (rc[i] > 0) rt[i] <- rp[i] / rc[i]
      rb <- K %*% rt

      Dh <- D
      Kh <- K
      Dh[1:200, ,] <- 0
      Kh[1:(m-3),] <- 0

    }

    if (ep %% tc == 0)
    {

      pi <- pisf(D, K, rb, df = df, max.iter = max.iter)$pi

      tt <- proc.time()[1]
      R[ep / tc] <- bj.evaluate.policy.fast(pi, ne)
#       print(R[ep / tc])
      time.eval <- time.eval + (proc.time()[1] - tt)

      T[ep / tc] <- (proc.time()[1] - start.time) - time.eval
    }

  }

  list(R = R, T = T)
}



# usar esse
# df: gamma do PISF
# tcc: # transicoes para plotar
bj.emsf.comp <- function(m = 10, alpha = 1, tcc = 100,
                    num.episodes = 3e4, epsilon = 0.15, tc = 1000, df = 0.999, max.iter = 200, ne = 1e6)
{

  C <- array(0, c(203, 203, 2)) # counting

  if (is.null(tcc)) tcc <- tc

  n <- 203
  na <- 2

  D  <- array(runif(n * m * 2), c(n, m, 2))
#   D  <- array(1, c(n, m, 2))
  for (u in 1:2) D[,,u] <- D[,,u] / apply(D[,,u], 1, sum)

  K  <- matrix(runif((m-3) * n), (m-3), n)
#   K  <- matrix(1, (m-3), n)
  K <- K / apply(K, 1, sum)

  K <- rbind(K, matrix(0, 3, n))
  K[m-2, 201] <- 1
  K[m-1, 202] <- 1
  K[m, 203] <- 1

  D[201:203,,1:2] <- 0
  D[201, m-2,1:2] <- 1
  D[202, m-1,1:2] <- 1
  D[203, m, 1:2] <- 1


  Dh <- D
  Kh <- K
  Dh[1:200, ,] <- 0
  Kh[1:(m-3),] <- 0

  rp <- array(0, n)
  rc <- array(0, n) # counting
  rb <- array(0, m)

  pi <- sample(1:2, n, TRUE)

  R <- array(0, num.episodes %/% tc)
  T <- array(0, num.episodes %/% tc)

  start.time <- proc.time()[1]
  time.eval <- 0
  for (ep in 1:num.episodes)
  {

    ph <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))
    dh <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))

    while (ph$sum.cards < 12) ph <- bj.draw.card(ph) # nothing to do otherwise

    s <- list(ph = ph, dh = dh, idx = 0)
    j <- bj.state.index(s)

    t <- NULL
    while (s$idx == 0)
    {

      i <- j
      if (runif(1) < epsilon) a <- sample(1:2, 1)
      else a <- pi[i]

      t <- bj.transition(s, a)

      s <- t$s
      j <- bj.state.index(s)

      C[i,j,a] <- C[i,j,a] + 1

      rp[j] <- rp[j] + t$r
      rc[j] <- rc[j] + 1

    }

   if (ep %% tcc == 0)
    {

      for (a in 1:na)
      {
        for (i in 1:nrow(C[,,a]))
        {
          for (j in 1:ncol(C[,,a]))
          {
            if (C[i,j,a] != 0)
            {
              w <- D[i,,a] * K[,j]
              g <- sum(w)
              if (g > 0)
              {
                w <- C[i,j,a] * w / g
                Dh[i,,a] <- Dh[i,,a] + w
                Kh[,j] <- Kh[,j] + w
              }
            }
          }
        }
      }

      for (u in 1:2)
      {
        for (i in 1:nrow(Dh[,,u]))
        {
          ss <- sum(Dh[i,,u])
          if (ss > 0) D[i,,u] <- (1 - alpha) *  D[i,,u] + alpha * Dh[i,,u] / ss
        }
      }

      K <- (1 - alpha) * K + alpha * Kh / apply(Kh, 1, sum)


      rt <- array(0, n)
      for (i in 1:length(rc)) if (rc[i] > 0) rt[i] <- rp[i] / rc[i]
      rb <- K %*% rt

      Dh <- D
      Kh <- K
      Dh[1:200, ,] <- 0
      Kh[1:(m-3),] <- 0

    }

    if (ep %% tc == 0)
    {

      pi <- pisf(D, K, rb, df = df, max.iter = max.iter)$pi

      tt <- proc.time()[1]
      R[ep / tc] <- bj.evaluate.policy.fast(pi, ne)
#       print(R[ep / tc])
      time.eval <- time.eval + (proc.time()[1] - tt)

      T[ep / tc] <- (proc.time()[1] - start.time) - time.eval
    }

  }

  list(R = R, T = T)
}




bj.qlearning <- function(alpha = 1e-1, num.episodes = 5e3, epsilon = 0.15, tc = 10, df = 0.999, ne = 1e6)
{
  n <- 203
  Q <- matrix(runif(n*2, 0, 1e-5), n, 2) # optimistic initialization
  Q[201:203,] <- 0

  pi <- apply(Q, 1, which.max)

  R <- array(0, num.episodes %/% tc)
  T <- array(0, num.episodes %/% tc)

  start.time <- proc.time()[1]
  time.eval <- 0
  for (ep in 1:num.episodes)
  {

    ph <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))
    dh <- bj.draw.card(list(sum.cards = 0, has.ace = FALSE))

    while (ph$sum.cards < 12) ph <- bj.draw.card(ph) # nothing to do otherwise

    s <- list(ph = ph, dh = dh, idx = 0)
    j <- bj.state.index(s)

    t <- NULL
    while (s$idx == 0)
    {

      i <- j
      if (runif(1) < epsilon) a <- sample(1:2, 1)
      else a <- pi[i]

      t <- bj.transition(s, a)

      s <- t$s
      j <- bj.state.index(s)

      Q[i,a] <- (1 - alpha) * Q[i,a] + alpha * (t$r + df * max(Q[j,]))

    }

   if (ep %% tc == 0)
   {
      pi <- apply(Q, 1, which.max)

      tt <- proc.time()[1]
      R[ep / tc] <- bj.evaluate.policy.fast(pi, ne)
      time.eval <- time.eval + (proc.time()[1] - tt)

      T[ep / tc] <- (proc.time()[1] - start.time) - time.eval
   }


  }

  list(R = R, T = T)
}


# num.avg: número de execuções (rodadas, repetições)
bj.emsf.comp.experiment <- function(m, alpha, tcc = 100, num.episodes = 3e4, epsilon = 0.15, tc = 500, num.avg = 50, load.previous = TRUE, dir = "./files/", ne = 1e6)
{

  if (is.null(tcc)) tcc <- tc

  prefix <- paste("bj_emsf_comp", num.episodes, epsilon, tc, m, alpha, tcc, sep= "_")
  file.return <- ps(dir, prefix, "_ret.txt")
  file.time   <- ps(dir, prefix, "_tim.txt")

  R <- matrix(0, num.episodes %/% tc, num.avg)
  T <- matrix(0, num.episodes %/% tc, num.avg)

  run <- 1
  # checks if a previous experiment was interrupted
  if (load.previous && file.exists(file.return) && file.exists(file.time))
  {
    R <- read.table(file.return)
    sums <- apply(R, 2, sum)
    while (run <= length(sums) && sums[run] != 0) run <- run + 1

    run2 <- 1
    T <- read.table(file.time)
    sums <- apply(T, 2, sum)
    while (run2 <= length(sums) && sums[run2] > 0) run2 <- run2 + 1

    run <- min(run,run2)
    if (run <= num.avg) print(paste("WARNING: starting experiments with bj_emsf at run", run))
    else print("WARNING: experiments with bj_emsf already finished")
  }

  if (run <= num.avg)
  {
    for (i in run:num.avg)
    {
      Y <- bj.emsf.comp(m, alpha, tcc, num.episodes, epsilon, tc, ne = ne)
      R[,i] <- Y$R
      T[,i] <- Y$T

      print(paste("Run", i, "  Return", mean(Y$R), "  Time", Y$T[length(Y$T)]))

      wt(R, file.return)
      wt(T, file.time)
    }
  }
}



# ml: maximum likelihood
bj.ml.experiment <- function(num.episodes = 3e4, epsilon = 0.15, tc = 500, num.avg = 50, load.previous = TRUE, dir = "./files/", ne = 1e6)
{

  prefix <- paste("bj_ml", num.episodes, epsilon, tc, sep= "_")
  file.return <- ps(dir, prefix, "_ret.txt")
  file.time   <- ps(dir, prefix, "_tim.txt")

  R <- matrix(0, num.episodes %/% tc, num.avg)
  T <- matrix(0, num.episodes %/% tc, num.avg)

  run <- 1
  # checks if a previous experiment was interrupted
  if (load.previous && file.exists(file.return) && file.exists(file.time))
  {
    R <- read.table(file.return)
    sums <- apply(R, 2, sum)
    while (run <= length(sums) && sums[run] != 0) run <- run + 1

    run2 <- 1
    T <- read.table(file.time)
    sums <- apply(T, 2, sum)
    while (run2 <= length(sums) && sums[run2] > 0) run2 <- run2 + 1

    run <- min(run,run2)
    if (run <= num.avg) print(paste("WARNING: starting experiments with bj_ml at run", run))
    else print("WARNING: experiments with bj_ml already finished")
  }

  if (run <= num.avg)
  {
    for (i in run:num.avg)
    {
      Y <- bj.ml(num.episodes, epsilon, tc, ne = ne)
      R[,i] <- Y$R
      T[,i] <- Y$T

      print(paste("Run", i, "  Return", mean(Y$R), "  Time", Y$T[length(Y$T)]))

      wt(R, file.return)
      wt(T, file.time)
    }
  }
}



bj.emsf.experiment <- function(m, alpha, beta, num.episodes, epsilon, tc, num.avg, load.previous = TRUE, dir = "./files/", ne = 1e6)
{

  prefix <- paste("bj_emsf", num.episodes, epsilon, tc, m, alpha, beta, sep= "_")
  file.return <- ps(dir, prefix, "_ret.txt")
  file.time   <- ps(dir, prefix, "_tim.txt")

  R <- matrix(0, num.episodes %/% tc, num.avg)
  T <- matrix(0, num.episodes %/% tc, num.avg)

  run <- 1
  # checks if a previous experiment was interrupted
  if (load.previous && file.exists(file.return) && file.exists(file.time))
  {
    R <- read.table(file.return)
    sums <- apply(R, 2, sum)
    while (run <= length(sums) && sums[run] != 0) run <- run + 1

    run2 <- 1
    T <- read.table(file.time)
    sums <- apply(T, 2, sum)
    while (run2 <= length(sums) && sums[run2] > 0) run2 <- run2 + 1

    run <- min(run,run2)
    if (run <= num.avg) print(paste("WARNING: starting experiments with bj_emsf at run", run))
    else print("WARNING: experiments with bj_emsf already finished")
  }

  if (run <= num.avg)
  {
    for (i in run:num.avg)
    {
      Y <- bj.emsf(m, alpha, beta, num.episodes, epsilon, tc, ne = ne)
      R[,i] <- Y$R
      T[,i] <- Y$T

      print(paste("Run", i, "  Return", mean(Y$R), "  Time", Y$T[length(Y$T)]))

      wt(R, file.return)
      wt(T, file.time)
    }
  }
}


# vamos usar esse
bj.qlearning.experiment <- function(alpha, num.episodes, epsilon, tc, num.avg, load.previous = TRUE, dir = "./files/", ne = 1e6)
{

  prefix <- paste("bj_qlearning", num.episodes, epsilon, tc, alpha, sep= "_")
  file.return <- ps(dir, prefix, "_ret.txt")
  file.time   <- ps(dir, prefix, "_tim.txt")

  R <- matrix(0, num.episodes %/% tc, num.avg)
  T <- matrix(0, num.episodes %/% tc, num.avg)

  run <- 1
  # checks if a previous experiment was interrupted
  if (load.previous && file.exists(file.return) && file.exists(file.time))
  {
    R <- read.table(file.return)
    sums <- apply(R, 2, sum)
    while (run <= length(sums) && sums[run] != 0) run <- run + 1

    run2 <- 1
    T <- read.table(file.time)
    sums <- apply(T, 2, sum)
    while (run2 <= length(sums) && sums[run2] > 0) run2 <- run2 + 1

    run <- min(run,run2)
    if (run <= num.avg) print(paste("WARNING: starting experiments with bj_qlearning at run", run))
    else print("WARNING: experiments with bj_qlearning already finished")
  }

  if (run <= num.avg)
  {
    for (i in run:num.avg)
    {
      Y <- bj.qlearning(alpha, num.episodes, epsilon, tc, ne = ne)
      R[,i] <- Y$R
      T[,i] <- Y$T

      print(paste("Run", i, "  Return", mean(Y$R), "  Time", Y$T[length(Y$T)]))

      wt(R, file.return)
      wt(T, file.time)
    }
  }
}

std <- function(x)
{
    apply(x, 1, sd)/sqrt(nrow(x))
}

plot.results <- function(num.episodes = 3e4, epsilon = 0.15, tc = 500, emsf.ms = 10, emsf.alpha = 1, emsf.beta, qlearning.alpha, dir = "./files/")
{
    R <- NULL
    S <- NULL
    T <- NULL
    legendas <- NULL;

    prefix <- paste("bj_ml", num.episodes, epsilon, tc, sep= "_")
    file.return <- ps(dir, prefix, "_ret.txt")
    file.time   <- ps(dir, prefix, "_tim.txt")

    result <- read.table(file.return)
    result.std <- std(result)
    result.mean <- apply(result, 1, mean)

    R <- cbind(R, result.mean)
    S <- cbind(S, result.std)
    T <- cbind(T, apply(read.table(file.time), 1, mean))

    legendas <- c(legendas, paste("ML"))

    prefix <- paste("bj_qlearning", num.episodes, epsilon, tc, qlearning.alpha, sep= "_")
    file.return <- ps(dir, prefix, "_ret.txt")
    file.time   <- ps(dir, prefix, "_tim.txt")

    result <- read.table(file.return)
    result.std <- std(result)
    result.mean <- apply(result, 1, mean)

    R <- cbind(R, result.mean)
    S <- cbind(S, result.std)
    T <- cbind(T, apply(read.table(file.time), 1, mean))

    legendas <- c(legendas, paste("QL"))

    for (i in 1:length(emsf.ms))
    {
        prefix <- paste("bj_emsf_comp", num.episodes, epsilon, tc, emsf.ms[i], emsf.alpha, tc, sep= "_")
        file.return <- ps(dir, prefix, "_ret.txt")
        file.time   <- ps(dir, prefix, "_tim.txt")

        result <- read.table(file.return)
        result.std <- std(result)
        result.mean <- apply(result, 1, mean)

        R <- cbind(R, result.mean)
        S <- cbind(S, result.std)
        T <- cbind(T, apply(read.table(file.time), 1, mean))

        legendas <- c(legendas, paste("QL"))
    }

    R <- R[1:10,]
    S <- S[1:10,]

    mp(seq(tc, num.episodes, length=(nrow(R))), R, R+S, R-S, t="l")
    leg(pos="bottomright", leg=legendas)
}
plot.results(num.episodes=5e3, epsilon=0.15, tc=100, qlearning.alpha=0.1, emsf.ms=c(10, 50, 100), emsf.alpha=1)

