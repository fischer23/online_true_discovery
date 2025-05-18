# This file contains the Algorithms used for simulations
# of the paper "Admissible online closed testing must employ e-values"

##### General SeqE-Guard algorithm
### Input:
# e_vals:        vector of e-values.
# idx_rejects:   vector containing indices of rejections. Maximum index should be smaller than length of e_vals.

### Output:      Simultaneous bound for the number of true discoveries in the rejection set defined by idx_rejects.

SeqE_Guard <- function(e_vals, idx_rejects) {
  M <- length(e_vals)
  A <- c()
  U <- c()
  bound <- rep(0, M)

  for (i in 1:M) {
    if (i %in% idx_rejects) {
      A <- c(A, i)
      if (prod(e_vals[c(A, U)]) >= 1 / alpha) {
        bound[i] <- ifelse(i > 1, bound[i - 1] + 1, 1)
        A <- A[-max(which(e_vals[A] == max(e_vals[A])))]
      } else {
        bound[i] <- ifelse(i > 1, bound[i - 1], 0)
      }
    } else if (e_vals[i] < 1) {
      U <- c(U, i)
      bound[i] <- ifelse(i > 1, bound[i - 1], 0)
    } else {
      bound[i] <- ifelse(i > 1, bound[i - 1], 0)
    }
  }
  return(bound)
}

##### Boosted SeqE-Guard algorithm
### Input:
# e_vals:        vector of e-values.
# idx_rejects:   vector containing indices of rejections. Maximum index should be smaller than length of e_vals.
# delta:         Parameter used for boosting the e-values. In case of the simple null H_i:X_i~N(mu_N,1)
#               vs. alternative H_i^A:X_i~N(mu_A,1) delta should be set to delta=mu_A-mu_N.

### Output:      (Boosted) Simultaneous bound for the number of true discoveries in the rejection set defined by idx_rejects.

boosted_SeqE_Guard <- function(e_vals, idx_rejects, delta) {
  M <- length(e_vals)
  A <- c()
  U <- c()
  bound <- rep(0, M)
  b <- rep(1, M)

  for (i in 1:M) {
    m_t <- max(max(c(e_vals[A], 0)), 1 / (alpha * prod(e_vals[c(A, U)])))
    if (m_t == Inf) {
      b[i] <- 1
    } else {
      b[i] <- max(1, e_boosted(m_t, delta))
    }
    e_vals[i] <- e_vals[i] * b[i]
    if (i %in% idx_rejects) {
      A <- c(A, i)
      if (prod(e_vals[c(A, U)]) >= 1 / alpha) {
        bound[i] <- ifelse(i > 1, bound[i - 1] + 1, 1)
        A <- A[-max(which(e_vals[A] == max(e_vals[A])))]
      } else {
        bound[i] <- ifelse(i > 1, bound[i - 1], 0)
      }
    } else if (e_vals[i] < 1) {
      U <- c(U, i)
      bound[i] <- ifelse(i > 1, bound[i - 1], 0)
    } else {
      bound[i] <- ifelse(i > 1, bound[i - 1], 0)
    }
  }
  return(bound)
}

##### hedged and boosted SeqE-Guard algorithm
### Input:
# e_vals:        vector of e-values.
# idx_rejects:   vector containing indices of rejections. Maximum index should be smaller than length of e_vals.
# delta:         Parameter used for boosting the e-values. In case of the simple null H_i:X_i~N(mu_N,1)
#               vs. alternative H_i^A:X_i~N(mu_A,1) delta should be set to delta=mu_A-mu_N.

### Output:      Simultaneous bound for the number of true discoveries in the rejection set defined by idx_rejects.
#               Before boosting the (GRO) e-values they are hedged according to the previous number of e-values larger than 1.


hedged_boosted_SeqE_Guard <- function(e_vals, idx_rejects, delta) {
  M <- length(e_vals)
  A <- c()
  U <- c()
  bound <- rep(0, M)
  b <- rep(1, M)

  taus <- c(1 / 2, (cumsum((e_vals[1:(M - 1)] > 1)) + 1 / 2) / (2:M))
  e_vals <- taus * e_vals + (1 - taus)

  for (i in 1:M) {
    m_t <- max(max(c(e_vals[A], 0)), 1 / (alpha * prod(e_vals[c(A, U)])))
    tau <- taus[i]
    if (m_t == Inf) {
      b[i] <- 1
    } else {
      b[i] <- max(1, e_boosted_hedged(m_t, delta, tau))
    }
    e_vals[i] <- e_vals[i] * b[i]
    if (i %in% idx_rejects) {
      A <- c(A, i)
      if (prod(e_vals[c(A, U)]) >= 1 / alpha) {
        bound[i] <- ifelse(i > 1, bound[i - 1] + 1, 1)
        A <- A[-max(which(e_vals[A] == max(e_vals[A])))]
      } else {
        bound[i] <- ifelse(i > 1, bound[i - 1], 0)
      }
    } else if (e_vals[i] < 1) {
      U <- c(U, i)
      bound[i] <- ifelse(i > 1, bound[i - 1], 0)
    } else {
      bound[i] <- ifelse(i > 1, bound[i - 1], 0)
    }
  }
  return(bound)
}



##### Mean SeqE-Guard algorithm
### Input:
# e_vals:        matrix of e-values. Each row contains the e-values for different hypotheses for one parameter.
# weights:       vector of weights with sum less or equal than 1. Length should be same as number of rows of e_vals.
# idx_rejects:   vector containing indices of rejections. Maximum index should be smaller than length of e_vals.

### Output:      Simultaneous bound for the number of true discoveries in the rejection set defined by idx_rejects.

SeqE_Guard_average <- function(e_vals, weights, idx_rejects) {
  M <- length(e_vals[1, ])
  A <- c()
  U <- c()
  bound <- rep(0, M)

  for (i in 1:M) {
    if (i %in% idx_rejects) {
      A <- c(A, i)
      if (length(c(A, U)) > 1) {
        if (sum(weights * rowProds(e_vals[, c(A, U)])) >= 1 / alpha) {
          bound[i] <- ifelse(i > 1, bound[i - 1] + 1, 1)
          A <- A[-max(which(e_vals[1, A] == max(e_vals[1, A])))]
        } else {
          bound[i] <- ifelse(i > 1, bound[i - 1], 0)
        }
      } else {
        if (sum(weights * e_vals[, c(A, U)]) >= 1 / alpha) {
          bound[i] <- ifelse(i > 1, bound[i - 1] + 1, 1)
          A <- A[-max(which(e_vals[1, A] == max(e_vals[1, A])))]
        } else {
          bound[i] <- ifelse(i > 1, bound[i - 1], 0)
        }
      }
    } else if (e_vals[1, i] < 1) {
      U <- c(U, i)
      bound[i] <- ifelse(i > 1, bound[i - 1], 0)
    } else {
      bound[i] <- ifelse(i > 1, bound[i - 1], 0)
    }
  }
  return(bound)
}
