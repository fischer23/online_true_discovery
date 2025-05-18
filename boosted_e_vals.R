# This file contains the functions for boosting the e-values as described in Section 3.5
# of the paper "Admissible online closed testing must employ e-values"

### General boosting method
## Input:
# m_t:           cutoff value for the e-value. For concrete calculation see equation (11)
# delta:         Parameter used for boosting the e-values. In case of the simple null H_i:X_i~N(mu_N,1)
#               vs. alternative H_i^A:X_i~N(mu_A,1) delta should be set to delta=mu_A-mu_N.

### Output:      boosting factor.

e_boosted <- function(m_t, delta) {
  b_factor <- function(b) {
    return(b * (1 - pnorm(delta / 2 - log(m_t / b) / delta)) + m_t * (1 - pnorm((log(m_t / b) + delta^2 / 2) / delta)) - 1)
  }
  return(uniroot(b_factor, lower = 0, upper = 100000000000)$root)
}

### General boosting method
## Input:
# m_t:           cutoff value for the e-value. For concrete calculation see equation (11)
# delta:         Parameter used for boosting the e-values. In case of the simple null H_i:X_i~N(mu_N,1)
#               vs. alternative H_i^A:X_i~N(mu_A,1) delta should be set to delta=mu_A-mu_N.
# tau:           factor that was previously used to hedge the e-value

### Output:      boosting factor for hedged e-value.

e_boosted_hedged <- function(m_t, delta, tau) {
  b_factor <- function(b) {
    s <- (tau - 1 + m_t / b) / tau
    return(b * tau * (1 - pnorm(delta / 2 - log(s) / delta)) + m_t * (1 - pnorm((log(s) + delta^2 / 2) / delta)) + b * (1 - tau) * pnorm((log(s) + delta^2 / 2) / delta) - 1)
  }
  return(uniroot(b_factor, lower = 0, upper = 1 / (1 - tau))$root)
}
