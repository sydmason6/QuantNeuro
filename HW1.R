##Exercise 1
#It would most likely not be statistically significant at p<0.05. However,whether or not this is significant would depend on the overall prevalance of HIV (i.e. the prior probability)


##Exercise 2
# Define prior probability (P(A))
prior_probability <- 0.1

# Define likelihood (P(B|A))
likelihood <- 0.95

# Define evidence (P(B))
evidence <- 1

# Calculate the posterior probability (P(A|B))
posterior_probability <- (likelihood * prior_probability) / evidence

# Print the result
cat("The posterior probability P(A|B) is:", posterior_probability, "\n")

#Range of prior probabilities
prior_range <- seq(0, 1, by=0.1)

# Calculate the posterior probability w/range of priors(P(A|B))
posterior_probability <- (likelihood * prior_range) / evidence

# Print the result
cat("The posterior probability P(A|B) is:", posterior_probability, "\n")
