#########################################################################################################
# Authors: Mark Zwart                                                                                   #
# Date: April 6, 2020 	                                                                                #
#                                                                                                       #
# Estimating cell density from colony counts from different dilutions. This is probably not the most    #
# elegant way to get the solution, but it works.							#
#                                                                                                       #
#########################################################################################################

# Example data
data.dil    = c( 10^-4, 10^-4, 10^-6, 10^-6) # Dilutions used 
data.counts = c(   320,   233,     5,     0) # Colony counts corresponding to dilutions

# Calculate the expected cell concentrations based on the plates with colonies only, to get a rough idea of
# the solution.
non.zero = which(data.counts > 0)
exp.cells = mean(data.counts[non.zero]/data.dil[non.zero])
log.exp.cells = log10(exp.cells)

# Determine the search range.
range.search = 1 	# Limit of parameter space for search for cell count, log scale
units = 0.001 		# Units to use in the parameter space, log scale
log.par.space = seq((log.exp.cells - range.search), (log.exp.cells + range.search), by = units)
par.space = 10^log.par.space
n.par.space = length(par.space)

# Array to store the results
results.array = array(NA, dim = c(n.par.space,2))
colnames(results.array) = c("Cells", "NLL")

# Loop to test all values in the parameter space, and for each value determine the negative log likelihood
# of the observations, assuming a Poisson error structure.

for(i in 1:n.par.space) {
	nll.now = sum(-dpois(x = (data.counts), lambda = par.space[i]*data.dil, log = TRUE))
	results.array[i,] = c(par.space[i], nll.now)	
} # end of i loop
	
# Order results based on NLL values and print the best result
order.results = order(results.array[,2])
estimate = results.array[order.results[1],1]
print(estimate) # This is the estimate of cell count.

# Plot cell count vs. NLL as a check, the Poisson-based solution and the mean solution
plot(x = log.par.space, y = results.array[,2], type = "l", xlab = "Log10 Cells", ylab = "NLL", col = "blue")
points(x = log.exp.cells, y = 0, pch = 2, col = "magenta") # Solution based on means for non-zero plates
points(x = log10(estimate), y = 0, pch = 16, col = "red") # Poisson solution

