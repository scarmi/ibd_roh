
# Implement this function according to the demography
# N is in number of mating pairs (couples)
# t=1 is present
# Larger values of t are more distant times in the past
get_N = function(t)
{
  # An example
  # The population size was 500 until 100 generations ago
  # Then it was 100 until the beginning of time
  if (t<=100)
    return(500)
  else
    return(100)
}

# "rates" is a vector of probabilities of the relationship between a husband and a wife (a mating pair).
# rates[1] is the probability they are full siblings
# rates[2] is the probability they are full first cousins
# etc., where the vector can be made arbitraily long.
# Note that we only model full relatives (i.e., no half-sibs)
# Also note that these events are mutually exclusive
# So, with probability 1-sum(rates), the spouses are unrelated in the 
# past length(rates) generations.
# The kinship coefficient implied by these rates should be equal
# to the kinship coefficient defined above for the time period 1-max_T_exact_model generations.
rates = c(0.1,0.2,0.3)

# For genrations "max_T_exact_model" and beyond, we will use a simpler approximation
# Implement this function according to the consanguinity rates
# This should return the kinship coefficient between spouses
get_q = function(t)
{
  # An example
  # The kinship coef was was 0.1 until 50 generations ago
  # Then it was 0.05 until the beginning of time
  # Note that the maximal kinship coef is 0.25 (when sibs mate each time)
  if (t<=50)
    return(0.1)
  else
    return(0.05)
}

# Demonstrate how to plot the distribution of segment lengths
ells = seq(3,50,by=0.1)

# res = calc_ell_dist(ells,get_N,rates,get_q)
ell_dist_within = res$dist_within
ell_dist_between = res$dist_between
ell_shared_within = res$shared_within
ell_shared_between = res$shared_between
ell_shared_within_inf = res$shared_within_inf
ell_shared_between_inf= res$shared_between_inf

ells_for_plot = (ells[1:(length(ells)-1)] + ells[2:length(ells)])/2

plot(ells_for_plot,ell_dist_between,log='y',type='l',xlab='Segment length',ylab='Probability',cex.axis=1.5,cex.lab=1.5,lwd=1.5)
lines(ells_for_plot,ell_dist_within,col='red')
legend('topright',c('Between (IBD)','Within (ROH)'),lwd=1.5,bty='n',col=c('black','red'))

plot(ells_for_plot,ell_shared_between,log='y',type='l',xlab='Segment length',ylab='Frac genome shared',cex.axis=1.5,cex.lab=1.5,lwd=1.5)
lines(ells_for_plot,ell_shared_within,col='red')
lines(ells_for_plot,ell_shared_between_inf,lty=2)
lines(ells_for_plot,ell_shared_within_inf,lty=2,col='red')
legend('topright',rep(c('Between (IBD)','Within (ROH)'),2),lwd=1.5,bty='n',col=rep(c('black','red'),2),lty=c(1,1,2,2))
