
# ell_seq: The sequence of segment lengths (ell), for which to return results
# N_func: A function that returns the population size for each t, in number of mating pairs
# rates: A vector indicating the consanguinity rates in the past few generations
# q_func: A function that returns the kinship coefficient between mating pairs beyond the max_T_exact_model past generations. Note that the kinship in the first few generations is now set automatically based on the "rates" argument.
# max_T_exact_model: Specify the number of generations during which we run the exact model with the above rates
# max_T_coal: This is the maximal generation to track when computing the distribution of the TMRCA
# Set a value based on the maximal population size (say, 10x of that size, to guarantee coalescence)
# num_runs: Number of times to run the exact Markov chain (for the distribution of TMRCA in the first few generations)

# bla bla

calc_ell_dist = function(ell_seq,N_func,rates,q_func,max_T_coal=10000,max_T_exact_model=20,num_runs=10000)
{

  n = length(rates)
  
  # This specifies the exact Markov chain
  transition_mat = matrix(0,n+3,n+3)
  # First state is coalescence
  transition_mat[1,] = c(1,rep(0,n+2))
  # Last state (n+3) is two unrelated individuals
  # Because the transition probabilities depend on N, they may change each generation, so we need to   specify them within the simulations below.
  # State 2 is two chromosomes in the same individual
  transition_mat[2,] = c(0,0,1,rep(0,n))
  
  # States 3-(n+2) (n states) are what happens when there are two chrs in the two individuals in a     mating pair
  prefactor = 1
  for (i in 1:n)
  {
    p = rates[i] / 4^(i-1)
    if (i>1)
    {
      prefactor = 1-sum(rates[1:(i-1)] * 4^(1-(1:(i-1))))
    } else {prefactor = 1}
    
    transition_mat[i+2,1:3] = c(1/4,1/4,1/2) * p / prefactor
    transition_mat[i+2,i+3] = 1 - p / prefactor
  }
  
  # This parameter is the number of times we run the exact Markov chain in order to obtain the distribution of the TMRCA. Probably don't need to change it.
  unrelated_state = n+3
  coal_state = 1
  coal_times = numeric()
  rapid_count = 0
  for (j in 1:num_runs)
  {
    state = 2 # Starting with two chrs in the same individual
    time = 1
    while (time<=max_T_exact_model)
    {
      if (state == unrelated_state)
      {
        curN = N_func(time)
        p = c(1/(4*curN),1/(4*curN),1/(2*curN),rep(0,n-1),1-1/curN)
      }
      else {
        p = transition_mat[state,]
      }
      state = sample.int(n+3,1,prob=p)
      time = time + 1
      if (state == coal_state)
      {
       rapid_count = rapid_count + 1
       coal_times = c(coal_times,time)
       break
      }
    }
  }
  
  rapid_mrca_dist = hist(coal_times,1:max(coal_times),plot=F)$counts/num_runs
  
  # This is the formula for the kinship coefficient. No need to change.
  kinship = sum(rates * 4^(-(1:n)))
  
  internal_q = function(t)
  {
    if (t<=max_T_exact_model)
      return(kinship)
    else
      return(q_func(t))
  }
  
  p_coal = numeric(max_T_coal)
  p_no_coal = 1
  for (t in 1:max_T_coal)
  {
    eff_N = N_func(t)*(1-3*internal_q(t))
    p_coal[t] = p_no_coal * 1/(4*eff_N)
    p_no_coal = p_no_coal * (1-1/(4*eff_N))
  }
  
  p_coal_between = p_coal
  p_coal_within = numeric()
  p_coal_within[1:max_T_exact_model] = rapid_mrca_dist
  ptemp = p_coal[(max_T_exact_model+1):max_T_coal]
  p_coal_within[(max_T_exact_model+1):max_T_coal] = ptemp * (1-rapid_count/num_runs) / sum(ptemp)
  
  norm_factor_between = 0
  for (t in 1:max_T_coal)
  {
    norm_factor_between = norm_factor_between + t*p_coal_between[t]
  }
  
  norm_factor_within = 0
  for (t in 1:max_T_coal)
  {
    norm_factor_within = norm_factor_within + t*p_coal_within[t]
  }
  
  # Given what was previously calculated (probability of coalescence in each generation), this function computes the probability of a random segment to have length between ell1 and ell2
  # within should be T for IBD (between individuals), F for ROH (within individuals)
  # Lengths are in centiMorgans
  p_ell = function(within,ell1,ell2)
  {
    ell1 = ell1/100
    ell2 = ell2/100
    if (within)
    {
      pc = p_coal_within
      nm = norm_factor_within
    } else {
      pc = p_coal_between
      nm = norm_factor_between
    }
    p = 0
    for (t in 1:max_T_coal)
    {
      p = p + t*pc[t] * (exp(-2*t*ell1) - exp(-2*t*ell2))
    }
    p = p / nm
    return(p)
  }
  
  chr_lens = c(286.279234,268.839622,223.361095,214.688476,204.089357,192.039918,187.220500,168.003442,  166.359329,181.144008,158.218650,174.679023,125.706316,120.202583,141.860238,134.037726,128.490529,117.708923,107.733846,108.266934,62.786478,74.109562)
  chr_lens = chr_lens/100
  
  f_ell = function(within,ell1,ell2)
  {
    ell1 = ell1/100
    ell2 = ell2/100
    if (within)
      pc = p_coal_within
    else
      pc = p_coal_between
    total_len = 0
    for (c in seq(1,22))
    {
      L = chr_lens[c]
      for (t in 1:max_T_coal)
      {
        total_len_given_t = exp(-2*t*ell1)*(L + 2*L*t*ell1 - 2*t*ell1^2) - 
                            exp(-2*t*ell2)*(L + 2*L*t*ell2 - 2*t*ell2^2)
        total_len = total_len + pc[t] * total_len_given_t
      }
    }
    f = total_len / sum(chr_lens)
    return(f)
  }
  
  f_ell_inf_chr = function(within,ell1,ell2)
  {
    ell1 = ell1/100
    ell2 = ell2/100
    if (within)
      pc = p_coal_within
    else
      pc = p_coal_between
    f = 0
    for (t in 1:max_T_coal)
    {
      f = f + pc[t] * ((1+2*t*ell1)*exp(-2*t*ell1) - (1+2*t*ell2)*exp(-2*t*ell2))
    }
    return(f)
  }
  
  get_dist = function(within, ell_seq)
  {
    cutoff = min(ell_seq)
    ell_dist = numeric(length(ell_seq)-1)
    for (i in 1:(length(ell_seq)-1))
    {
      ell_dist[i] = p_ell(within,ell_seq[i],ell_seq[i+1])
    }
    ell_dist = ell_dist / p_ell(within,cutoff,max(ell_seq))
    return(ell_dist)
  }
  
  get_f = function(within, ell_seq)
  {
    cutoff = min(ell_seq)
    ell_shared = numeric(length(ell_seq)-1)
    for (i in 1:(length(ell_seq)-1))
    {
      ell_shared[i] = f_ell(within,ell_seq[i],ell_seq[i+1])
    }
    return(ell_shared)
  }
  
  get_f_inf = function(within, ell_seq)
  {
    cutoff = min(ell_seq)
    ell_shared = numeric(length(ell_seq)-1)
    for (i in 1:(length(ell_seq)-1))
    {
      ell_shared[i] = f_ell_inf_chr(within,ell_seq[i],ell_seq[i+1])
    }
    return(ell_shared)
  }
  
  ell_dist_within = get_dist(T,ell_seq)
  ell_dist_between = get_dist(F,ell_seq)
  
  ell_shared_within = get_f(T,ell_seq)
  ell_shared_between = get_f(F,ell_seq)
  
  ell_shared_within_inf = get_f_inf(T,ell_seq)
  ell_shared_between_inf = get_f_inf(F,ell_seq)
  
  return(list(dist_within=ell_dist_within, dist_between=ell_dist_between, shared_within = ell_shared_within, shared_between = ell_shared_between, shared_within_inf = ell_shared_within_inf, shared_between_inf = ell_shared_between_inf))
}

# ****************************************************

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

res = calc_ell_dist(ells,get_N,rates,get_q)
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
