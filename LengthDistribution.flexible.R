# ell_seq: The sequence of segment lengths (ell), for which to return results
# N_func: A function that returns the population size for each t, in number of mating pairs
# rates: A vector indicating the consanguinity rates in the past few generations
# q_func: A function that returns the kinship coefficient between mating pairs beyond the max_T_exact_model past generations. Note that the kinship in the first few generations is now set automatically based on the "rates" argument.
# max_T_exact_model: Specify the number of generations during which we run the exact model with the above rates
# max_T_coal: This is the maximal generation to track when computing the distribution of the TMRCA
# Set a value based on the maximal population size (say, 10x of that size, to guarantee coalescence)
# num_runs: Number of times to run the exact Markov chain (for the distribution of TMRCA in the first few generations)
# old.kinship: kinship estimate to use for t>max_T_exact_model
# ibdne: a flag to indicate whether the historic Ne estimates are from IBDNe or not
calc_ell_dist = function(ell_seq,N_func,rates,q_func,max_T_coal=10000,max_T_exact_model=20,num_runs=10000,old.kinship,ibdne)
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
    time = 0
    while (time<max_T_exact_model)
    {
      if (state == unrelated_state)
      {
        curN = N_func(time,historic.ne[[pop]])/2 ### need to divide the Ne by 2 since N is actually the number of couples 
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
  
  # This will generate the frequency of coal times at 1,2,3,...,max_T_exact_model
  exact_mrca_dist = hist(coal_times,0:max_T_exact_model,plot=F)$counts/num_runs
  
  # This is the formula for the kinship coefficient. No need to change.
  kinship = sum(rates * 4^(-(1:n)))
  
  internal_q = function(t)
  {
    if (t<=max_T_exact_model)
      return(kinship)
    else
      return(q_func(t,old.kinship)) 
  }
  
  p_coal = numeric(max_T_coal)
  p_no_coal = 1
  for (t in 1:max_T_coal)
  {
    if(ibdne){
      eff_N = (N_func(t,historic.ne[[pop]])/2) ### IBDNe is actually estimating N*(1-3q) already, so we remove this multiplication by (1-3*internal_q(t)); the division by 2 is because we want the number of mating couples rather than the Ne
    } else {
      eff_N = (N_func(t,historic.ne[[pop]])/2)*(1-3*internal_q(t)) ### Need to divide Ne by 2 to get N, which is the number of mating pairs
    }
    p_coal[t] = p_no_coal * 1/(4*eff_N)
    p_no_coal = p_no_coal * (1-1/(4*eff_N))
  }
  
  p_coal_between = p_coal
  p_coal_within = numeric()
  p_coal_within[1:max_T_exact_model] = exact_mrca_dist
  ptemp = p_coal[(max_T_exact_model+1):max_T_coal]
  # The sum of everything during t>max_T_exact_model should be (1-rapid_count/num_runs), i.e., the     probability of no coalescence during the exact Markov chain simulations
  p_coal_within[(max_T_exact_model+1):max_T_coal] = ptemp/sum(ptemp) * (1-rapid_count/num_runs)
  
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
  
    return(list(ell_dist_within=ell_dist_within, ell_dist_between=ell_dist_between, shared_within = ell_shared_within, shared_between = ell_shared_between, shared_within_inf = ell_shared_within_inf, shared_between_inf = ell_shared_between_inf))
}
