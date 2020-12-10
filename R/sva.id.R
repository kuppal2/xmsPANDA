sva.id <-
function(dat, n.u.treatment, lm.fm, B = 500, sv.sig = 0.05)
{
  #message("Number of complete variables (and samples) used in SVD")
  #message(dim(dat))
  n = ncol(dat)
  ncomp = n.u.treatment
  # message("Number of treatment groups (in svd.id): ", ncomp)
  # should be true for either case and can be used later
  
  if (ncomp > 1) {
    formula1 = paste('t(dat)~', as.character(lm.fm$lm.formula)[2], sep = '')
    fit_lmAll = stats::lm(eval(parse(text = formula1)))
    res = t(stats::residuals(fit_lmAll))
  } else {
    res = dat
  }
  # center each peptide around zero (subtract its mean across samples)
  res_center = t(scale(t(res), center = TRUE, scale = FALSE))
  
  uu = svd(t(res_center)) # NEED a WRAPPER for t(). the diag is min(n, m)
  temp = uu$u
  uu$u = uu$v
  uu$v = temp
  
  s0 = uu$d
  s0 = s0 ^ 2
  dstat = s0 / sum(s0)
  ndf = length(dstat)
  dstat0 = matrix(0, nrow = B, ncol = ndf) # num samples
  
  #  message("Starting Bootstrap.....")
  # Bootstrap procedure determines the number of significant eigertrends...
  for (ii in seq_len(B)) {
    if (ii %% 50 == 0) {
      #  message('Iteration ', ii)
    }
    res0 = t(apply(res, 1, sample, replace = FALSE)) # regression
    # center each peptide around zero (subtract its mean across samples)
    # note: not changing matrix itself, only centerig what we pass to svd
    res0_center = t(scale(t(res0), center = TRUE, scale = FALSE))
    uu0 = svd(res0_center)
    temp = uu0$u  # why did tom do this??
    uu0$u = uu0$v
    uu0$v = temp
    
    ss0 = uu0$d  # no need for diag....
    ss0 = ss0 ^ 2
    dstat0[ii, ] = ss0 / sum(ss0) # Tk0 in Matlab
  }
  
  # yuliya: check p-values here, Tom had mean value...
  psv = rep(1, n)
  for (ii in seq_len(ndf)) {
    # should this be compared to a mean? ie dstat0[ii,] ?
    posGreater = dstat0[, ii] > dstat[ii]
    psv[ii] = sum(posGreater) / B
  }
  
  # p-values for trends have to be in the monotonically increasing order,
  # set equal to previous one if not the case
  for (ii in 2:ndf) {
    if (psv[(ii - 1)] > psv[ii]) {
      psv[ii] = psv[(ii - 1)]
    }
  }
  nsv = sum(psv <= sv.sig)
  return(list(n.sv = nsv, p.sv = psv))
}
