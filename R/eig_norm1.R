eig_norm1 <-
function (m, treatment, prot.info, write_to_file = "",plot.trend=FALSE)
{
  
  if (is.factor(treatment)) {
    n.treatment = 1
    n.u.treatment = length(unique(treatment))[1]
  }
  else {
    n.treatment = dim(treatment)[2]
    n.u.treatment = dim(unique(treatment))[1]
  }
 # print(n.u.treatment)
  
  m = as.matrix(m)
  grpFactors = treatment
  nGrpFactors = n.treatment
  if (nGrpFactors > 1) {
    ugrps = unique(grpFactors)
    udims = dim(ugrps)
    grp = NULL
    for (ii in seq_len(udims[1])) {
      pos = grpFactors[, 1] == ugrps[ii, 1]
      for (jj in seq(2, udims[2])) {
        pos = pos & grpFactors[, jj] == ugrps[ii, jj]
      }
      grp[pos] = rep(ii, sum(pos))
    }
    grp = as.factor(grp)
  }
  else {
    grp = treatment
  }
  nobs = array(NA, c(nrow(m), length(unique(grp))))
  message("Treatment groups: ", grp)
  
  for (ii in seq_len(nrow(m))) {
    for (jj in seq_len(length(unique(grp)))) {
      nobs[ii, jj] = sum(!is.na(m[ii, grp == unique(grp)[jj]]))
    }
  }
  present.min = matrixStats::rowMins(nobs)
  ii = present.min == 0
  pmiss = rbind(m[ii, ])
  rownames(pmiss) = prot.info[ii, 1]
  present = prot.info[!(prot.info[, 1] %in% rownames(pmiss)),
  ]
  pres = m[!(prot.info[, 1] %in% rownames(pmiss)), ]
  rownames(pres) = prot.info[!(prot.info[, 1] %in% rownames(pmiss)),
                             1]
  #   message("Selecting complete data")
  nobs = array(NA, nrow(pres))
  numiter = nrow(pres)
  nobs = matrixStats::rowSums2(!is.na(pres))
  iii = nobs == ncol(pres)
  complete = rbind(pres[iii, ])
  if (write_to_file != "") {
    utils::write.table(complete, file = write_to_file, append = FALSE,
                       quote = FALSE, sep = "\t", eol = "\n", na = "NaN",
                       dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape",
                                                                                  "double"))
  }
  if (n.u.treatment > 1) {
    #   message("Got 2+ treatment grps")
    lm.fm = makeLMFormula(treatment, "TREATS")
    TREATS = treatment
    TREATS = data.frame(treatment)
    if (is.factor(treatment)) {
      colnames(TREATS) = "TREATS"
    }
    else {
      colnames(TREATS) = colnames(treatment)
    }
    attach(TREATS)
    mod.c = stats::model.matrix(lm.fm$lm.formula, data = TREATS,
                                eval(parse(text = lm.fm$lm.params)))
    Y.c = as.matrix(complete)
    formula1 = paste("t(Y.c)~", as.character(lm.fm$lm.formula)[2],
                     sep = "")
    TREATS = treatment
    fit_lmAll = stats::lm(eval(parse(text = formula1)))
    R.c = stats::residuals(fit_lmAll)
  }
  else {
    message("Got 1 treatment grp")
    mod.c = as.numeric(t(treatment))
    R.c = t(as.matrix(complete))
    TREATS = treatment
  }
  #    message("Computing SVD, estimating Eigentrends...")
  R.c_center = scale(R.c, center = TRUE, scale = FALSE)
  my.svd = svd(R.c_center)
  temp = my.svd$u
  my.svd$u = my.svd$v
  my.svd$v = temp
  numcompletepep = dim(complete)[1]
  #   message("Number of treatments: ", n.u.treatment)
  
  set.seed(123)
  h.c = sva.id(complete, n.u.treatment, lm.fm = lm.fm)$n.sv
  message("Number of bias trends automatically detected ",h.c)
  complete_center = scale(t(complete), center = TRUE, scale = FALSE)
  
  
  # n.u.treatment
  toplot1 = svd(complete_center)
  temp = toplot1$u
  toplot1$u = toplot1$v
  toplot1$v = temp
  if(plot.trend==TRUE){
    message("Preparing to plot...")
    graphics::par(mfcol = c(3, 2))
    graphics::par(mar = c(2, 2, 2, 2))
    plot.eigentrends(toplot1, "Raw Data")
    plot.eigentrends(my.svd, "Residual Data")
  }
  d = my.svd$d
  ss = d^2
  Tk = signif(ss/sum(ss) * 100, 2)
  retval = list(m = m, treatment = treatment, my.svd = my.svd,
                pres = pres, n.treatment = n.treatment, n.u.treatment = n.u.treatment,
                h.c = h.c, present = present, prot.info = prot.info,
                complete = complete, toplot1 = toplot1, Tk = Tk, ncompl = numcompletepep,
                grp = grp)
  return(retval)
}
