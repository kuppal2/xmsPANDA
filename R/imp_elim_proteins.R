imp_elim_proteins <-
function(nosib.miss, treatment){
  # Impute 1 variable proteins with 1 grp completely missing
  #
  # INPUT: proteins - 1 variable proteins with 1 grp missing completely (NaNs)
  #        grps - grouping information for observations, can be 2+ grps here,
  #               but not sure what to do with more then 2 groups if 2+ gs are
  #               missing completely. So only 1 grp is completely missing.
  # OUTPUT: prs - imputed proteins, same size as proteins parameter
  #nosib.miss <= nosib.miss
  #treatment <= treatment
  tlvls = unique(treatment)
  proteins = nosib.miss
  ng = length(unique(treatment))
  for (i in 1:nrow(proteins)) {
    pr = as.vector(proteins[i,])
    #% find number of missing values in ALL groups
    miss = array(NA, c(1, ng))
    for (j in 1:ng) miss[j] = sum(is.na(pr[treatment==unique(treatment)[j]]))
    #    pos = miss==0
    pos = miss==min(miss)  # Tom 092510
    present_groups = pr[treatment==unique(treatment)[pos]]
    
    # compute mean and stdev from one of the present groups, make sure no NaNs are used
    pospres = !is.na(present_groups)
    presvals = present_groups[pospres]
    pepmean = mean(presvals)
    pepstd =  sd(presvals)
    if(is.na(pepstd)) next;
    
    #% imputing only COMPLETELY missing variables here
    for (j in 1:ng) {
      if   (!pos[j]) { #% imute only the ones not at pos complete
        imppos = is.na(pr[treatment==tlvls[j]])  #% should be all in this group, but double check
        imppepmean = pepmean - 6* pepstd
        imppepstd = pepstd
        tt = imppepmean - 3 * imppepstd
        kk = 0
        while (tt < 0 && kk < 10){  # added kk counter - tom 092510
          offset = .25
          gap = imppepstd * offset
          imppepstd = imppepstd * (1- offset)
          imppepmean = imppepmean + 3 * gap
          tt = imppepmean - 3 * imppepstd
          kk = kk + 1
        }
        imp_tmp = rnorm(length(imppos), imppepmean, pepstd)
        pr[treatment==tlvls[j]] = imp_tmp
      }
    }
    proteins[i,] = pr
  }
  
  # tom - this routine gives some nearly-blank rows
  # to avoid singularities, i'm going to scan this to see which protein rows
  # have all blanks in one row and remove them
  #  proteins <= proteins
  xx = (!is.na(proteins)) %*% model.matrix(~treatment-1)
  notblank.idx = rep(TRUE, nrow(xx))
  for(jj in 1:ncol(xx)) notblank.idx = notblank.idx & xx[,jj]
  #  proteins = proteins[blank.idx,,drop=FALSE]
  
  return(list(proteins=proteins, notblank.idx=notblank.idx))
}
