get_scatter_plots_vold1 <-
function (X = NA, Y = NA, feature_table_file, parentoutput_dir, 
                             class_labels_file, group_by_mat_file = NA, scatterplot.col.opt = "journal", 
                             alphacol = 0.3, newdevice = TRUE, cex.plots = 0.6, replace.by.NA = FALSE, 
                             pairedanalysis = FALSE, filename = "", ylabel = "Response", 
                             alphabetical.order = FALSE, name = NA, add.jitter = TRUE, 
                             add.pvalues = TRUE, xlabel = "Predictor", ellipse = FALSE, 
                             ypos.adj.factor = 0.5, group_by_mat = NA, cor.method = "pearson", multiple.facets=TRUE,
                             ...) 
{
  suppressMessages(library(ggpubr))
  if (is.na(X[1]) == TRUE) {
    data_matrix <- read.table(feature_table_file, sep = "\t", 
                              header = TRUE)
  }
  else {
    data_matrix <- X
  }
  if (is.na(Y[1]) == TRUE) {
    classlabels <- read.table(class_labels_file, sep = "\t", 
                              header = TRUE)
  }
  else {
    classlabels <- Y
  }
  if (is.na(group_by_mat_file[1]) == FALSE) {
    group_by_mat <- read.table(group_by_mat_file, sep = "\t", 
                               header = TRUE)
  }
  class_labels_levels <- c("A")
  sample.col.opt = scatterplot.col.opt
  if (sample.col.opt == "default") {
    col_vec <- c("#CC0000", "#AAC000", "blue", 
                 "mediumpurple4", "mediumpurple1", "blueviolet", 
                 "cornflowerblue", "cyan4", "skyblue", 
                 "darkgreen", "seagreen1", "green", 
                 "yellow", "orange", "pink", "coral1", 
                 "palevioletred2", "red", "saddlebrown", 
                 "brown", "brown3", "white", "darkgray", 
                 "aliceblue", "aquamarine", "aquamarine3", 
                 "bisque", "burlywood1", "lavender", 
                 "khaki3", "black")
  }
  else {
    if (sample.col.opt == "topo") {
      col_vec <- topo.colors(length(class_labels_levels), 
                             alpha = alphacol)
    }
    else {
      if (sample.col.opt == "heat") {
        col_vec <- heat.colors(length(class_labels_levels), 
                               alpha = alphacol)
      }
      else {
        if (sample.col.opt == "rainbow") {
          col_vec <- rainbow(length(class_labels_levels), 
                             start = 0, end = alphacol)
        }
        else {
          if (sample.col.opt == "terrain") {
            col_vec <- cm.colors(length(class_labels_levels), 
                                 alpha = alphacol)
          }
          else {
            if (sample.col.opt == "colorblind") {
              if (length(class_labels_levels) < 9) {
                col_vec <- c("#0072B2", "#E69F00", 
                             "#009E73", "#56B4E9", "#D55E00", 
                             "#CC79A7", "#E64B35FF", 
                             "grey57")
              }
              else {
                col_vec <- c("#0072B2", "#E69F00", 
                             "#009E73", "#56B4E9", "#D55E00", 
                             "#CC79A7", "#E64B35B2", 
                             "#4DBBD5B2", "#00A087B2", 
                             "#3C5488B2", "#F39B7FB2", 
                             "#8491B4B2", "#91D1C2B2", 
                             "#DC0000B2", "#7E6148B2", 
                             "#374E55B2", "#DF8F44B2", 
                             "#00A1D5B2", "#B24745B2", 
                             "#79AF97B2", "#6A6599B2", 
                             "#80796BB2", "#0073C2B2", 
                             "#EFC000B2", "#868686B2", 
                             "#CD534CB2", "#7AA6DCB2", 
                             "#003C67B2", "grey57")
              }
            }
            else {
              check_brewer <- grep(pattern = "brewer", 
                                   x = sample.col.opt)
              if (length(check_brewer) > 0) {
                sample.col.opt_temp = gsub(x = sample.col.opt, 
                                           pattern = "brewer.", replacement = "")
                col_vec <- colorRampPalette(brewer.pal(10, 
                                                       sample.col.opt_temp))(length(class_labels_levels))
              }
              else {
                if (sample.col.opt == "journal") {
                  col_vec <- c("#0072B2", "#E69F00", 
                               "#009E73", "#56B4E9", 
                               "#D55E00", "#CC79A7", 
                               "#E64B35FF", "#3C5488FF", 
                               "#F39B7FFF", "#8491B4FF", 
                               "#91D1C2FF", "#DC0000FF", 
                               "#B09C85FF", "#5F559BFF", 
                               "#808180FF", "#20854EFF", 
                               "#FFDC91FF", "#B24745FF", 
                               "#374E55FF", "#8F7700FF", 
                               "#5050FFFF", "#6BD76BFF", 
                               "#E64B3519", "#4DBBD519", 
                               "#631879E5", "grey75")
                  if (length(class_labels_levels) < 8) {
                    col_vec <- c("#0072B2", "#E69F00", 
                                 "#009E73", "#56B4E9", 
                                 "#D55E00", "#CC79A7", 
                                 "grey75")
                  }
                  else {
                    if (length(class_labels_levels) <= 
                        28) {
                      col_vec <- c("#0072B2", "#E69F00", 
                                   "#009E73", "#56B4E9", 
                                   "#D55E00", "#CC79A7", 
                                   "#E64B35FF", "#3C5488FF", 
                                   "#F39B7FFF", "#8491B4FF", 
                                   "#91D1C2FF", "#DC0000FF", 
                                   "#B09C85FF", "#5F559BFF", 
                                   "#808180FF", "#20854EFF", 
                                   "#FFDC91FF", "#B24745FF", 
                                   "#374E55FF", "#8F7700FF", 
                                   "#5050FFFF", "#6BD76BFF", 
                                   "#8BD76BFF", "#E64B3519", 
                                   "#9DBBD0FF", "#631879E5", 
                                   "#666666", "grey75")
                    }
                    else {
                      colfunc <- colorRampPalette(c("#0072B2", 
                                                    "#E69F00", "#009E73", 
                                                    "#56B4E9", "#D55E00", 
                                                    "#CC79A7", "grey75"))
                      col_vec <- colfunc(length(class_labels_levels))
                      col_vec <- col_vec[sample(col_vec)]
                    }
                  }
                }
                else {
                  if (length(sample.col.opt) == 1) {
                    col_vec <- rep(sample.col.opt, length(class_labels_levels))
                  }
                  else {
                    colfunc <- colorRampPalette(sample.col.opt)
                    col_vec <- colfunc(length(class_labels_levels))
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  rm(X)
  rm(Y)
  cnames <- colnames(data_matrix)
  cnames <- tolower(cnames)
  cnames <- tolower(cnames)
  check_names <- grep(cnames, pattern = "^name$")
  if (length(check_names) > 0) {
    if (check_names == 1) {
      check_names1 <- grep(cnames, pattern = "^mz$")
      check_names2 <- grep(cnames, pattern = "^time$")
      if (length(check_names1) < 1 & length(check_names2) < 
          1) {
        mz <- seq(1.00001, nrow(data_matrix) + 1, 1)
        time <- seq(1.01, nrow(data_matrix) + 1, 1)
        check_ind <- gregexpr(cnames, pattern = "^name$")
        check_ind <- which(check_ind > 0)
        X <- as.data.frame(data_matrix)
        Name <- as.character(X[, check_ind])
        X <- cbind(mz, time, X[, -check_ind])
        names_with_mz_time = cbind(Name, mz, time)
        names_with_mz_time <- as.data.frame(names_with_mz_time)
        write.table(names_with_mz_time, file = "Stage1/Name_mz_time_mapping.txt", 
                    sep = "\t", row.names = FALSE)
        X <- as.data.frame(X)
      }
      else {
        if (length(check_names1) > 0 & length(check_names2) > 
            0) {
          check_ind <- gregexpr(cnames, pattern = "^name$")
          check_ind <- which(check_ind > 0)
          Name <- as.character(X[, check_ind])
          X <- X[, -check_ind]
          names_with_mz_time = cbind(Name, X$mz, X$time)
          colnames(names_with_mz_time) <- c("Name", 
                                            "mz", "time")
          names_with_mz_time <- as.data.frame(names_with_mz_time)
          X <- as.data.frame(X)
          write.table(names_with_mz_time, file = "Stage1/Name_mz_time_mapping.txt", 
                      sep = "\t", row.names = FALSE)
        }
      }
    }
  }
  else {
    check_names1 <- grep(cnames[1], pattern = "^mz$")
    check_names2 <- grep(cnames[2], pattern = "^time$")
    if (length(check_names1) < 1 || length(check_names2) < 
        1) {
      stop("Invalid feature table format. First two columns should be mz and time. Please check example files.")
    }
  }
  dir.create(parentoutput_dir,showWarnings = FALSE)
  setwd(parentoutput_dir)
  data_matrix <- X
  data_m <- data_matrix[, -c(1:2)]
  data_m <- as.matrix(data_m)
  mzvec <- data_matrix[, 1]
  timevec <- data_matrix[, 2]
  goodfeats <- data_m
  rm(data_m)
  file_ind <- 0
  boxplots_fname <- paste("scatterplots", filename, ".pdf", 
                          sep = "")
  if (newdevice == TRUE) {
    pdf(boxplots_fname)
  }
  par(mfrow = c(1, 1), family = "sans", cex = cex.plots)
  class_vec <- classlabels[, 2]
  lapply(1:dim(goodfeats)[1], function(m) {
    if (m%%9 == 0) {
      file_ind <- file_ind + 1
      boxplots_fname <- paste("scatterplots_file", 
                              file_ind, ".tiff", sep = "")
    }
    round_mzval <- sprintf("%.4f", mzvec[m])
    round_timeval <- sprintf("%.1f", timevec[m])
    if (is.na(name) == TRUE) {
      if (length(check_names) > 0) {
        if (check_names == 1) {
          mzname <- as.character(names_with_mz_time[m, 
                                                    1])
        }
        else {
          mzname <- paste("mz_time: ", round_mzval, 
                          "_", round_timeval, sep = "")
        }
      }
      else {
        mzname <- paste("mz_time: ", round_mzval, 
                        "_", round_timeval, sep = "")
      }
    }
    else {
      mzname = as.character(name[m])
    }
    cur_d <- new("list")
    feat_vec <- {
    }
    temp_dm <- cbind(colnames(goodfeats), class_vec, as.vector(t(goodfeats[m, 
                                                                           ])))
    temp_dm <- as.data.frame(temp_dm)
    colnames(temp_dm) <- c("SID", "Class", "Feature")
    temp_dm$Class <- as.numeric(as.character(temp_dm$Class))
    temp_dm$Feature <- as.numeric(as.character(temp_dm$Feature))
    par(mfrow = c(1, 1), family = "sans", cex = cex.plots)
    w <- 0.1
    par(omd = c(0, 1 - w, 0, 1), cex.main = 0.7)
    if (is.na(group_by_mat[1]) == FALSE) {
      colnames(group_by_mat) <- c("SID", "GroupBy")
      temp_dm <- merge(temp_dm, group_by_mat, by = "SID")
    }
    temp_dm <- as.data.frame(temp_dm)
    s1 = summary(temp_dm$Feature)
    sadj = (s1[5] - s1[3]) * ypos.adj.factor
    if (is.na(group_by_mat) == TRUE) {
      p <- ggscatter(temp_dm, y = "Class", x = "Feature", 
                     title = mzname, xlab = xlabel, ylab = ylabel, 
                     palette = col_vec[1], col = col_vec[1], shape = 20, 
                     size = 3, add = "reg.line", add.params = list(color = "black", 
                                                                   fill = "lightgray"), conf.int = TRUE, 
      ) + theme(plot.title = element_text(hjust = 0.5, 
                                          size = 10)) + stat_cor(method = cor.method, label.x = 3, 
                                                                 label.y = max(temp_dm$Feature + (sadj)))
    }
    else {
      temp_dm$GroupBy <- as.factor(temp_dm$GroupBy)
      
      if(multiple.facets==TRUE){
        p <- ggscatter(temp_dm, y = "Class", x = "Feature", 
                       title = mzname, xlab = xlabel, ylab = ylabel, 
                       add = "reg.line", color = "GroupBy", facet.by="GroupBy",
                       palette = col_vec, fullrange = TRUE, shape = 20, 
                       size = 3, conf.int = FALSE, ) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                       size = 10)) + stat_cor(method = cor.method, aes(color = GroupBy))
      }else{
        
        
        p <- ggscatter(temp_dm, y = "Class", x = "Feature", 
                       title = mzname, xlab = xlabel, ylab = ylabel, 
                       add = "reg.line", color = "GroupBy", 
                       palette = col_vec, fullrange = TRUE, shape = 20, 
                       size = 3, conf.int = FALSE, ) + theme(plot.title = element_text(hjust = 0.5, 
                                                                                       size = 10)) + stat_cor(method = cor.method, aes(color = GroupBy))
      }
    }
    print(p)
  })
  if (newdevice == TRUE) {
    dev.off()
  }
}
