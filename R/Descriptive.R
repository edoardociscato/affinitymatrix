#' Print affinity matrix
#'
#' This function prints the estimated affinity matrix in tex style and stores it in a txt file.
#'
#' For the sake of clarity, we refer to "men" as the observations on one side of the market and to "women" as the observations on the other side.
#'
#' @param res A list corresponding to the output of \code{\link{estimate.affinity.matrix}}, \code{\link{estimate.affintiy.matrix.lowrank}} or \code{\link{estimate.affinity.matrix.unipartite}}.
#' @param path A string indicating the path where to save the txt file. Defaults to current path.
#' @param name A string indicating the name of the txt file. Defaults to \code{"affinity_matrix"}.
#' @param labels_x A vector of strings indicating the names of men's matching variables. Defaults to \code{"Trait k"} for every \code{k} matching variable.
#' @param labels_y A vector of strings indicating the names of women's matching variables. Defaults to \code{"Trait k"} for every \code{k} matching variable.
#' @param pr A probability indicating the two-tailed significance level required for an estimated parameter to be printed in boldface. Defaults to 0.05 and can be set to 0 to avoid printing any estimate in boldface.
#'
#' @return The function saves a txt file in tex style in \code{path}.
#'
#' @export
print.affinity.matrix = function(res,
                                 path = getwd(),
                                 name = "affinity_matrix",
                                 labels_x = paste0("Trait ",1:Kx),
                                 labels_y = paste0("Trait ",1:Ky),
                                 pr = 0.05) {

  # Prepare table
  z_stat = stats::qnorm(1-pr/2)
  A = res$Aopt; sdA = res$sdA
  Kx = dim(A)[1]; Ky = dim(A)[2]
  tabular = matrix("", nrow=2*Kx+1, ncol=Ky+1)
  for (i in 1:(2*Kx+1)) {
    for (j in 1:(Ky+1)) {
      k = i%/%2; l = j-1
      if (i==1 & j==1) num = "\t&"
      if (i==1 & 1<j & j<Ky+1) num = paste0(labels_y[l],"\t&")
      if (i==1 & j==Ky+1) num = paste0(labels_y[l],"\t\\\\\\hline\n")
      if (i>1 & j==1 & i%%2==0) num = paste0(labels_x[k],"\t&")
      if (i>1 & 1<j & j<Ky+1 & i%%2==0) {
        if(abs(A[k,l])/sdA[k,l]>z_stat) {
          num = sprintf("\\textbf{%0.2f}\t&", A[k,l])
        } else {
          num = sprintf("%0.2f\t&", A[k,l])
        }
      }
      if (i>1 & j==Ky+1 & i%%2==0) {
        if(abs(A[k,l])/sdA[k,l]>z_stat) {
          num = sprintf("\\textbf{%0.2f}\t\\\\\n", A[k,l])
        } else {
          num = sprintf("%0.2f\t\\\\\n", A[k,l])
        }
      }
      if (i>1 & j==1 & i%%2==1) num = "\t&"
      if (i>1 & 1<j & j<Ky+1 & i%%2==1) num = sprintf("(%0.2f) \t&", sdA[k,l])
      if (i>1 & j==Ky+1 & i%%2==1) num = sprintf("(%0.2f) \t\\\\\n", sdA[k,l])
      tabular[i,j] = num
    }
  }

  # Save
  utils::write.table(tabular,
                     file = paste0(path,"/",name,".txt"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE)

}

#' Print the diagonal of the affinity matrix
#'
#' This function prints the estimated diagonal elements of the affinity matrix in tex style and stores them in a txt file.
#'
#' For the sake of clarity, we refer to "men" as the observations on one side of the market and to "women" as the observations on the other side.
#'
#' @param res A list corresponding to the output of \code{\link{estimate.affinity.matrix}}, \code{\link{estimate.affintiy.matrix.lowrank}} or \code{\link{estimate.affinity.matrix.unipartite}}.
#' @param path A string indicating the path where to save the txt file. Defaults to current path.
#' @param name A string indicating the name of the txt file. Defaults to \code{"diagonal_affinity_matrix"}.
#' @param labels A vector of strings indicating the names of the matching variables. Defaults to \code{"Trait k"} for every \code{k} matching variable.
#' @param pr A probability indicating the two-tailed significance level required for an estimated parameter to be printed in boldface. Defaults to 0.05 and can be set to 0 to avoid printing any estimate in boldface.
#'
#' @return The function saves a txt file in tex style in \code{path}.
#'
#' @export
print.diagonal = function(res,
                          path = getwd(),
                          name = "diagonal_affinity_matrix",
                          labels = paste0("Trait ",1:K),
                          pr = 0.05) {

  # Prepare table
  z_stat = stats::qnorm(1-pr/2)
  A = diag(res$Aopt); sdA = diag(res$sdA)
  K = length(A)
  tabular = matrix("", nrow=3, ncol=K)
  for (i in 1:3) {
    for (j in 1:K) {
      k = j
      if (i==1 && j<K) num = paste0(labels[k],"\t&")
      if (i==1 && j==K) num = paste0(labels[k],"\t\\\\\\hline\n")
      if (i==2 && j<K) {
        if(abs(A[k])/sdA[k]>z_stat) {
          num = sprintf("\\textbf{%0.2f}\t&", A[k])
        } else {
          num = sprintf("%0.2f\t&", A[k])
        }
      }
      if (i==2 && j==K) {
        if(abs(A[k])/sdA[k]>z_stat) {
          num = sprintf("\\textbf{%0.2f}\t\\\\\n", A[k])
        } else {
          num = sprintf("%0.2f\t\\\\\n", A[k])
        }
      }
      if (i==3 && j<K) num = sprintf("(%0.2f)\t&", sdA[k])
      if (i==3 && j==K) num = sprintf("(%0.2f)\t\\\\\n", sdA[k])
      tabular[i,j] = num
    }
  }

  # Save
  utils::write.table(tabular,
                     file = paste0(path,"/",name,".txt"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE)

}

#' Print summaries of rank tests
#'
#' This function prints the summaries of the first \code{n_tests} rank tests in tex style and stores them in a txt file.
#'
#' For the sake of clarity, we refer to "men" as the observations on one side of the market and to "women" as the observations on the other side.
#'
#' @param res A list corresponding to the output of \code{\link{estimate.affinity.matrix}}, \code{\link{estimate.affintiy.matrix.lowrank}} or \code{\link{estimate.affinity.matrix.unipartite}}.
#' @param path A string indicating the path where to save the txt file. Defaults to current path.
#' @param name A string indicating the name of the txt file. Defaults to \code{"rank_tests"}.
#' @param pr A probability indicating the significance level required to pass a rank test. Defaults to 0.05.
#' @param n_tests An integer indicating the number of tests to print. The function prints the first \code{n_tests} rank tests. Defaults to \code{min(nrow(Y),nrow(X))-1}.
#'
#' @return The function saves a txt file in tex style in \code{path}.
#'
#' @export
print.test = function(res,
                      path = getwd(),
                      name = "rank_tests",
                      pr = .05,
                      n_tests = K-1) {

  # Prepare table
  K = length(res$lambda)
  n_tests = min(n_tests + 1, K)
  lambda = res$lambda/sum(res$lambda); U = res$U; V = res$V
  K = length(lambda)
  tabular = matrix("", nrow=4, ncol=n_tests)
  for (i in 1:4) {
    for (j in 1:n_tests) {
      k = j-1
      if (i==1 && j==1) num = "$H_0$: $rk(A)=k$\t&"
      if (i==1 && 1<j && j<n_tests) num = paste0("$k=",k,"$\t&")
      if (i==1 && j==n_tests) num = paste0("$k=",k,"$\t\\\\\\hline\n")
      if (i==2 && j==1) num = paste0("$\\chi^2$\t&")
      if (i==2 && 1<j && j<n_tests) num = sprintf("%0.2f\t&", res$rank.tests[[k]]$chi2)
      if (i==2 && j==n_tests) num = sprintf("%0.2f\t\\\\\n", res$rank.tests[[k]]$chi2)
      if (i==3 && j==1) num = paste0("$df$\t&")
      if (i==3 && 1<j && j<n_tests) num = sprintf("%d\t&", res$rank.tests[[k]]$df)
      if (i==3 && j==n_tests) num = sprintf("%d\t\\\\\n", res$rank.tests[[k]]$df)
      if (i==4 && j==1) num = "Rejected?\t&"
      if (i==4 && 1<j && j<n_tests) {
        if (res$rank.tests[[k]]$chi2>stats::qchisq(1-pr, df=res$rank.tests[[k]]$df)) {
          num = "Yes\t&"
        } else {num = "No\t&"}
      }
      if (i==4 && j==n_tests) {
        if (res$rank.tests[[k]]$chi2>stats::qchisq(1-pr, df=res$rank.tests[[k]]$df)) {
          num = "Yes\\\\\n"
        } else {num = "No\\\\\n"}
      }
      tabular[i,j] = num
    }
  }

  # Save
  utils::write.table(tabular,
                     file = paste0(path,"/",name,".txt"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE)

}

#' Print summary of saliency analysis
#'
#' This function prints the results from the saliency analysis in tex style and stores them in a txt file.
#'
#' For the sake of clarity, we refer to "men" as the observations on one side of the market and to "women" as the observations on the other side.
#'
#' @param res A list corresponding to the output of \code{\link{estimate.affinity.matrix}}, \code{\link{estimate.affintiy.matrix.lowrank}} or \code{\link{estimate.affinity.matrix.unipartite}}.
#' @param path A string indicating the path where to save the txt file. Defaults to current path.
#' @param name_x A string indicating the name of the txt file storing the results for men. Defaults to \code{"saliency_analysis_x"}.
#' @param name_y A string indicating the name of the txt file storing the results for women. Defaults to \code{"saliency_analysis_y"}.
#' @param ncol_x An integer indicating the number of singular vector to print for men. The function prints the first \code{ncol_x} singular vectors. Defaults to \code{ncol(U)}.
#' @param ncol_y An integer indicating the number of singular vector to print for women. The function prints the first \code{ncol_y} singular vectors. Defaults to \code{ncol(V)}.
#' @param labels_x A vector of strings indicating the names of men's matching variables. Defaults to \code{"Trait k"} for every \code{k} matching variable.
#' @param labels_y A vector of strings indicating the names of women's matching variables. Defaults to \code{"Trait k"} for every \code{k} matching variable.
#' @param pr A probability indicating the two-tailed significance level required for an estimated parameter to be printed in boldface. Defaults to 0.05 and can be set to 0 to avoid printing any estimate in boldface.
#'
#' @return The function saves a txt file in tex style in \code{path}.
#'
#' @export
print.saliency = function(res,
                          path = getwd(),
                          name_x = "saliency_analysis_x",
                          name_y = "saliency_analysis_y",
                          ncol_x = Kx,
                          ncol_y = Ky,
                          labels_x = paste0("Trait ",1:Kx),
                          labels_y = paste0("Trait ",1:Ky),
                          pr = .05) {

    # Normalize s.t. largest weight in every column (in absolut terms) is positive
    lambda = res$lambda/sum(res$lambda);
    U = res$U; V = res$V; K = length(lambda)
    change_sign_x = kronecker(t(apply(U,2,function(x) sign(x[which.max(abs(x))]))),matrix(1,nrow=Kx,ncol=1)) # must be the same for both sides
    change_sign_y = kronecker(t(apply(U,2,function(x) sign(x[which.max(abs(x))]))),matrix(1,nrow=Ky,ncol=1))
    U = change_sign_x*U; V = change_sign_y*V

    # Confidence intervals
    df.bootstrap = res$df.bootstrap
    lambdaCI = matrix(0,nrow=K,ncol=2);
    for (k in 1:K) lambdaCI[k,] = stats::quantile(df.bootstrap[,Kx*Ky+k],c(pr/2, 1-pr/2))
    UCI = matrix(0,nrow=Kx*K,ncol=2);
    for (k in 1:(Kx*K)) UCI[k,] = stats::quantile(df.bootstrap[,Kx*Ky+K+k],c(pr/2, 1-pr/2))
    VCI = matrix(0,nrow=Ky*K,ncol=2);
    for (k in 1:(Ky*K)) VCI[k,] = stats::quantile(df.bootstrap[,Kx*Ky+K+Kx*K+k],c(pr/2, 1-pr/2))
    Ulb = matrix(UCI[,1],nrow=Kx); Uub = matrix(UCI[,2],nrow=Kx)
    Vlb = matrix(VCI[,1],nrow=Ky); Vub = matrix(VCI[,2],nrow=Ky)
    testU = sign(Ulb*Uub)==1; testV = sign(Vlb*Vub)==1

    # Prepare table (men)
    ncol_x = min(ncol_x+1, K+1)
    tabular_m = matrix("", nrow=Kx+2, ncol=ncol_x)
    for (i in 1:(Kx+2)) {
      for (j in 1:ncol_x) {
        k = i-1; l = j-1
        if (i==1 & j==1) num_m = "\t&"
        if (i==1 & 1<j & j<ncol_x) num_m = paste0("Index ",l,"\t&")
        if (i==1 & j==ncol_x) num_m = paste0("Index ",l,"\t\\\\\\hline\n")
        if (1<i & i<Kx+2 & j==1) num_m = paste0(labels_x[k],"\t&")
        if (1<i & i<Kx+2 & 1<j & j<ncol_x) {
          if(testU[k,l]>.1) num_m = sprintf("\\textbf{%0.2f}\t&", U[k,l]) else num_m = sprintf("%0.2f\t&", U[k,l])
        }
        if (1<i & i<Kx+2 & j==ncol_x) {
          if(testU[k,l]>.1) num_m = sprintf("\\textbf{%0.2f}\t\\\\\n", U[k,l]) else num_m = sprintf("%0.2f\t\\\\\n", U[k,l])
        }
        if (i==Kx+2 & j==1) num_m = "\\hline Index share\t&"
        if (i==Kx+2 & 1<j & j<ncol_x) num_m = sprintf("%0.2f\t&", lambda[l])
        if (i==Kx+2 & j==ncol_x) num_m = sprintf("%0.2f\t\\\\\n", lambda[l])
        tabular_m[i,j] = num_m
    }
  }

  # Prepare table (women)
  ncol_y = min(ncol_y+1, K+1)
  tabular_f = matrix("", nrow=Ky+2, ncol=ncol_y)
  for (i in 1:(Ky+2)) {
    for (j in 1:ncol_y) {
      k = i-1; l = j-1
      if (i==1 & j==1) num_f = "\t&"
      if (i==1 & 1<j & j<ncol_y) num_f = paste0("Index ",l,"\t&")
      if (i==1 & j==ncol_y) num_f = paste0("Index ",l,"\t\\\\\\hline\n")
      if (1<i & i<Ky+2 & j==1) num_f = paste0(labels_y[k],"\t&")
      if (1<i & i<Ky+2 & 1<j & j<ncol_y) {
        if(testV[k,l]>.1) num_f = sprintf("\\textbf{%0.2f}\t&", V[k,l]) else num_f = sprintf("%0.2f\t&", V[k,l])
      }
      if (1<i & i<Ky+2 & j==ncol_y) {
        if(testV[k,l]>.1) num_f = sprintf("\\textbf{%0.2f}\t\\\\\n", V[k,l]) else num_f = sprintf("%0.2f\t\\\\\n", V[k,l])
      }
      if (i==Ky+2 & j==1) num_f = "\\hline Index share\t&"
      if (i==Ky+2 & 1<j & j<ncol_y) num_f = sprintf("%0.2f\t&", lambda[l])
      if (i==Ky+2 & j==ncol_y) num_f = sprintf("%0.2f\t\\\\\n", lambda[l])
      tabular_f[i,j] = num_f
    }
  }

  # Save
  utils::write.table(tabular_m,
                     file = paste0(path,"/",name_x,".txt"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
  utils::write.table(tabular_f,
                     file = paste0(path,"/",name_y,".txt"),
                     row.names = FALSE, col.names = FALSE, quote = FALSE)

}

#' Print cross validation summary
#'
#' This function prints a plot where it is shown the covariance mismatch as an estimated function of the rank restriction parameter lambda - the result of the cross validation exercise. The function is expected to be convex in lambda and the chosen lambda is the unique minimum.
#'
#' For the sake of clarity, we refer to "men" as the observations on one side of the market and to "women" as the observations on the other side.
#'
#' @param res A list corresponding to the output of \code{\link{estimate.affinity.matrix}}, \code{\link{estimate.affintiy.matrix.lowrank}} or \code{\link{estimate.affinity.matrix.unipartite}}.
#' @param path A string indicating the path where to save the plot. Defaults to current path.
#' @param name A string indicating the name of the plot. Defaults to \code{"affinity_matrix"}.
#'
#' @return The function saves a png file in \code{path}.
#'
#' @export
print.CV = function(res,
                    path = getwd(),
                    name = "cross_validation_plot") {

  # Plot cross-validation results
  df_mean = aggregate(res$df.cross.validation, list(res$df.cross.validation$lambda), mean)
  ggplot2::ggplot() + ggplot2::geom_line(data = df_mean, ggplot2::aes(lambda, cov_err), size = 1) +
    ggplot2::labs(x="Lambda", y="Covariance mismatch") + ggplot2::scale_color_viridis_d() + ggplot2::theme_minimal() +
    ggplot2::geom_vline(xintercept = df_mean$lambda[which.min(df_mean$cov_err)], color = "red", size = 1, linetype = "dotted")
  ggplot2::ggsave(paste0(path, "/",name,".png"), width=8, height=4)

}

#' Print affinity matrix
#'
#' This function prints a series of plots, one for each orthogonal sorting dimension. In the k-th plot, the correlation between a man's observed matching variable and the man's k-th matching factor is plotted on the x-axis; the correlation between a woman's observed matching variable and the woman's k-th matching factor is plotted on the y-axis. In addition, the user can supply additional variables stored in the matrix \code{Z} that were not previously used in the estimation ("outcome variables"). The function prints the correlation between the outcome variable and the man's k-th matching factor on the x-axis, while the correlation between the outcome variable and the woman's k-th matching factor is on the y-axis.
#'
#' For the sake of clarity, we refer to "men" as the observations on one side of the market and to "women" as the observations on the other side.
#'
#' @param res A list corresponding to the output of \code{\link{estimate.affinity.matrix}}, \code{\link{estimate.affintiy.matrix.lowrank}} or \code{\link{estimate.affinity.matrix.unipartite}}.
#' @param Z A matrix Z with additional variables that were not previously used in the estimation. The i-th row of \code{Z} must contain information on the couple formed by the i-th row of \code{X} and the i-th row of \code{Y}, so that \code{nrow(Z)=nrow(X)}. Defaults to an empty matrix: \code{Z} is optional.
#' @param path A string indicating the path where to save the plots. Defaults to current path.
#' @param name A string indicating a common root for plot names. Defaults to \code{"correlation_plot"}.
#' @param labels_x A vector of strings indicating the names of men's matching variables. Defaults to \code{"Trait k"} for every \code{k} matching variable.
#' @param labels_y A vector of strings indicating the names of women's matching variables. Defaults to \code{"Trait k"} for every \code{k} matching variable.
#' @param labels_z A vector of strings indicating the names of the outcome variables. Defaults to \code{"Outcome k"} for every \code{k} outcome variable.
#' @param ndims An integer indicating the number of orthogonal matching dimensions that will be plotted. The function plots the first \code{ndims} dimensions. Defaults to all dimensions unless the latter are more than 10, in which case only the first 10 are plotted.
#' @param pr A probability indicating the two-tailed significance level required for a matching or outcome variable to be displayed in a plot. In order to avoid having too many variables plotted at the same time, the function only selects those whose correlation with the matching factor is significantly different from zero (in a two-tailed test) at the \code{pr} level. Defaults to 0.02 and can be set to 1 to print all variables.
#' @param color_arrows A string or a vector of strings containing color names for the arrows. All matching variables are assigned the first color given in the vector, while all outcome variables are assigned the second color. See \code{\link[ggplot2]{ggplot}}. Defaults to \code{"black"} and \code{"red"} respectively.
#' @param size_arrows A positive real number or a vector containing the size of the arrows. All matching variables are assigned the first size given in the vector, while all outcome variables are assigned the second size. See \code{\link[ggplot2]{ggplot}}. Defaults to 0.5 for both.
#' @param font_labels A string or a vector of strings containing font types for the labels. All matching variables are assigned the first font type given in the vector, while all outcome variables are assigned the second font type. See \code{\link[ggplot2]{ggplot}}. Defaults to \code{"bold"} and \code{"italic"} respectively.
#' @param label_x_axis A string containing a root for all x-axis names in different plots. Defaults to \code{"First partner"}.
#' @param label_y_axis A string containing a root for all y-axis names in different plots. Defaults to \code{"Second partner"}.
#'
#' @return The function saves one separate png file for each of the \code{ndims} dimensions in \code{path}.
#'
#' @export
print.correlations = function(res,
                              Z = matrix(0, nrow = N, ncol = 0),
                              path = getwd(),
                              name = "correlation_plot",
                              labels_x = paste0("Trait ",1:Kx),
                              labels_y = paste0("Trait ",1:Ky),
                              labels_z = if (Kz>0) paste0("Outcome ",1:Kz) else c(),
                              ndims = min(Kx,Ky,10),
                              pr = 0.02,
                              color_arrows = c("black", "red"),
                              size_arrows = .5,
                              font_labels = c("bold", "italic"),
                              label_x_axis = "First partner",
                              label_y_axis = "Second partner") {

  # Indices
  X = res$X; Y = res$Y; N = nrow(X)
  Kx = ncol(X); Ky = ncol(Y);
  Z = as.matrix(Z); Kz = ncol(Z)
  U = res$U; V = res$V
  if (Kx>=Ky) {
    change_sign_x = kronecker(t(apply(U,2,function(x) sign(x[which.max(abs(x))]))),matrix(1,nrow=Kx,ncol=1)) # must be the same for both sides
    change_sign_y = change_sign_x[1:Ky,1:Ky]
  } else {
    change_sign_y = kronecker(t(apply(V,2,function(x) sign(x[which.max(abs(x))]))),matrix(1,nrow=Ky,ncol=1)) # must be the same for both sides
    change_sign_x = change_sign_y[1:Kx,1:Kx]
  }
  U = change_sign_x*U
  V = change_sign_y*V
  Ix = X%*%U; Iy = Y%*%V

  # Compute correlations and plot them
  K = max(Kx,Ky); ndims = min(ndims, K)
  if (Kx>=Ky) labels = labels_x else labels = labels_y
  if (length(color_arrows)==1) color_arrows = rep(color_arrows, 2)
  if (length(size_arrows)==1) size_arrows = rep(size_arrows, 2)
  if (length(font_labels)==1) font_labels = rep(font_labels, 2)
  for (d in 1:ndims) {
    cor = data.frame(index = 1:(K+Kz),
                     label = c(labels, labels_z))
    IxX = Hmisc::rcorr(cbind(Ix[,d],X))
    IyY = Hmisc::rcorr(cbind(Iy[,d],Y))
    IxZ = Hmisc::rcorr(cbind(Ix[,d],Z))
    IyZ = Hmisc::rcorr(cbind(Iy[,d],Z))
    cor$x = c(IxX$r[1,-1], rep(0, max(Ky-Kx,0)), IxZ$r[1,-1])
    cor$px = c(IxX$P[1,-1], rep(1, max(Ky-Kx,0)), IxZ$P[1,-1])
    cor$y = c(IyY$r[1,-1], rep(0, max(Kx-Ky,0)), IyZ$r[1,-1])
    cor$py = c(IyY$P[1,-1], rep(1, max(Kx-Ky,0)), IyZ$P[1,-1])
    cor$filter = cor$px<pr | cor$py<pr
    cor$origin = 0;
    cor$color = c(rep(color_arrows[1],K),rep(color_arrows[2],Kz))
    cor$size = c(rep(size_arrows[1],K),rep(size_arrows[2],Kz))
    cor$Variable = c(rep("Matching",K),rep("Outcome",Kz))
    cor$font = c(rep(font_labels[1],K),rep(font_labels[2],Kz))
    cor$hjust = .5 + cor$x/2 #ifelse(cor$x>0, 0, 1)
    cor$vjust = ifelse(cor$y>0, 1, 0)
    circle_dat = data.frame(tt = seq(0,2*pi,length.out = 5000))
    circle_dat$xx = cos(circle_dat$tt); circle_dat$yy = sin(circle_dat$tt)
    ggplot2::ggplot() + ggplot2::coord_fixed(ratio = 1) + ggplot2::theme_minimal() +
      ggplot2::theme(axis.text=ggplot2::element_text(size=8),
                     axis.title=ggplot2::element_text(size=8),
                     legend.title=ggplot2::element_text(size=8),
                     legend.text=ggplot2::element_text(size=8)) +
      ggplot2::geom_path(data = circle_dat, ggplot2::aes(xx,yy), color="grey") +
      ggplot2::geom_vline(xintercept = 0, size=.5, color="grey") + ggplot2::geom_hline(yintercept = 0, size=.5, color="grey") +
      ggplot2::geom_segment(data = cor[cor$filter,],
                            ggplot2::aes(x=origin,y=origin,colour=Variable,size=size,xend=x,yend=y),
                            size = cor$size[cor$filter], color = cor$color[cor$filter],
                            arrow = grid::arrow(length=grid::unit(cor$size[cor$filter]*0.015, "npc"), type="closed")) +
      ggrepel::geom_text_repel(data = cor[cor$filter,],
                               ggplot2::aes(x=x,y=y,label=label,hjust=hjust,vjust=vjust,fontface=font),
                               size=2, segment.size=0.25, seed=1) +
      ggplot2::labs(x=paste0(label_x_axis,", dim. ",d), y=paste0(label_y_axis,", dim. ",d))
    ggplot2::ggsave(paste0(path, "/",name,"_dim_",d,".png"), width=8, height=4)
  }

}

