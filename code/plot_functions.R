library(fields)
panel_image <- function(X, yleg=NULL, xleg=NULL, title=NULL, legend = "SSI(mW/m2/nm)\n", zlim = range(X, na.rm = T),
                        ylab = "Wavelength (nm)", xlab = "Time (day)", text = F, col = tim.colors(64),
                        xbreaks = 10, ybreaks = 10,
                        save=F, file=NULL, width=12, height=9, units = 'in', res = 300, pointsize=10,
                        legend.strip = T, legend.side = 3-2*horizontal, horizontal = F, ...) {
  d1=nrow(X)
  d2=ncol(X)
  if (is.null(title)) {
    title = paste0(deparse(substitute(X)),"\n\n")
  }
  if (legend.strip == F) {
    par(cex=1, cex.axis = 1, cex.lab  = 1)
    image(matrix(t(X)[,d1:1], nrow = d2, ncol = d1), axes=FALSE, xlab="", ylab="",
          main=title, zlim = zlim, col = tim.colors())
    if (!is.null(xleg)) {
      ind = c(1+0:(xbreaks-1)*floor(d2/xbreaks), d2)
      axis(3, at = ind/d2, labels = round(xleg[ind], 1), hadj = 0.5, padj = 1)
    }
    axis(3, at = 0.5, labels = xlab, tick = F)
    if (!is.null(yleg)) {
      ind = c(d1,1+(ybreaks-1):0*floor(d1/ybreaks))
      axis(2, at = rev(ind/d1), labels = round(yleg[ind],1), las=2, hadj = 0.8)
    }
    axis(2, at = 0.5, labels = ylab, tick = F, padj = -3)
    if (text==T) {
      m =  matrix(t(X)[,d1:1], nrow = d2, ncol = d1)
      for (i in 1:d2) {
        for (j in 1:d1) {
          text((i-1.5)/(d2-1), (j-1)/(d1-1), format(m[i, j], digits=2), col = "black", cex = 1, pos = 4)
        }
      }
    }      
  } else {
    par(cex=1, cex.axis = 1, cex.lab  = 1)
    image.plot(matrix(t(X)[,d1:1], nrow = d2, ncol = d1), axes=FALSE, xlab="", ylab="", horizontal = horizontal,
               legend.args = list(text=legend, side=legend.side, line=1.2), main=title, zlim = zlim,
               col = col, ...)
    if (!is.null(xleg)) {
      ind = c(1+0:(xbreaks-1)*floor(d2/xbreaks), d2)
      axis(3, at = ind/d2, labels = round(xleg[ind], 1), hadj = 0.5, padj = 1)
    }
    axis(3, at = 0.5, labels = xlab, tick = F)
    #axis(1, at = seq(0, 1, 0.1), labels = rep("",11))
    if (!is.null(yleg)) {
      ind = c(d1,1+(ybreaks-1):0*floor(d1/ybreaks))
      axis(2, at = rev(ind/d1), labels = round(yleg[ind],1), las=2, hadj = 0.8)
    }
    axis(2, at = 0.5, labels = ylab, tick = F, padj = -3)
    if (text==T) {
      m =  matrix(t(X)[,d1:1], nrow = d2, ncol = d1)
      for (i in 1:d2) {
        for (j in 1:d1) {
          text((i-1.5)/(d2-1), (j-1)/(d1-1), format(m[i, j], digits=2), col = "black", cex = 1, pos = 4)
        }
      }
    }
  }
  if (save == T) {
    dev.copy(png, file=paste0(file, ".png"), width=width, height=height, units=units, res=res, pointsize=pointsize)
    dev.off()
  }
}

row_plot <- function(x, i, x_imp=NULL, xt=NULL,  
                     xlab = "Time", ylab = "SSI", 
                     save=F, file=NULL, width=8, height=6, units = 'in', res = 300, pointsize=10) {
  x_coordinates = which(apply(x, 2, function(t) all(is.na(t))))
  par(cex=1.2, cex.axis = 1.2, cex.lab  = 1.2)
  plot(x[i,], type = "l", lty = 1, xlab = xlab, ylab = ylab)
  for (j in seq_along(x_coordinates)) {
    rect(xleft = x_coordinates[j] - 0.5, xright = x_coordinates[j] + 0.5, ybottom = min(x[i,], na.rm=T), ytop = max(x[i,], na.rm=T), col = "gray", border = NA)
  } # downtime
  
  if (!is.null(xt)) {
    lines(xt[i,], lty = 2)
    points(which(is.na(x[i, ])), xt[i, which(is.na(x[i, ]))], pch = 16, cex = 0.8, col = "orange") # missing
  }
  if (!is.null(x_imp)) {
    lines(x_imp[i,], col = "blue") 
  }
  if (save == T) {
    dev.copy(png, file=paste0(file, ".png"), width=width, height=height, units=units, res=res, pointsize=pointsize)
    dev.off()
  }
}

weights <- function(x, t, lambda2=0, lambda3=0, alpha=0, beta, a, b) {
  # - beta: p x r
  p = nrow(beta)
  beta = rbind(-1, beta) # (p+1) x r
  r = ncol(b)
  d2 = ncol(x)
  
  Ome_prime = apply(x, 2, function(v) all(is.na(v))) # The indicator set of columns that are entirely missing
  
  Ome_p_1 = Ome_prime # downtime left p neighbor indicator set. Length: n-p
  for (i in which(Ome_prime==1)){
    l = max(1, i-p)
    #h = min(length(Ome_prime), i+p)
    #Ome_prime[l:h] = 1 # also include column within the p neighborhood
    Ome_p_1[l:i] = 1
  }
  S1 = which(Ome_p_1 == 1)

  Ome_p_2 = Ome_prime # downtime right p neighbor indicator set. Length: n-p
  for (i in which(Ome_prime==1)){
    #l = max(1, i-p)
    h = min(d2, i+p)
    #Ome_prime[l:h] = 1 # also include column within the p neighborhood
    Ome_p_2[i:h] = 1
  }
  S2 = which(Ome_p_2==1)

  S = which(Ome_prime==1)
  
  comp <- function(S) {
    return(setdiff(1:d2, S))
  }
  
  crossprod_Gamma <- function(t,i) {
    # - t: t \in [1,d2]
    # - i: i \in [1,p] or [-p,-1]
    l = max(intersect(c(0, comp(union(S1, S2))), 0:(t-1))) # l \in [0,d2-1]
    u = min(intersect(c(comp(union(S1, S2)), d2+1), (t+1):(d2+1))) # u \in [2,d2+1]
    if (max(i, p-(t-l-1), 0) > min(u-t-1, p+i, p)) {
      return(0)
    }
    rowSums(matrix(sapply(max(i, p-(t-l-1), 0):min(u-t-1, p+i, p) + 1, function (j) beta[j,] * beta[j-i,]), nrow = r))
  }
  
  W = matrix(0, nrow = r, ncol = ncol(x))
  
  if (t %in% setdiff(S, 1:p)) {
    for (i in 1:p) {
      W[, t+i] = - alpha * diag(t(a) %*% a + alpha * diag(colSums(beta^2), nrow = r))^(-1) * crossprod_Gamma(t,i) 
      W[, t-i] = W[, t+i]
    }
    W[, t] = diag(t(a) %*% a + alpha * diag(colSums(beta^2), nrow = r))^(-1) * diag(matrix(t(a) %*% a)) # weight on b_t in the last iteration
  } else if (t %in% intersect(S1, union(1:p, comp(S2)))) {
    l = max(intersect(c(0, comp(S1)), 0:(t-1))) # l \in [0,d2-1]
    for (i in 1:p) {
      W[, t+i] = - alpha * diag(t(a) %*% a + lambda3 * diag(r) + alpha * diag(colSums(matrix(beta[(p-(t-l-1)):p + 1,]^2)), nrow = r))^(-1) * crossprod_Gamma(t,i) 
      W[, t-i] = - alpha * diag(t(a) %*% a + lambda3 * diag(r) + alpha * diag(colSums(matrix(beta[(p-(t-l-1)):p + 1,]^2)), nrow = r))^(-1) * crossprod_Gamma(t,-i)
    }
    W[, t] = diag(t(a) %*% a + lambda3 * diag(r) + alpha * diag(colSums(matrix(beta[(p-(t-l-1)):p + 1,]^2)), nrow = r))^(-1) * diag(matrix(t(a) %*% a)) # weight on b_t in the last iteration
  } else if (t %in% setdiff(setdiff(S2, S), 1:p)) {
    u = min(intersect(c(comp(S2), d2+1), (t+1):(d2+1))) # u \in [2,d2+1]
    for (i in 1:p) {
      W[, t+i] = - alpha * diag(t(a) %*% a + alpha * diag(colSums(matrix(beta[0:(u-t-1) + 1,]^2)), nrow = r))^(-1) * crossprod_Gamma(t,i)
      W[, t-i] = - alpha * diag(t(a) %*% a + alpha * diag(colSums(matrix(beta[0:(u-t-1) + 1,]^2)), nrow = r))^(-1) * crossprod_Gamma(t,-i)
    }
    W[, t] = diag(t(a) %*% a + alpha * diag(colSums(matrix(beta[0:(u-t-1) + 1,]^2)), nrow = r))^(-1) * diag(matrix(t(a) %*% a)) # weight on b_t in the last iteration
  } else {
    W[, t] = diag(t(a) %*% a + lambda2 * diag(nrow = r))^(-1) * diag(matrix(t(a) %*% a))
  }
  
  return(W)
  
}

err <- function(x_imp, x=NULL, o=is.na(x), xt, type=c("msre", "rmse", "rse")) { 
  # - o: the indicator matrix of NA's in x
  type = match.arg(type)
  o[which(o==0)] = NA # Does not calculate err on the obverved entries
  if (type == "rmse") { # relative mean squared error
    return(apply((o*(x_imp - xt))^2, 1, mean, na.rm = T) / apply((xt*o)^2, 1, sum, na.rm = T))
  }
  if (type == "rse") { # relative squared error
    return(apply((o*(x_imp - xt))^2, 1, sum, na.rm = T) / apply((xt*o)^2, 1, sum, na.rm = T))
  }
  if (type == "msre") { # mean relative squared error
    return(apply((o*(x_imp - xt)/xt)^2, 1, mean, na.rm = T))
  }
}
