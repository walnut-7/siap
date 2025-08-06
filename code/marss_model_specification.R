marss_model <- function(d1, r = d1, p, model = c("simple", "low-rank")){
  model <- match.arg(model)
  if (model == "simple") {
    r = d1
    # ----------- state process ----------
    B <- rbind(cbind(matrix(0, r*(p-1), r), diag(r*(p-1))), matrix(0, nrow = r, ncol = p*r)) 
    b_ls = as.list(B)
    b_par <- sapply(paste0("b", 1:p), function(i) paste(i, sep = ",", 1:r)) 
    b_ls[sapply(p:1, function(i) ((i-1)*r + (1:r)-1)*r*p + r*(p-1) + 1:r)] <- as.vector(b_par)
    B <- matrix(b_ls, r*p, r*p) # see Notion notes.
    
    U <- "zero"
    
    Q <- matrix(0, nrow = r*p, r*p)
    q_ls <- as.list(rep(0, r*p))
    q_ls[r*(p-1) + 1:r] <- as.list(rep("q", r))
    Q[cbind(1:(r*p), 1:(r*p))] <- q_ls
    Q <- matrix(Q, r*p, r*p)
    
    ## ------- observation process -------
    Z <- matrix(0, d1, d1*p)
    Z[cbind(1:d1, 1:d1)] <- 1
    
    A <- "zero" # zero mean
    
    R <- matrix(0, d1, d1)
    
    xi <- "zero" # mean of the initial state, X_0
    Lam <- "diagonal and equal" # "diagonal and equal", "identity" # variance of the initial state
    
    model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,x0=xi,V0=Lam,tinitx=0)
    
    return(model.list)
  }
  if (model == "low-rank"){
    # ----------- state process ----------
    B <- rbind(cbind(matrix(0, r*(p-1), r), diag(r*(p-1))), matrix(0, nrow = r, ncol = p*r)) 
    b_ls = as.list(B)
    b_par <- sapply(paste0("b", 1:p), function(i) paste(i, sep = ",", 1:r)) 
    b_ls[sapply(p:1, function(i) ((i-1)*r + (1:r)-1)*r*p + r*(p-1) + 1:r)] <- as.vector(b_par)
    B <- matrix(b_ls, r*p, r*p) # see Notion notes.
    
    U <- "zero"
    
    Q <- matrix(0, nrow = r*p, r*p)
    q_ls <- as.list(rep(0, r*p))
    q_ls[r*(p-1) + 1:r] <- as.list(rep("q", r))
    Q[cbind(1:(r*p), 1:(r*p))] <- q_ls
    Q <- matrix(Q, r*p, r*p)
    
    ## ------- observation process -------
    Z <- matrix(0, d1, r*p)
    z_ls <- as.list(as.vector(sapply(paste0("z", 1:d1), function(i) paste(i, sep = "," ,1:r))))
    z_ls[1] <- 1 # IC
    z <- matrix(z_ls, byrow = T, nrow = d1, ncol = r) 
    Z[1:d1, 1:r] <- z
    Z <- matrix(Z, nrow = d1, r*p)
    
    A <- "zero" # zero mean
    
    rr_ls <- as.list(rep("r", d1))
    rr <- matrix(0, nrow = d1, ncol = d1)
    rr[cbind(1:d1, 1:d1)] <- rr_ls
    R <- matrix(rr, nrow = d1, ncol = d1)

    xi <- "zero" # mean of the initial state, X_0
    Lam <- "diagonal and equal" # "diagonal and equal", "identity" # variance of the initial state
    
    model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,x0=xi,V0=Lam,tinitx=1)
    
    return(model.list)
  }
}