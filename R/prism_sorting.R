# pivot to prism
pivot_prism <- function(data, dv, id, iv1 = NULL, iv2) {
  I <- levels(data[[id]])
  N_I <- length(I)
  if(!is.null(iv1)) {
    L1 <- levels(data[[iv1]])
    N_L1 <- length(L1)
  }
  L2 <- levels(data[[iv2]])
  N_L2 <- length(L2)
  if(!is.null(iv1)) {
    M <- 
      matrix(nrow = N_L2+1, ncol = N_L1*N_I+1, dimnames = list(c(iv1, L2), c(iv2, rep(I, times = N_L1))))
    M[1,] <- c(iv1, rep(L1, each = N_I))
    M[,1] <- c(iv1, L2)
    for(i in 2:ncol(M)) {
      for(j in 2:nrow(M)) {
        J <- data[data[iv2] == row.names(M)[j] & data[id] == colnames(M)[i] & data[iv1] == M[1,i],][[dv]]
        if(length(J) != 0) {
          M[j,i] <- mean(J)
        }
      }
    }
  } else {
    M <- 
      matrix(nrow = N_L2, ncol = N_I+1, dimnames = list(L2, c(iv2, I)))
    M[,1] <- L2
    for(i in 2:ncol(M)) {
      for(j in 1:nrow(M)) {
        J <- data[data[iv2] == row.names(M)[j] & data[id] == colnames(M)[i],][[dv]]
        if(length(J) != 0) {
          M[j,i] <- mean(J)
          M[[j]] <- as.double(M[[j]])
        }
      }
    }
  }
  return(as.data.frame(M))
}

# make a list of dataframes filtered by each level of a variable
list_by_var <- function(data, var) {
  L <- levels(data[[var]])
  N_L <- length(L)
  out <- vector("list", N_L)
  for(k in 1:N_L) {
    out[[k]] <- data[data[var] == L[k],]
    out[[k]] <- droplevels(out[[k]])
    names(out)[k] <- L[k]
  }
  return(out)
}

# sort a dataframe by a numerical variable stored as a factor
sortby_num_fac <- function(dataframe, factor) {
  dataframe[order(as.numeric(as.character(dataframe[[factor]]))),]
}

pivot_prism_vec <- function(data, dv, id, iv, vars, num_fac = TRUE) {
  L <- length(data)
  M <- length(vars)
  list_out <- vector("list", L*M)
  for(l in 1:L) {
  for(m in 1:M) {
    list_out[[((l*M)-M)+m]] <- 
        pivot_prism(data = data[[l]], dv = vars[m], id = id, iv1 = NULL, iv2 = iv)
    if(num_fac) {
      list_out[[((l*M)-M)+m]] <- 
        sortby_num_fac(list_out[[((l*M)-M)+m]], iv)
    }
    names(list_out)[((l*M)-M)+m] <- paste(names(gen.list)[l], vars[m], sep = "_")
  }
  }
  return(list_out)
}

