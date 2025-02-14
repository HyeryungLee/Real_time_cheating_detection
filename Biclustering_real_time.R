################################################################
##### Real-Time Cheating Detection on Mixed-Format Tests #######
################ Using a Biclustering Approach #################
############ Hyeryung Lee & Walter P. Vispoel (2025) ###########
################################################################

#### Library
library(QUBIC)
library(DescTools)

#### Data List (Yor own data)
key_data # vector of answer key for each item
resp_data # matrix of response answer data collected at a specific time point
time_data # matrix of response time data collected at a specific time point

time_point # a specific time point at which cheating detection is conducted 
# Identify the last item where the cumulative response time remains within time_point
find_column_N <- function(row) {
  cumsum_row <- unlist(cumsum(replace(row, is.na(row), 0))) 
  if (cumsum_row[1] > time_point) {  
    return(0)  
  }
  max(which(cumsum_row <= time_point), na.rm = TRUE)  
}
max_column <- apply(time_data, 1, find_column_N) 


#### Input Matrix Generation 
K <- ncol(time_data)
N <- nrow(time_data)
mat <- matrix(0, nrow = N, ncol= K)
for (i in 1:K){
  # Use response answer for fast responses
  out <- which(time_data[,i] < median(time_data[, i], na.rm = T)/2)
  mat[out,i] <- resp_data[out,i]
  # Convert missing values to zero
  mat[which(mat[,i]==79), i] <- 0
  miss <- which(is.na(resp_data[, i]))
  mat[miss, i] <- 0 
}

# Customized QUBIC function : exclude discretization 
qubic_no_disc <- function (x, r = 1L, q = 0.06, c = 0.95, o = 100, f = 1, 
                           k = max(ncol(x)%/%20, 2), type = "default", P = FALSE, 
                           C = FALSE, verbose = TRUE, 
                           weight = NULL, seedbicluster = NULL) 
{x_d <- x
return(qubiclust_d(x_d, c, o, f, k, type, P, C, verbose, 
                   weight, seedbicluster))
}

#### Cheating detection using QUBIC
res <- tryCatch(
  qubic_no_disc(mat, P = T,  c = 1, k = 4, verbose = F), error = function(e) {NULL}
)
if(length(res) > 0){
  biclusters <- list()
  for (h in 1:res@Number) { # Collect biclusters
    bicluster <- list(
      rows = which(res@RowxNumber[, h] == T),
      cols = which(res@NumberxCol[h, ] == T)
    )
    # Filter biclusters with > 80% correct responses
    if (length(which(apply(mat[bicluster$rows, bicluster$cols], 2, Mode) == 
      key_data[bicluster$cols])) < length(bicluster$cols) * 0.8)  
      
    { biclusters[[h]] <- NULL
    } else {
      biclusters[[h]] <- bicluster # Collect filtered biclusters
    }
  }
  # Calculate p-values for biclusters and filter by p < 0.01
  if (length(biclusters) > 0) {
    p_mat <- matrix(list(), nrow = length(biclusters), ncol = 4)
    
    # Compute ability-based strata for examinees 
    score_all <- table(which(resp_data==key_data, arr.ind=T)[,1])
    breaks <- unique(quantile(score_all, seq(0, 1, 0.1)))
    strata <- cut(score_all, breaks = breaks, include.lowest = T)
    
    for (z in 1:length(biclusters)) {
      if (is.null(biclusters[[z]]) || 
          length(biclusters[[z]]$rows) == 0 || 
          length(biclusters[[z]]$cols) == 0) {
        p_mat[z,] <- NA
      } else {
        rowss <- biclusters[[z]]$rows
        colss <- biclusters[[z]]$cols
        p_mat[[z, 2]] <- length(colss)
        p_mat[[z, 3]]<- length(rowss)
        
        ori_pat <- apply(mat[rowss, colss], 2, Mode) # Pattern within a bicluster
        pattern <- rep(NA, ncol(mat))  
        pattern[colss] <- ori_pat

        column_probs <- sapply(colss, function(col) {
          length(which(cur_resp[, col] == pattern[col]))/(length(which(max_column >= col)))/2
        })
        pattern_prob <- prod(c(column_probs), na.rm = T) # Probability of the same pattern
        observed_size <- length(rowss)
        
        # Calculate the p-value using Binomial distribution
        p_value <- pbinom(observed_size - 1, size = N, prob = pattern_prob, lower.tail = F)
        p_mat[[z,1]] <- p_value
        
        # Ability-based filtering
        strata_prop <- matrix(NA, nrow=length(unique(strata)), ncol=2)
        strata_prop[,2] <- unique(strata)
        for(st in unique(strata)){
          matching_rows <- which(apply(resp_data[which(strata == st), colss], 1, function(x) identical(as.numeric(as.vector(x)), ori_pat)))
          strata_prop[which(unique(strata)==st),1] <- length(matching_rows)/ length(which(strata==st))
        }
        rare_cases <- intersect(strata[rowss], unique(strata)[which(strata_prop[,1] < 0.1)])
        p_mat[[z, 4]] <- rowss[which(strata[rowss] %in% rare_cases)] # Flag examinees meeting ability-based filtering criteria
        if (length(p_mat[[z, 4]]) < 2) {
          p_mat[[z, 4]] <- NA
        }
      }
    }
    person <- integer()
    signi <- which((p_mat[,1]) < 0.05) # Filter biclusters by p-value < 0.05
    
    # Flag examinees from biclusters with p-value < 0.05
    for (c in signi) {
      person_current <- biclusters[[c]]$rows
      person <- union(person, person_current)
    }
  }
}
flagged_examinees <- unique(c(person, unlist(p_mat[, 4]))) # Flag all possible cheaters  
flagged_examinees
