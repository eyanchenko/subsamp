### Comparing graph sub-sampling methods for CP divide-and-conquer
### Eric Yanchenko
### February 6, 2023
library(pdfCluster)
library(Rcpp)
library(doParallel)
library(dplyr)
library(igraph)
library(LaplacesDemon)
library(pROC)
library(RColorBrewer)
library(reticulate)
library(ggplot2)
library(ggpubr)

# Generate SBM networks
generateA <- function(n, p11, p12, p22, prop=0.50){
  m = as.integer(round(n^2*prop^2, digits=0))
  A11 = matrix(rbinom(m, 1, p11), ncol=as.integer(round(n*prop)))
  A11[lower.tri(A11)] <- 0
  diag(A11) <- 0
  A11 = A11 + t(A11)
  
  m = as.integer(round(n^2*(1-prop)^2))
  A22 = matrix(rbinom(m, 1, p22), ncol=as.integer(round(n*(1-prop))), nrow = as.integer(round(n*(1-prop))))
  A22[lower.tri(A22)] <- 0
  diag(A22) <- 0
  A22 = A22 + t(A22)
  
  m = as.integer(round(n^2*prop*(1-prop)))
  A12 <- matrix(rbinom(m, 1, p12), ncol=as.integer(round(n*(1-prop))), nrow = as.integer(round(n*prop)))
  A21 <- t(A12)
  
  A = cbind(rbind(A11, A21), rbind(A12, A22))
  return(A)
}


## NODE
# Random node, uniform
node_unif <- function(g, no_nodes){
  n = vcount(g)
  idx = sample(1:n, no_nodes, replace = FALSE)
  return(idx)
}

# Random node, degree
node_deg <- function(g, no_nodes, probs){
  
  idx = rcat(no_nodes, probs)
  idx = unique(idx)
  m = length(unique(idx))
  
  # Ensure that we have unique nodes
  while( m < no_nodes){
    probs[idx] <- 0
    idx = c(idx, rcat(no_nodes-m, probs))
    idx = unique(idx)
    m = length(idx)
  }
  
  return(idx[1:no_nodes])
  
}

## EDGE
# Random edge, uniform
edge_unif <- function(g, no_nodes){
  gg <- as_edgelist(g)
  m  <- dim(gg)[1]
  
  # Nodes that are sampled
  idx <- c()
  
  # Keep sampling edges until we have no_node nodes sampled
  while( length(idx) < no_nodes ){
    id  = sample(1:m, 1) # Probability of sampling the same edge twice is very low
    # Even if it happens, no new nodes will be added to set idx
    # since these nodes were already included
    idx = c(idx, gg[id, ])
    idx = unique(idx)
    
  }
  
  return(idx)
}

## EXPLORATION
# Breadth-first search
bfs_nodes <- function(g, no_nodes){
  # Randomly select root
  n = vcount(g)
  strt = sample(1:n,1)
  # Keep first no_nodes that are traversed
  idx = as.numeric(bfs(g, strt)$order[1:no_nodes])
  
  return(idx)
}

# Depth-first search
dfs_nodes <- function(g, no_nodes){
  # Randomly select root
  n = vcount(g)
  strt = sample(1:n,1)
  # Keep first no_nodes that are traversed
  idx = as.numeric(dfs(g, strt)$order[1:no_nodes])
  
  return(idx)
}

# Random node-neighbor
random_node_neigh <- function(g, no_nodes){
  gg <- as_edgelist(g)
  n  <- vcount(g)
  
  # Nodes that are sampled
  idx <- c()
  
  # Keep sampling until we have no_node nodes sampled
  while( length(idx) < no_nodes ){
    # Random seed node
    id = sample(1:n, 1)
    
    # Find all neighbors of id
    neighs = unique(as.vector(gg[gg[,1]== id | gg[,2]==id, ]))
    idx = c(idx, neighs)
    idx = unique(idx)
    
  }
  
  if(length(idx) > no_nodes){
    idx = sample(idx, no_nodes)
  }
  
  return(idx)
}

# Random walk
rand_walk <- function(g, no_nodes){
  # Random starting location
  n = vcount(g)
  strt = sample(1:n, 1)
  
  idx = as.numeric(random_walk(graph=g, start=strt, steps=no_nodes-1))
  
  idx = unique(idx)
  m = length(idx)
  
  while(m < no_nodes){
    strt = sample(1:n, 1)
    idx = c(idx, as.numeric(random_walk(graph=g, start=strt, steps=no_nodes-m-1)))
    idx = unique(idx)
    m = length(idx)
  }
  
  
  return(idx)
}



### Borgatti test statistic

### Function to find the core-periphery assignments which maximize the Borgatti and Everett (2000) metric
### Written in Rcpp

cppFunction('IntegerVector borgattiCpp(IntegerMatrix A){

  double n = A.nrow();
  double m = n*(n-1)/2;
  int k = 0;
  
  
  // Initialize vector 
  
  IntegerVector C(n);
  
  
  for(int i = 0; i < n; i++){
    
    if(rand()%10 > 5){
      C(i) = 1;
      k++;
    }else{
      C(i) = 0;
    }
     
  }
  
  int kk = k;
  
  
  double sum_A = 0; 
  double sum_CP = 0.5*k*(k-1)+(n-k)*k; 
  double sum_ACP = 0;
  
  for (int i = 1; i < n; i++){
      for(int j = 0; j < i; j++){
        sum_A   += A(i,j);
        sum_ACP += A(i,j) * (C(i) + C(j) - C(i)*C(j)); 
      }
    }
  
  double obj_last = (sum_ACP - sum_A / m * sum_CP) / 
    (sqrt(sum_A - sum_A/m*sum_A) * sqrt(sum_CP - sum_CP/m*sum_CP) );
  
                   
  int ind = 1;
  int n_iter = 0;
  int max_iter = 100;
  
  IntegerVector Ctest = clone(C);
  
  IntegerVector v(n);
  for(int i=0; i<n;i++){v(i) = i;}
  
  double obj_new = 0;
  
  while(ind > 0 & n_iter < max_iter){
  
  ind = 0;
  /* Randomly shuffle node order */
  int N = n;
  for(int a=0;a<n;a++){
    int b = round(rand()%N);
    int c = v(a); v(a)=v(b); v(b)=c;
  }
  
  
  
   for(int i = 0; i < n; i++){
    
     Ctest = clone(C);
     
     Ctest(v(i)) = 1 - Ctest(v(i));
     
    /* Update size of core for test vector, kk */
    
     if(Ctest(v(i))==0){
       kk = k - 1;
     }else{
       kk = k + 1;
     }

     double sum_CP_test = 0.5*kk*(kk-1) + (n-kk)*kk; 
     double sum_ACP_test = sum_ACP;
    
    /* Update value of objective function */
    
     for(int j = 0; j < n; j++){
        sum_ACP_test += A(v(i),j) * (Ctest(v(i)) + Ctest(j) - Ctest(v(i))*Ctest(j)); 
        sum_ACP_test -= A(v(i),j) * (C(v(i)) + C(j) - C(v(i))*C(j)); 
     }

    
     obj_new = (sum_ACP_test - sum_A / m * sum_CP_test) / 
      (sqrt(sum_A - sum_A/m*sum_A) * sqrt(sum_CP_test - sum_CP_test/m*sum_CP_test) );

    
     if(obj_new > obj_last){
       C(v(i)) = Ctest(v(i));
       sum_ACP = sum_ACP_test;
       ind++;
       obj_last = obj_new;
       k = kk;
     }
    
    
    }
    
    n_iter++;
    
  }
  
  return( C );
}')



#' @title Core-Periphery sub-sampling algorithm (using mclapply)
#' @description This function finds the CP labels of
#' a large network using a sub-sampling heuristic
#' @param A adjacency matrix
#' @param q proportion of network to sub-sample
#' @param B number of iterations to run
#' @return Proportion of sub-samples assigned to the core
#' @export

ssCP <- function(A, q=0.10, B=1000, sampler=node_unif, no_cores = 9, ...){
  start = proc.time()
  n = dim(A)[1]
  g <- graph_from_adjacency_matrix(A)
  
  apply_fun <- function(i){
    # sub-sample nodes
    idx <- sampler(g, round(n*q), ...)
    AA   <- A[idx, idx]
    
    # re-sample until non-empty matrix is returned
    while(sum(AA) == 0){
      idx <- sampler(g, round(n*q), ...)
      AA   <- A[idx, idx]
    }
    
    # Run Borgatti and Everett algorithm
    start = proc.time()
    C <- as.logical( borgattiCpp(AA) )
    end = proc.time()
    
    # return indices of nodes which were assigned to the core
    return(list(core=idx[C], time=as.numeric((end-start)[3])))
  }
  
  
  out <- mclapply(1:B, apply_fun, mc.cores = no_cores)

  
  core = out[[1]]$core
  sub_time = out[[1]]$time/B
  
  for(i in 2:B){
    core = c(core, out[[i]]$core)
    sub_time = sub_time + out[[i]]$time/B
  }
  
  # Count number of times each node is assigned to the core
  props = numeric(n)
  
  for(i in 1:n){
    props[i] = sum(core==i)
  }
  
  end = proc.time()
  total_time = as.numeric((end-start)[3])

  # return proportions
  return(list(props = (props / B), total_time = total_time, sub_time = sub_time))
  
}


# Set colors
myColors <- brewer.pal(8,"Set1")
myColors[6] <- myColors[8]
myColors[8] <- "black"

######### SIMULATION EXPERIMENTS
mc_cores = 9
n.iter = 100

########## Fixed n, increasing rho
n=5000
p22 = 0.001
p11.seq <- seq(0.002, 0.02, 0.001)
prop = 0.01
Ctrue = numeric(n)
Ctrue[1:(n*prop)] <- 1
B = 1000

df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(p11.seq)),
             p11  = 0, 
             AUC = 0, 
             timeT = 0, # total run time
             timeS = -1, # sub-sample run time
             acc = -1)  # accuracy (only for full algorithm)

idx = 1
for(j in 1:length(p11.seq)){
  p11 <- p11.seq[j]
  p12 <- p11/2
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, p22, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- p11
    
    
    out <- ssCP(A, q=100/n, sampler=node_unif, B = B, no_cores = mc_cores)
    df[idx,4] <- auc(Ctrue, out$props)[1]
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time

    deg = colSums(A)
    out <- ssCP(A, q=100/n, sampler=node_deg, probs=deg, B = B, no_cores = mc_cores)
    df[idx+1,4] <- auc(Ctrue, out$props)[1]
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time

    out <- ssCP(A, q=100/n, sampler=edge_unif, B = B, no_cores = mc_cores)
    df[idx+2,4] <- auc(Ctrue, out$props)[1]
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=bfs_nodes, B = B, no_cores = mc_cores)
    df[idx+3,4] <- auc(Ctrue, out$props)[1]
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time

    out <- ssCP(A, q=100/n, sampler=dfs_nodes, B = B, no_cores = mc_cores)
    df[idx+4,4] <- auc(Ctrue, out$props)[1]
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time

    out <- ssCP(A, q=100/n, sampler=random_node_neigh, B = B, no_cores = mc_cores)
    df[idx+5,4] <- auc(Ctrue, out$props)[1]
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time

    out <- ssCP(A, q=100/n, sampler=rand_walk, B = B, no_cores = mc_cores)
    df[idx+6,4] <- auc(Ctrue, out$props)[1]
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    df[idx+7,5] <- as.numeric(system.time(props <- borgattiCpp(A))[3])
    df[idx+7,4] <- auc(Ctrue, props)[1] # 0.85
    df[idx+7,7] <- mean(Ctrue==props) # 0.70

    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_n5000_vary_p.RData")
  }
  print(p11)
}


load("sims_n5000_vary_p.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))

df_plot <- df %>% group_by(method, p11) %>% 
  summarize(AUC=mean(AUC), timeT=mean(timeT), timeS=mean(timeS), acc=mean(acc))

pM1 <- ggplot(df_plot, aes(x=p11, y=AUC, color=method, linetype=method))+
  ylim(0.49, 1.00)+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=0.0025, y=0.99, label="(7)", size=6)+
  xlab(bquote(p[11]))
pM1

pT1 <- ggplot(df_plot, aes(x=p11, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=0.0025, y=10, label="(7)", size=6)+
  xlab(bquote(p[11]))+
  ylab("Total time (s)")
pT1

df_plot <- df_plot[df_plot$method!="FULL", ]
pS1 <- ggplot(df_plot, aes(x=p11, y=timeS, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=0.0025, y=0.0050, label="(7)", size=6)+
  xlab(bquote(p[11]))+
  ylab("Sub-graph detection time (s)")
pS1



########## Fixed rho, increasing n (fixed B)
n.seq = seq(1000, 10000, 1000)
p22 = 0.001
p11 = 0.004
p12 = 0.002
prop = 0.01
B = 1000

df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(n.seq)),
             n  = 0, 
             AUC = 0, 
             timeT = 0, # total run time
             timeS = -1, # sub-sample run time
             acc = -1)  # accuracy (only for full algorithm)

idx = 1
for(n in n.seq){
  Ctrue = numeric(n)
  Ctrue[1:(n*prop)] <- 1
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, p22, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- n
    
    
    out <- ssCP(A, q=100/n, sampler=node_unif, B = B, no_cores = mc_cores)
    df[idx,4] <- auc(Ctrue, out$props)[1]
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- ssCP(A, q=100/n, sampler=node_deg, probs=deg, B = B, no_cores = mc_cores)
    df[idx+1,4] <- auc(Ctrue, out$props)[1]
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=edge_unif, B = B, no_cores = mc_cores)
    df[idx+2,4] <- auc(Ctrue, out$props)[1]
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=bfs_nodes, B = B, no_cores = mc_cores)
    df[idx+3,4] <- auc(Ctrue, out$props)[1]
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=dfs_nodes, B = B, no_cores = mc_cores)
    df[idx+4,4] <- auc(Ctrue, out$props)[1]
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=random_node_neigh, B = B, no_cores = mc_cores)
    df[idx+5,4] <- auc(Ctrue, out$props)[1]
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=rand_walk, B = B, no_cores = mc_cores)
    df[idx+6,4] <- auc(Ctrue, out$props)[1]
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    df[idx+7,5] <- as.numeric(system.time(props <- borgattiCpp(A))[3])
    df[idx+7,4] <- auc(Ctrue, props)[1]
    df[idx+7,7] <- mean(Ctrue==props)
    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_n5000_vary_n_fixedB.RData")
  }
  print(n)
}


load("sims_n5000_vary_n_fixedB.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))
# Only keep n<=8000
df <- df[df$n<=8000, ]

df_plot <- df %>% group_by(method, n) %>% 
  summarize(AUC=mean(AUC), timeT=mean(timeT), timeS=mean(timeS), acc=mean(acc))

pM2 <- ggplot(df_plot, aes(x=n, y=AUC, color=method, linetype=method))+
  ylim(0.49, 1.00)+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1200, y=0.99, label="(8)", size=6)
pM2

pT2 <- ggplot(df_plot, aes(x=n, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1200, y=35, label="(8)", size=6)+
  ylab("Total time (s)")
pT2

df_plot <- df_plot[df_plot$method!="FULL", ]
pS2 <- ggplot(df_plot, aes(x=n, y=timeS, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1200, y=0.0045, label="(8)", size=6)+
  ylab("Sub-graph detection time (s)")
pS2


########## Fixed rho, increasing n (increasing B)
n.seq = seq(1000, 10000, 1000)
p22 = 0.001
p11 = 0.004
p12 = 0.002
prop = 0.01


df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(n.seq)),
             n  = 0, 
             AUC = 0, 
             timeT = 0, # total run time
             timeS = -1, # sub-sample run time
             acc = -1)  # accuracy (only for full algorithm)

idx = 1
for(n in n.seq){
  B = n / 2
  Ctrue = numeric(n)
  Ctrue[1:(n*prop)] <- 1
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, p22, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- n
    
    
    out <- ssCP(A, q=100/n, sampler=node_unif, B = B, no_cores = mc_cores)
    df[idx,4] <- auc(Ctrue, out$props)[1]
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- ssCP(A, q=100/n, sampler=node_deg, probs=deg, B = B, no_cores = mc_cores)
    df[idx+1,4] <- auc(Ctrue, out$props)[1]
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=edge_unif, B = B, no_cores = mc_cores)
    df[idx+2,4] <- auc(Ctrue, out$props)[1]
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time

    out <- ssCP(A, q=100/n, sampler=bfs_nodes, B = B, no_cores = mc_cores)
    df[idx+3,4] <- auc(Ctrue, out$props)[1]
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=dfs_nodes, B = B, no_cores = mc_cores)
    df[idx+4,4] <- auc(Ctrue, out$props)[1]
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=random_node_neigh, B = B, no_cores = mc_cores)
    df[idx+5,4] <- auc(Ctrue, out$props)[1]
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=rand_walk, B = B, no_cores = mc_cores)
    df[idx+6,4] <- auc(Ctrue, out$props)[1]
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    # df[idx+7,5] <- as.numeric(system.time(props <- borgattiCpp(A))[3])
    # df[idx+7,4] <- auc(Ctrue, props)[1]
    # df[idx+7,7] <- mean(Ctrue==props)
    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_n5000_vary_n_varyB.RData")
  }
  print(n)
}


load("sims_n5000_vary_n_varyB.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))
df <- df[df$n<=8000, ]

df_plot <- df %>% group_by(method, n) %>% 
  summarize(AUC=mean(AUC), timeT=mean(timeT), timeS=mean(timeS), acc=mean(acc))

# Load results from Setting 2 for FULL
load("sims_n5000_vary_n_fixedB.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))
df <- df[df$n<=8000, ]

df_plot2 <- df %>% group_by(method, n) %>% 
  summarize(AUC=mean(AUC), timeT=mean(timeT), timeS=mean(timeS), acc=mean(acc))

df_plot[df_plot$method=="FULL", ] = df_plot2[df_plot2$method=="FULL", ]


pM3 <- ggplot(df_plot, aes(x=n, y=AUC, color=method, linetype=method))+
  ylim(0.49, 1.00)+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1200, y=0.99, label="(9)", size=6)
pM3

pT3 <- ggplot(df_plot, aes(x=n, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1200, y=35, label="(9)", size=6)+
  ylab("Total time (s)")
pT3

df_plot <- df_plot[df_plot$method!="FULL", ]
pS3 <- ggplot(df_plot, aes(x=n, y=timeS, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1200, y=0.0035, label="(9)", size=6)+
  ylab("Sub-graph detection time (s)")
pS3



########## Fixed n, rho, varying core size

n=5000
p22 = 0.001
p11 = 0.004
p12 = 0.002
prop.seq = c(0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.20, 0.30)

B = 1000

df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(prop.seq)),
             prop  = 0, 
             AUC = 0, 
             timeT = 0, # total run time
             timeS = -1, # sub-sample run time
             acc = -1)  # accuracy (only for full algorithm)

idx = 1
for(prop in prop.seq){
  
  Ctrue = numeric(n)
  Ctrue[1:(n*prop)] <- 1
  
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, p22, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- prop
    
    
    out <- ssCP(A, q=100/n, sampler=node_unif, B = B, no_cores = mc_cores)
    df[idx,4] <- auc(Ctrue, out$props)[1]
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- ssCP(A, q=100/n, sampler=node_deg, probs=deg, B = B, no_cores = mc_cores)
    df[idx+1,4] <- auc(Ctrue, out$props)[1]
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=edge_unif, B = B, no_cores = mc_cores)
    df[idx+2,4] <- auc(Ctrue, out$props)[1]
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=bfs_nodes, B = B, no_cores = mc_cores)
    df[idx+3,4] <- auc(Ctrue, out$props)[1]
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=dfs_nodes, B = B, no_cores = mc_cores)
    df[idx+4,4] <- auc(Ctrue, out$props)[1]
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=random_node_neigh, B = B, no_cores = mc_cores)
    df[idx+5,4] <- auc(Ctrue, out$props)[1]
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=rand_walk, B = B, no_cores = mc_cores)
    df[idx+6,4] <- auc(Ctrue, out$props)[1]
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    df[idx+7,5] <- as.numeric(system.time(props <- borgattiCpp(A))[3])
    df[idx+7,4] <- auc(Ctrue, props)[1]
    df[idx+7,7] <- mean(Ctrue==props)
    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_n5000_vary_prop.RData")
  }
  print(prop)
}


load("sims_n5000_vary_prop.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))


df_plot <- df %>% group_by(method, prop) %>% summarize(AUC=mean(AUC), timeT=mean(timeT), timeS=mean(timeS))

pM4 <- ggplot(df_plot, aes(x=prop, y=AUC, color=method, linetype=method))+
  ylim(0.49, 1.00)+
  geom_line(linewidth=1.5)+
  xlab(expression(alpha))+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=0.015, y=0.99, label="(10)", size=6)
pM4

pT4 <- ggplot(df_plot, aes(x=prop, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  xlab(expression(alpha))+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=0.015, y=19, label="(10)", size=6)+
  ylab("Total time (s)")
pT4

df_plot <- df_plot[df_plot$method!="FULL", ]
pS4 <- ggplot(df_plot, aes(x=prop, y=timeS, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  xlab(expression(alpha))+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=0.015, y=0.004, label="(10)", size=6)+
  ylab("Sub-graph detection time (s)")
pS4


########## Fix all but B

n = 5000
p22 = 0.001
p11 = 0.004
p12 = 0.002
prop = 0.01
B.seq = c(100, 250, 500, 1000, 2500, 5000, 10000)

Ctrue = numeric(n)
Ctrue[1:(n*prop)] <- 1

df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(B.seq)),
             B  = 0, 
             AUC = 0, 
             timeT = 0, # total run time
             timeS = -1, # sub-sample run time
             acc = -1)  # accuracy (only for full algorithm)

idx = 1
for(B in B.seq){

  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, p22, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- B
    
    
    out <- ssCP(A, q=100/n, sampler=node_unif, B = B, no_cores = mc_cores)
    df[idx,4] <- auc(Ctrue, out$props)[1]
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- ssCP(A, q=100/n, sampler=node_deg, probs=deg, B = B, no_cores = mc_cores)
    df[idx+1,4] <- auc(Ctrue, out$props)[1]
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=edge_unif, B = B, no_cores = mc_cores)
    df[idx+2,4] <- auc(Ctrue, out$props)[1]
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=bfs_nodes, B = B, no_cores = mc_cores)
    df[idx+3,4] <- auc(Ctrue, out$props)[1]
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=dfs_nodes, B = B, no_cores = mc_cores)
    df[idx+4,4] <- auc(Ctrue, out$props)[1]
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=random_node_neigh, B = B, no_cores = mc_cores)
    df[idx+5,4] <- auc(Ctrue, out$props)[1]
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time
    
    out <- ssCP(A, q=100/n, sampler=rand_walk, B = B, no_cores = mc_cores)
    df[idx+6,4] <- auc(Ctrue, out$props)[1]
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    # df[idx+7,5] <- as.numeric(system.time(props <- borgattiCpp(A))[3])
    # df[idx+7,4] <- auc(Ctrue, props)[1]
    # df[idx+7,7] <- mean(Ctrue==props)
    
    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_n5000_vary_B.RData")
  }
  print(B)
}


load("sims_n5000_vary_B.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))
df_plot <- df %>% group_by(method, B) %>% 
  summarize(AUC=mean(AUC), timeT=mean(timeT), timeS=mean(timeS), acc=mean(acc))


# Load results from Setting 1 for FULL
load("sims_n5000_vary_p.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))

df_plot2 <- df %>% group_by(method, p11) %>% 
  summarize(AUC=mean(AUC), timeT=mean(timeT), timeS=mean(timeS), acc=mean(acc))
df_plot2 <- df_plot2[as.logical((df_plot2$method=="FULL")*(df_plot2$p11==0.004)), ]

df_plot$AUC[df_plot$method=="FULL"] <- df_plot2$AUC
df_plot$timeT[df_plot$method=="FULL"] <- df_plot2$timeT

pM5 <- ggplot(df_plot, aes(x=B, y=AUC, color=method, linetype=method))+
  ylim(0.49, 1.00)+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=800, y=0.99, label="(11)", size=6)
pM5

pT5 <- ggplot(df_plot, aes(x=B, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=800, y=100, label="(11)", size=6)+
  ylab("Total time (s)")
pT5

df_plot <- df_plot[df_plot$method!="FULL", ]
pS5 <- ggplot(df_plot, aes(x=B, y=timeS, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=800, y=0.0045, label="(11)", size=6)+
  ylab("Sub-graph detection time (s)")
pS5



########## Fix all but q

n = 5000
p22 = 0.001
p11 = 0.004
p12 = 0.002
prop = 0.01
B = 1000

q.seq = c(50, 100, 250, 500) / n

Ctrue = numeric(n)
Ctrue[1:(n*prop)] <- 1

df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(q.seq)),
             q  = 0, 
             AUC = 0, 
             timeT = 0, # total run time
             timeS = -1, # sub-sample run time
             acc = -1)  # accuracy (only for full algorithm)

idx = 1
for(q in q.seq){
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, p22, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- q
    
    out <- ssCP(A, q=q, sampler=node_unif, B = B, no_cores = mc_cores)
    df[idx,4] <- auc(Ctrue, out$props)[1]
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- ssCP(A, q=q, sampler=node_deg, probs=deg, B = B, no_cores = mc_cores)
    df[idx+1,4] <- auc(Ctrue, out$props)[1]
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- ssCP(A, q=q, sampler=edge_unif, B = B, no_cores = mc_cores)
    df[idx+2,4] <- auc(Ctrue, out$props)[1]
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- ssCP(A, q=q, sampler=bfs_nodes, B = B, no_cores = mc_cores)
    df[idx+3,4] <- auc(Ctrue, out$props)[1]
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- ssCP(A, q=q, sampler=dfs_nodes, B = B, no_cores = mc_cores)
    df[idx+4,4] <- auc(Ctrue, out$props)[1]
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- ssCP(A, q=q, sampler=random_node_neigh, B = B, no_cores = mc_cores)
    df[idx+5,4] <- auc(Ctrue, out$props)[1]
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time
    
    out <- ssCP(A, q=q, sampler=rand_walk, B = B, no_cores = mc_cores)
    df[idx+6,4] <- auc(Ctrue, out$props)[1]
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    # df[idx+7,5] <- as.numeric(system.time(props <- borgattiCpp(A))[3])
    # df[idx+7,4] <- auc(Ctrue, props)[1]
    # df[idx+7,7] <- mean(Ctrue==props)
    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_n5000_vary_q.RData")
  }
  print(q)
}


load("sims_n5000_vary_q.RData")

df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))
df_plot <- df %>% group_by(method, q) %>% 
  summarize(AUC=mean(AUC), timeT=mean(timeT), timeS=mean(timeS), acc=mean(acc))


# Load results from Setting 1 for FULL
load("sims_n5000_vary_p.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))

df_plot2 <- df %>% group_by(method, p11) %>% 
  summarize(AUC=mean(AUC), timeT=mean(timeT), timeS=mean(timeS), acc=mean(acc))
df_plot2 <- df_plot2[as.logical((df_plot2$method=="FULL")*(df_plot2$p11==0.004)), ]

df_plot$AUC[df_plot$method=="FULL"] <- df_plot2$AUC
df_plot$timeT[df_plot$method=="FULL"] <- df_plot2$timeT


pM6 <- ggplot(df_plot, aes(x=q*n, y=AUC, color=method, linetype=method))+
  ylim(0.49, 1.00)+
  xlab("qn")+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=70, y=0.99, label="(12)", size=6)
pM6

pT6 <- ggplot(df_plot, aes(x=q*n, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  xlab("qn")+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=70, y=27, label="(12)", size=6)+
  ylab("Total time (s)")
pT6

df_plot <- df_plot[df_plot$method!="FULL", ]
pS6 <- ggplot(df_plot, aes(x=q*n, y=timeS, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  xlab("qn")+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=70, y=0.09, label="(12)", size=6)+
  ylab("Sub-graph detection time (s)")
pS6



ggarrange(pM1, pM2, pM3, pM4, pM5, pM6, ncol=2, nrow=3, common.legend = T, legend="bottom")

ggsave(file="cp_sims.pdf",
       height = 8,
       width = 8.5,
       units="in")


ggarrange(pT1, pT2, pT3, pT4, pT5, pT6, ncol=2, nrow=3, common.legend = T, legend="bottom")

ggsave(file="~cp_sims_timeT.pdf",
       height = 8,
       width = 8.5,
       units="in")

ggarrange(pS1, pS2, pS3, pS4, pS5, pS6, ncol=2, nrow=3, common.legend = T, legend="bottom")

ggsave(file="cp_sims_timeS.pdf",
       height = 8,
       width = 8.5,
       units="in")






