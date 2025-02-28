library(pdfCluster)
library(doParallel)
library(dplyr)
library(igraph)
library(LaplacesDemon)
library(pROC)
library(reticulate)
library(ggplot2)
library(RandPro)
library(RColorBrewer)
library(ggpubr)

# Generate SBM networks
generateA <- function(n, pin, pout, props){
  
  K <- length(props)
  
  P <- matrix(0, nrow=n, ncol=n)
  
  idx <- 1
  # Block probabilities
  for(k in 1:K){
    intra <- idx:(idx + n*props[k]-1)
    inter <- (1:n)[-intra]
    P[intra, intra] <- pin  # intra-community probs.
    P[intra, inter] <- pout # inter-community probs.
    P[inter, intra] <- pout
    
    idx =  idx + n*props[k]
  }
  
  A <- matrix(rbinom(n^2, 1, P), ncol=n, nrow=n)
  
  A[lower.tri(A, diag=TRUE)] <- 0
  A = A + base::t(A)
  
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

#' @title PACE algorithm
#' @description This function finds the community labels
#' @param A adjacency matrix
#' @param q proportion of network to sub-sample
#' @param no_comms pre-specified number of communities
#' @param B number of iterations to run
#' @return Proportion of sub-samples assigned to the core
#' @export

pace <- function(A, q=0.10, B=1000, no_comms=2, clustering = cluster_fast_greedy, sampler=node_unif, no_cores = 9, ...){
  start1 = proc.time()
  n = dim(A)[1]
  g <- as.undirected(graph_from_adjacency_matrix(A))
  
  apply_fun <- function(i){
    # sub-sample nodes
    idx <- sampler(g, round(n*q), ...)
    idx <- sort(idx)
    AA   <- A[idx, idx]
    
    # re-sample until non-empty matrix is returned
    while(sum(AA) == 0){
      idx <- sampler(g, round(n*q), ...)
      idx <- sort(idx)
      AA   <- A[idx, idx]
    }
    
    # Run clustering algorithm
    gg <- as.undirected(graph_from_adjacency_matrix(AA))
    start = proc.time()
    clust <- clustering(gg)
    end = proc.time()

    # First column is the nodes that were sampled
    # Second column is their group membership
    ret_mat = matrix(0, nrow = length(idx), ncol=2)
    ret_mat[,1] <- idx
    # ret_mat[,2] <- clust
    ret_mat[,2] <- clust$membership
    
    return(list(ret_mat=ret_mat, time=as.numeric((end-start)[3])))
  }
  
  out <- mclapply(1:B, apply_fun, mc.cores = no_cores)
  
  Chat = matrix(0, nrow = n, ncol = n)
  NN   = matrix(0, nrow = n, ncol = n)
  
  sub_time = 0
  
  for(b in 1:B){
    K = max(out[[b]]$ret_mat[,2])
    sub_time = sub_time + out[[b]]$time/B
    
    for(k in 1:K){
      idx = out[[b]]$ret_mat[,1][out[[b]]$ret_mat[,2] == k]
      Chat[idx, idx] = Chat[idx, idx] + 1 
    }
    
    idx = out[[b]]$ret_mat[,1]
    NN[idx, idx] = NN[idx, idx] + 1 
  }
  # Set tuning parameter tau to be 40th%tile of N
  tau = quantile(NN[upper.triangle(NN)], 0.4)
  
  # Compute final Chat as number of times assigned to same cluster divided by number of times in same sub-sampler
  Chat = Chat / NN 
  Chat[NN < tau] <- 0
  Chat[NN==0] <- 0
  
  # # Project Chat to a lower dimension using random projection matrix (Gaussian)
  # # s = 10*floor(log(n))
  # # R <- form_matrix(n, s, JLT = F)
  # # Cproj = Chat %*% R / sqrt(s)
  # 
  # # Lastly, must cluster on Chat to get final membership vector
  # # labels = kmeans(Cproj, no_comms)$cluster
  # 
  
  labels = kmeans(Chat, no_comms)$cluster
  
  end1 = proc.time()
  total_time = as.numeric((end1-start1)[3])
  
  return(list(labels=labels, total_time=total_time, sub_time=sub_time))
}

# Set colors
myColors <- brewer.pal(8,"Set1")
myColors[6] <- myColors[8]
myColors[8] <- "black"


######### SIMULATION EXPERIMENTS
mc_cores = 9
n.iter = 100
########## Fixed n, increasing p11 - p12

n=5000
p12 = 0.01
p11.seq = seq(0.02, 0.10, 0.005)
prop = c(0.75, 0.25)
Ctrue = c(rep(1, n*prop[1]), rep(2, n*prop[2]))
B = 1000



df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(p11.seq)),
             p11  = 0, 
             ARI = 0, 
             timeT = 0, # total run time
             timeS = -1) # sub-sample run time

idx = 1
for(j in 1:length(p11.seq)){
  p11 <- p11.seq[j]
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- p11
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_unif)
    df[idx,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_deg, probs=deg)
    df[idx+1,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=edge_unif)
    df[idx+2,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=bfs_nodes)
    df[idx+3,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=dfs_nodes)
    df[idx+4,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=random_node_neigh)
    df[idx+5,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time

    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=rand_walk)
    df[idx+6,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    G <- as.undirected(graph_from_adjacency_matrix(A))
    df[idx+7,5] <- as.numeric(system.time(out <- cluster_fast_greedy(G)$membership)[3])
    df[idx+7,4] <- adj.rand.index(Ctrue, out)
    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_pace_n5000_vary_p.RData")
  }
  print(p11)
}



load("sims_pace_n5000_vary_p.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))

df_plot <- df %>% group_by(method, p11) %>% 
  summarize(ARI=mean(ARI), timeT=mean(timeT), timeS=mean(timeS))

pM1 <- ggplot(df_plot, aes(x=p11, y=ARI, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=0.025, y=0.99, label="(1)", size=6)+
  xlab(bquote(p[11]))
pM1

pT1 <- ggplot(df_plot, aes(x=p11, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=0.025, y=30, label="(1)", size=6)+
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
  annotate(geom="text", x=0.025, y=0.013, label="(1)", size=6)+
  xlab(bquote(p[11]))+
  ylab("Sub-graph detection time (s)")
pS1


########## Fixed p11 - p12, increasing n (fixed B)

n.seq = seq(1000, 8000, 1000)
p12 = 0.01
p11 = 0.04
prop = c(0.75, 0.25)
Ctrue = c(rep(1, n*prop[1]), rep(2, n*prop[2]))

B = 1000

df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(n.seq)),
             n  = 0, 
             ARI = 0, 
             timeT = 0, # total run time
             timeS = -1) # sub-sample run time

idx = 1
for(n in n.seq){
  Ctrue = c(rep(1, n*prop[1]), rep(2, n*prop[2]))
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- n
    
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_unif)
    df[idx,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_deg, probs=deg)
    df[idx+1,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=edge_unif)
    df[idx+2,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=bfs_nodes)
    df[idx+3,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=dfs_nodes)
    df[idx+4,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=random_node_neigh)
    df[idx+5,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=rand_walk)
    df[idx+6,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    G <- as.undirected(graph_from_adjacency_matrix(A))
    df[idx+7,5] <- as.numeric(system.time(out <- cluster_fast_greedy(G)$membership)[3])
    df[idx+7,4] <- adj.rand.index(Ctrue, out)
    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_pace_vary_n.RData")
  }
  print(n)
}


load("sims_pace_vary_n.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))

# Just go up to n=8000
df_plot <- df %>% group_by(method, n) %>% 
  summarize(ARI=mean(ARI), timeT=mean(timeT), timeS=mean(timeS))

pM2 <- ggplot(df_plot, aes(x=n, y=ARI, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1000, y=0.99, label="(2)", size=6)
pM2

pT2 <- ggplot(df_plot, aes(x=n, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1000, y=40, label="(2)", size=6)+
  ylab("Total time (s)")
pT2

df_plot <- df_plot[df_plot$method!="FULL", ]
pS2 <- ggplot(df_plot, aes(x=n, y=timeS, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1200, y=0.010, label="(2)", size=6)+
  ylab("Sub-graph detection time (s)")
pS2



########## Fixed p11 - p12, increasing n (increasing B) 

n.seq = seq(1000, 8000, 1000)
p12 = 0.01
p11 = 0.04
prop = c(0.75, 0.25)
Ctrue = c(rep(1, n*prop[1]), rep(2, n*prop[2]))

df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(n.seq)),
             n  = 0, 
             ARI = 0, 
             timeT = 0, # total run time
             timeS = -1) # sub-sample run time

idx = 1
for(n in n.seq){
  B = n / 2
  Ctrue = c(rep(1, n*prop[1]), rep(2, n*prop[2]))
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- n
    
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_unif)
    df[idx,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_deg, probs=deg)
    df[idx+1,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=edge_unif)
    df[idx+2,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=bfs_nodes)
    df[idx+3,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=dfs_nodes)
    df[idx+4,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=random_node_neigh)
    df[idx+5,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=rand_walk)
    df[idx+6,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    # G <- as.undirected(graph_from_adjacency_matrix(A))
    # df[idx+7,5] <- as.numeric(system.time(out <- cluster_fast_greedy(G)$membership)[3])
    # df[idx+7,4] <- adj.rand.index(Ctrue, out)
    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_pace_vary_n_increaseB.RData")
  }
  print(n)
}


load("sims_pace_vary_n_increaseB.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))

df_plot <- df %>% group_by(method, n) %>% 
  summarize(ARI=mean(ARI), timeT=mean(timeT), timeS=mean(timeS))

# Load results from Setting 2 for FULL
load("sims_pace_vary_n.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))

df_plot2 <- df %>% group_by(method, n) %>% 
  summarize(ARI=mean(ARI), timeT=mean(timeT), timeS=mean(timeS))

df_plot[df_plot$method=="FULL", ] = df_plot2[df_plot2$method=="FULL", ]


pM3 <- ggplot(df_plot, aes(x=n, y=ARI, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1100, y=0.99, label="(3)", size=6)
pM3

pT3 <- ggplot(df_plot, aes(x=n, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1100, y=152, label="(3)", size=6)+
  ylab("Total time (s)")
pT3

df_plot <- df_plot[df_plot$method!="FULL", ]
pS3 <- ggplot(df_plot, aes(x=n, y=timeS, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=1100, y=0.0075, label="(3)", size=6)+
  ylab("Sub-graph detection time (s)")
pS3



########## Fixed p11 - p12, n, increasing proportion of nodes in community 1 
n = 5000
B = 1000
p12 = 0.01
p11 = 0.04
prop.seq = c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(prop.seq)),
             prop  = 0, 
             ARI = 0, 
             timeT = 0, # total run time
             timeS = -1) # sub-sample run time

idx = 1
for(prop in prop.seq){
  Ctrue = c(rep(1, n*prop), rep(2, n - n*prop))
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, c(prop, 1-prop))
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- prop
    
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_unif)
    df[idx,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_deg, probs=deg)
    df[idx+1,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=edge_unif)
    df[idx+2,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=bfs_nodes)
    df[idx+3,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=dfs_nodes)
    df[idx+4,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=random_node_neigh)
    df[idx+5,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=rand_walk)
    df[idx+6,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    G <- as.undirected(graph_from_adjacency_matrix(A))
    df[idx+7,5] <- as.numeric(system.time(out <- cluster_fast_greedy(G)$membership)[3])
    df[idx+7,4] <- adj.rand.index(Ctrue, out)
    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_pace_vary_prop.RData")
  }
  print(prop)
}


load("sims_pace_vary_prop.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))


df_plot <- df %>% group_by(method, prop) %>% summarize(ARI=mean(ARI), timeT=mean(timeT), timeS=mean(timeS))

pM4 <- ggplot(df_plot, aes(x=prop, y=ARI, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  xlab(expression(alpha))+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=0.52, y=0.82, label="(4)", size=6)
pM4

pT4 <- ggplot(df_plot, aes(x=prop, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  xlab(expression(alpha))+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=0.52, y=9, label="(4)", size=6)+
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
  annotate(geom="text", x=0.5, y=0.0025, label="(4)", size=6)+
  ylab("Sub-graph detection time (s)")
pS4


########## Fixed n and p, increasing B 

n = 5000
B.seq = 1000
p12 = 0.01
p11 = 0.04
B.seq = c(100, 250, 500, 1000, 2500, 5000, 10000)
prop = c(0.75, 0.25)
Ctrue = c(rep(1, n*prop[1]), rep(2, n*prop[2]))


df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(B.seq)),
             B  = 0, 
             ARI = 0, 
             timeT = 0, # total run time
             timeS = -1) # sub-sample run time

idx = 1
for(B in B.seq){
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- B
    
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_unif)
    df[idx,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_deg, probs=deg)
    df[idx+1,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=edge_unif)
    df[idx+2,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=bfs_nodes)
    df[idx+3,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=dfs_nodes)
    df[idx+4,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=random_node_neigh)
    df[idx+5,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time
    
    out <- pace(A, q = 250/n, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=rand_walk)
    df[idx+6,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    # G <- as.undirected(graph_from_adjacency_matrix(A))
    # df[idx+7,5] <- as.numeric(system.time(out <- cluster_fast_greedy(G)$membership)[3])
    # df[idx+7,4] <- adj.rand.index(Ctrue, out)
    
    idx  =  idx + 8
    print(sim)
    save(df, file="sims_pace_n5000_vary_B.RData")
  }
  print(B)
}



load("sims_pace_n5000_vary_B.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))
df_plot <- df %>% group_by(method, B) %>% 
  summarize(ARI=mean(ARI), timeT=mean(timeT), timeS=mean(timeS))


# Load results from Setting 1 for FULL
load("sims_pace_n5000_vary_p.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))

df_plot2 <- df %>% group_by(method, p11) %>% 
  summarize(ARI=mean(ARI), timeT=mean(timeT), timeS=mean(timeS))
df_plot2 <- df_plot2[as.logical((df_plot2$method=="FULL")*(df_plot2$p11==0.04)), ]

df_plot$ARI[df_plot$method=="FULL"] <- df_plot2$ARI
df_plot$timeT[df_plot$method=="FULL"] <- df_plot2$timeT

pM5 <- ggplot(df_plot, aes(x=B, y=ARI, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=100, y=0.85, label="(5)", size=6)
pM5

pT5 <- ggplot(df_plot, aes(x=B, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=500, y=80, label="(5)", size=6)+
  ylab("Total time (s)")
pT5

df_plot <- df_plot[df_plot$method!="FULL", ]
pS5 <- ggplot(df_plot, aes(x=B, y=timeS, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=500, y=0.0025, label="(5)", size=6)+
  ylab("Sub-graph detection time (s)")
pS5



########## Fixed n and p, increasing q 

n=5000
p12 = 0.01
p11 = 0.05
q.seq = seq(100, 500, 100) / n
prop = c(0.75, 0.25)
Ctrue = c(rep(1, n*prop[1]), rep(2, n*prop[2]))

B = 1000

df =  tibble(iter = 0, 
             method = rep(c("RN", "DN", "RE", "BFS", "DFS", "RNN", "RW", "FULL"), n.iter*length(q.seq)),
             q  = 0, 
             ARI = 0, 
             timeT = 0, # total run time
             timeS = -1) # sub-sample run time

idx = 1
for(q in q.seq){
  for(sim in 1:n.iter){
    
    A <- generateA(n, p11, p12, prop)
    df[idx:(idx+7), 1] <- sim
    df[idx:(idx+7), 3] <- q
    
    out <- pace(A, q = q, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_unif)
    df[idx,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx,5] <- out$total_time
    df[idx,6] <- out$sub_time
    
    deg = colSums(A)
    out <- pace(A, q = q, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=node_deg, probs=deg)
    df[idx+1,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+1,5] <- out$total_time
    df[idx+1,6] <- out$sub_time
    
    out <- pace(A, q = q, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=edge_unif)
    df[idx+2,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+2,5] <- out$total_time
    df[idx+2,6] <- out$sub_time
    
    out <- pace(A, q = q, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=bfs_nodes)
    df[idx+3,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+3,5] <- out$total_time
    df[idx+3,6] <- out$sub_time
    
    out <- pace(A, q = q, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=dfs_nodes)
    df[idx+4,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+4,5] <- out$total_time
    df[idx+4,6] <- out$sub_time
    
    out <- pace(A, q = q, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=random_node_neigh)
    df[idx+5,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+5,5] <- out$total_time
    df[idx+5,6] <- out$sub_time
    
    out <- pace(A, q = q, B = B, no_comms = 2, clustering = cluster_fast_greedy, sampler=rand_walk)
    df[idx+6,4] <- adj.rand.index(Ctrue, out$labels)
    df[idx+6,5] <- out$total_time
    df[idx+6,6] <- out$sub_time
    
    # G <- as.undirected(graph_from_adjacency_matrix(A))
    # df[idx+7,5] <- as.numeric(system.time(out <- cluster_fast_greedy(G)$membership)[3])
    # df[idx+7,4] <- adj.rand.index(Ctrue, out)   

    idx  =  idx + 8
    print(sim)
    save(df, file="sims_pace_n5000_vary_q.RData")
  }
  print(q)
}


load("sims_pace_n5000_vary_q.RData")

df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))
df_plot <- df %>% group_by(method, q) %>% 
  summarize(ARI=mean(ARI), timeT=mean(timeT), timeS=mean(timeS))


# Load results from Setting 1 for FULL
load("sims_pace_n5000_vary_p.RData")
df$method <- factor(df$method, levels=c("BFS", "DFS", "DN", "RE", "RN", "RNN", "RW", "FULL"))

df_plot2 <- df %>% group_by(method, p11) %>% 
  summarize(ARI=mean(ARI), timeT=mean(timeT), timeS=mean(timeS))
df_plot2 <- df_plot2[as.logical((df_plot2$method=="FULL")*(df_plot2$p11==0.05)), ]

df_plot$ARI[df_plot$method=="FULL"] <- df_plot2$ARI
df_plot$timeT[df_plot$method=="FULL"] <- df_plot2$timeT


pM6 <- ggplot(df_plot, aes(x=q*n, y=ARI, color=method, linetype=method))+
  xlab("qn")+
  geom_line(linewidth=1.5)+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=110, y=0.90, label="(6)", size=6)
pM6

pT6 <- ggplot(df_plot, aes(x=q*n, y=timeT, color=method, linetype=method))+
  geom_line(linewidth=1.5)+
  xlab("qn")+
  scale_colour_manual(values = myColors)+
  guides(color=guide_legend(title="Sub-sampler"), linetype=guide_legend(title="Sub-sampler"))+
  theme_bw()+
  theme(text = element_text(size = 16))+
  annotate(geom="text", x=110, y=32, label="(6)", size=6)+
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
  annotate(geom="text", x=110, y=0.02, label="(6)", size=6)+
  ylab("Sub-graph detection time (s)")
pS6

ggarrange(pM1, pM2, pM3, pM4, pM5, pM6, ncol=2, nrow=3, common.legend = T, legend="bottom")

ggsave(file="pace_sims.pdf",
       height = 8,
       width = 8,
       units="in")


ggarrange(pT1, pT2, pT3, pT4, pT5, pT6, ncol=2, nrow=3, common.legend = T, legend="bottom")

ggsave(file="pace_sims_timeT.pdf",
       height = 8,
       width = 8,
       units="in")

ggarrange(pS1, pS2, pS3, pS4, pS5, pS6, ncol=2, nrow=3, common.legend = T, legend="bottom")

ggsave(file="pace_sims_timeS.pdf",
       height = 8,
       width = 8,
       units="in")









