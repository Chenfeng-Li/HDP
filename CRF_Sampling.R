###* Load the target data $X_{ji}$ and encode the conditional functions *###
#### Loading the given data
raw.data<-read.csv("final_data.csv")
attach(raw.data)

## pre-processing data
X1 = logArea[group == 'immuno']
X2 = logArea[group == 'cryo']
X = logArea
all_i.1 = length(X1)
all_i.2 = length(X2)
all_i = length(X)

## set up this HDP
tau.phi = 1
tau.x = 1
alpha0 = 1
gam = 1

## compute the conditional density x_ji|x_-ji
f_k = function(i, x_list, k, k_list){ 
  x = x_list[i]
  k_list_cut = k_list[-i]
  x_list_cut = x_list[-i]
  n.k = sum(k_list_cut == k)
  x_ji_list = x_list_cut[k_list_cut == k]
  factor = sqrt((tau.x^2*n.k+tau.phi^2)/(tau.x^2*(n.k+1)+tau.phi^2))
  power = tau.x^2*x^2+(tau.x^2*sum(x_ji_list))^2/(tau.x^2*n.k+tau.phi^2)-
    (tau.x^2*(x+sum(x_ji_list)))^2/(tau.x^2*(n.k+1)+tau.phi^2)
  result = (tau.x/sqrt(2*pi))*factor*exp(-0.5*power)
  return(result)
}

f_k_new = function(i, x_list){
  x = x_list[i]
  factor = sqrt(1/(tau.x^2+tau.phi^2))
  power = (tau.x^2*tau.phi^2*x^2)/(tau.x^2+tau.phi^2)
  result = (tau.x*tau.phi/sqrt(2*pi))*factor*exp(-0.5*power) 
  return(result)
}


###* Likelihood function conditional on new t *###
P_t_new = function(i, x_list, t_list, k_list){ 
  K = max(k_list)
  k_list1 = k_list[1:all_i.1]
  k_list2 = k_list[(all_i.1+1):all_i]
  t_list1 = t_list[1:all_i.1]
  t_list2 = t_list[(all_i.1+1):all_i]
  m = length(unique(t_list1)) + length(unique(t_list2)) 
  prod = 0
  for (k in 1:K){
    rec1 = t_list1[k_list1 == k]
    val1 = length(unique(rec1))
    rec2 = t_list2[k_list2 == k]
    val2 = length(unique(rec2))
    prod = prod + (val1+val2)*f_k(i, x_list, k, k_list)
  }
  result = (prod+gam*f_k_new(i, x_list))/(m+gam)
  return(result)
}


###* Sampling t and k from the full conditional functions (one step) *###
update_t = function(i, x_list, t_list, k_list){ 
  if (i <= all_i.1){
    x = x_list[i]
    k_list1 = k_list[1:all_i.1] 
    k_list_cut1 = k_list1[-i] 
    t_list1 = t_list[1:all_i.1] 
    t_list_cut1 = t_list1[-i] 
    t_unique1 = unique(t_list_cut1) 
    n1 = length(t_unique1)
    T1 = max(abs(t_unique1))
    prob = rep(0, n1+1)
    for (j in 1:n1){
      n = sum(t_list_cut1 == t_unique1[j])
      k = k_list_cut1[t_list_cut1 == t_unique1[j]][1]
      f = f_k(i, x_list, k, k_list)
      prob[j] = n*f
    }
    prob[n1+1] = alpha0 * P_t_new(i, x_list, t_list, k_list)
    prob = prob/sum(prob)
    new_t = sample(c(t_unique1,T1+1), 1, replace = TRUE, prob = prob) 
    if (new_t == T1+1){
      ## sample new_k
      K = max(k_list)
      k_list1 = k_list[1:all_i.1]
      k_list2 = k_list[(all_i.1+1):all_i]
      t_list1 = t_list[1:all_i.1]
      t_list2 = t_list[(all_i.1+1):all_i]
      k_unique1 = unique(k_list)
      n.k1 = length(k_unique1)
      K1 = max(k_unique1)
      prob_k = rep(0, n.k1+1) 
      for (j in 1:n.k1){
        rec1 = t_list1[k_list1 == k_unique1[j]]
        val1 = length(unique(rec1))
        rec2 = t_list2[k_list2 == k_unique1[j]]
        val2 = length(unique(rec2))
        prob_k[j] = (val1+val2)*f_k(i, x_list, k_unique1[j], k_list)
      }
      prob_k[n.k1+1] = gam * f_k_new(i, x_list)
      prob_k = prob_k/sum(prob_k)
      new_k = sample(c(k_unique1,K1+1), 1, replace = TRUE, prob = prob_k)
    }else{
      new_k = k_list_cut1[t_list_cut1 == new_t][1]
    }
  }else{
    #prob = rep(0, T+1)
    x = x_list[i]
    k_list2 = k_list[(all_i.1+1):all_i] 
    k_list_cut2 = k_list2[-(i-all_i.1)] 
    t_list2 = t_list[(all_i.1+1):all_i] 
    t_list_cut2 = t_list2[-(i-all_i.1)] 
    t_unique2 = unique(t_list_cut2)
    n2 = length(t_unique2)
    T2 = max(abs(t_unique2))
    prob = rep(0, n2+1)
    for (j in 1:n2){
      n = sum(t_list_cut2 == t_unique2[j])
      k = k_list_cut2[t_list_cut2 == t_unique2[j]][1]
      f = f_k(i, x_list, k, k_list)
      prob[j] = n*f
    }
    prob[n2+1] = alpha0 * P_t_new(i, x_list, t_list, k_list)
    prob = prob/sum(prob)
    new_t = sample(c(t_unique2, -T2-1), 1, replace = TRUE, prob = prob) 
    if (new_t == -T2-1){
      ## sample new_k
      K = max(k_list)
      k_list1 = k_list[1:all_i.1]
      k_list2 = k_list[(all_i.1+1):all_i] 
      t_list1 = t_list[1:all_i.1]
      t_list2 = t_list[(all_i.1+1):all_i] 
      k_unique2 = unique(k_list)
      n.k2 = length(k_unique2)
      K2 = max(k_unique2)
      prob_k = rep(0, n.k2+1)
      for (j in 1:n.k2){
        rec1 = t_list1[k_list1 == k_unique2[j]]
        val1 = length(unique(rec1))
        rec2 = t_list2[k_list2 == k_unique2[j]]
        val2 = length(unique(rec2))
        prob_k[j] = (val1+val2)*f_k(i, x_list, k_unique2[j], k_list)
      }
      prob_k[n.k2+1] = gam * f_k_new(i, x_list)
      prob_k = prob_k/sum(prob_k)
      new_k = sample(c(k_unique2,K2+1), 1, replace = TRUE, prob = prob_k)
    }else{
      new_k = k_list_cut2[t_list_cut2 == new_t][1]
    }
  }
  t_list[i] = new_t
  k_list[i] = new_k
  return(c(t_list, k_list))
}


###* MCMC sampling for 1500 iterations *###
set.seed(123)
## Initialize
t_list = rep(1, all_i) 
t_list[(all_i.1+1):all_i] = rep(-1, all_i.2) 
k_list = rep(1, all_i)
x_list = X
## save the sampling result
t_trace = t_list
k_trace = k_list
N = 1500
for (n in 1:N){
  for (i in 1:all_i){
    tk_list = update_t(i, x_list, t_list, k_list) 
    t_list = tk_list[1:all_i]
    k_list = tk_list[(all_i+1):(2*all_i)]
  }
  t_trace = cbind(t_trace, t_list)
  k_trace = cbind(k_trace, k_list)
}

###* Sample result of k_jt_ji *###
par(mfrow=c(1,1))
plot(k_trace[1,], type='l', main = "trace plot of k_11", xlab='iteration', ylab='k_11')

par(mfrow=c(1,1))
plot(k_trace[(all_i.1+1),], type='l',main = "trace plot of k_21",xlab="iteration", ylab='k_21')


###* Boxplot of the observations X_j from the two groups j=1,2 *###
###* Check the clustering result *###
par(mfrow=c(1,1))
boxplot(X1,X2, main="original data",xlab="cluster",xaxt="n")
axis(1, at = c(1,2),label=c("immuno","cryo"), las=1)


###* Cluster distributions of X in sample iteration = 1500, 1400,1300,1200 *###
par(mfrow=c(2,2))
col_num = 1501
boxplot(X[k_trace[,col_num] == 31], X[k_trace[,col_num] == 32], X[k_trace[,col_num] == 33],
        main="iteration=1500",xlab="cluster",xaxt="n")
axis(1, at = c(1,2,3), las=1)
col_num = 1401
boxplot(X[k_trace[,col_num] == 34], X[k_trace[,col_num] == 31], X[k_trace[,col_num] == 36],
        main="iteration=1400",xlab="cluster",xaxt="n")
axis(1, at = c(1,2,3), las=1)
col_num = 1301
boxplot(X[k_trace[,col_num] == 34], X[k_trace[,col_num] == 31], X[k_trace[,col_num] == 35],
        main="iteration=1300",xlab="cluster",xaxt="n")
axis(1, at = c(1,2,3), las=1)
col_num = 1201
boxplot(X[k_trace[,col_num] == 29], X[k_trace[,col_num] == 31],main="iteration=1200",
        xlab="cluster",xaxt="n")
axis(1, at = c(1,2), las=1)


###* Number of cluster at each iteration *###
par(mfrow=c(1,1))
l = ncol(k_trace) 
num_list = rep(0, l) 
for (i in 1:l){
  num_list[i] = length(unique(k_trace[,i]))
}
plot(1:l, num_list, 'l', main="number of clusters at each iteration",
     xlab='Iteration', ylab='number of clusters')


###* Visualize the cluster result through heatmap.*###
###* Remark: points 1,...,71 ($X_{1i}, i=1,2,...,71$) come from group 1 *###
###* and points 72,...,119 ($X_{2i}, i=1,2,...,48$) come from group 2 *###
simu_mat = matrix(0, all_i, all_i) 
for (i in 1:all_i){
  for (j in 1:all_i){
    simu_mat[i,j] = sum(k_trace[i,] == k_trace[j,]) / l
  }
}
par(mfrow=c(1,2))
heatmap(simu_mat, Rowv = NA, Colv = NA, main="original order")
heatmap(simu_mat, main="altered order")
