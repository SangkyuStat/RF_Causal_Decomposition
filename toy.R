##############################################################################################################

# Title: Toy example version of the simulation (just one replication of the Monte Carlo simulation). 

# The original simulation contains the data that is not shareable (the data is also very huge), 
# and also the simulation data and setting were much more complex, so we had to use the HPC with 
# a lot much computation time more than a week on HPC. 

# This version is more simpler (less number of tuning parameter girds), and favorable to users who wants 
# to understand the structure and codes of our algorithms introduced in the paper. 

##############################################################################################################
library(grf)
library(dplyr)
library(caret)
library(survival)
library(rlist)
##############################################################################################################
set.seed(2026)

n_total <- 500
sample_ratio <- 0.8
scale <- sqrt(sample_ratio)
n_sample <- n_total*sample_ratio

X1 <- sample(20:100, n_total, replace = T) # A1
X2 <- sample(1:2, n_total, replace = T) %>% as.factor() # A1

X3 <- sample(1:2, n_total, replace = T) %>% as.factor() # A2
X4 <- sample(0:1, n_total, replace = T) %>% as.factor() # A2

X5 <- rbinom(n_total, 1, 0.6) %>% as.factor() 

M1 <- sample(1:4, n_total, replace = T) %>% as.factor() # insurance
M2 <- sample(1:3, n_total, replace = T) %>% as.factor() # stage

full_data <- data.frame(X1, X2, X3, X4, X5, M1, M2)

full_data_x <- model.matrix(~., full_data[,-5])[,-1] # remove one variable to generate Y that is not related with the specific variable

rep.max <- 101
res <- list()

set.seed(seed_num*100000)
beta <- c(-0.05, 0.5, -0.5, 0.5, 1, 1.2, 1.5, -1, 1) # true beta vector, we are assuming that M1 is a significant predictor while M2 is not
####### 
err <- runif(length(beta), -0.05, 0.05)
beta_censor <- beta + err

u_vec <- runif(n_total)
lambda <- 0.1
timecens <- - log(u_vec)/(lambda*exp(full_data_x%*%beta))
u_vec_censor <- runif(n_total)
lambda_censor <- 0.2
censor <- -log(u_vec_censor)/(lambda*exp(full_data_x%*%beta_censor))

DiedCancSpecCens <- as.numeric(I(timecens <= censor))
timecens <- apply(cbind(timecens, censor), 1, min)

# mean(DiedCancSpecCens)
####### 

for(rep in 1:rep.max){
  # cat(rep, "\n")
  set.seed(seed_num*10000 + rep)
  if(rep == 1){
    idx <- sample(1:n_total, size = n_total)
  } else {
    idx <- sample(1:n_total, size = n_sample)
  }
  
  no_missing_total_male_new <- cbind(full_data[idx, ], timecens[idx], DiedCancSpecCens[idx])
  
  Y_total <- no_missing_total_male_new$timecens
  X_total <- no_missing_total_male_new[, c("X1", "X2", "X3", "X4", "X5", "M1", "M2")]
  D_total <- no_missing_total_male_new$DiedCancSpecCens
  
  number_tree <- 400
  small_tree <- 400
  
  for(i in 1:2){
    if(i == 1){
      set.seed(2000 + i + rep)
      d_total = data.frame(Y_total, X_total, D_total)
      
      # f1
      data = d_total[, !(names(d_total) %in% c("X5", "Y_total", "D_total"))]
      
      dummy = dummyVars(~ X2, data = data)
      dummy = predict(dummy, data)
      data = as.matrix(cbind(d_total$X1, dummy))
      fit = probability_forest(X = data, Y = d_total$X5, seed = 2000+i+rep,
                               num.trees = small_tree, honesty = honest_input,
                               num.threads = 16)
      
      f_hat_r1_c = predict(fit, newdata = data)$predictions[,2]
      f_hat_r0_c = 1-f_hat_r1_c
      
      # f2
      data = d_total[, !(names(d_total) %in% c("X5", "Y_total", "D_total"))]
      
      dummy = dummyVars(~ X2 + X3 + X4, data = data)
      dummy = predict(dummy, data)
      data = as.matrix(cbind(d_total$X1, dummy))
      fit = probability_forest(X = data, Y = d_total$X5, seed = 2001+i+rep,
                               num.trees = small_tree, honesty = honest_input,
                               num.threads = 16)
      
      f_hat_r1_x0 =  predict(fit, newdata = data)$predictions[,2]
      f_hat_r0_x0 = 1-f_hat_r1_c
      
      # f3
      data = d_total[, !(names(d_total) %in% c("X5", "Y_total", "D_total"))]
      
      dummy = dummyVars(~ X2 + X3 + X4 + M1, data = data)
      dummy = predict(dummy, data)
      data = as.matrix(cbind(d_total$X1, dummy))
      fit = probability_forest(X = data, Y = d_total$X5, seed = 2002+i+rep,
                               num.trees = small_tree, honesty = honest_input,
                               num.threads = 16)
      
      f_hat_r1_x0_m1 =  predict(fit, newdata = data)$predictions[,2]
      f_hat_r2_x0_m1 = 1-f_hat_r1_x0_m1
      
      # f4
      data = d_total[, !(names(d_total) %in% c("X5", "Y_total", "D_total"))]
      
      dummy = dummyVars(~ X2 + X3 + X4 + M1 + M2, data = data)
      dummy = predict(dummy, data)
      data = as.matrix(cbind(d_total$age, dummy))
      fit = probability_forest(X = data, Y = d_total$X5, seed = 2003+i+rep,
                               num.trees = small_tree, honesty = honest_input,
                               num.threads = 16)
      
      f_hat_r1_x0_m1_m2 = predict(fit, newdata = data)$predictions[,2]
      f_hat_r2_x0_m1_m2 = 1-f_hat_r1_x0_m1_m2
      
      # omega(r0 = 0, r1 = 1, r2 = 0)
      w1 = (f_hat_r2_x0_m1_m2 / f_hat_r2_x0_m1) * (f_hat_r1_x0_m1/f_hat_r1_x0) * (f_hat_r0_x0/f_hat_r0_c)
      w1 = w1/sum(w1)
      
      # omega(r0 = 0, r1 = 1, r2 = 1)
      w2 = (f_hat_r1_x0_m1_m2 / f_hat_r1_x0_m1) * (f_hat_r1_x0_m1/f_hat_r1_x0) * (f_hat_r0_x0/f_hat_r0_c)
      w2 = w2/sum(w2)
      
      # omega(r0 = 0, r1 = 0, r2 = 1)
      w3 = (f_hat_r1_x0_m1_m2 / f_hat_r1_x0_m1) * (f_hat_r2_x0_m1/f_hat_r0_x0) * (f_hat_r0_x0/f_hat_r0_c)
      w3 = w3/sum(w3)
      
      # omega(r0 = 0, r1 = 0, r2 = 0)
      w4 = (f_hat_r2_x0_m1_m2 / f_hat_r2_x0_m1) * (f_hat_r2_x0_m1/f_hat_r0_x0) * (f_hat_r0_x0/f_hat_r0_c)
      w4 = w4/sum(w4)
    } else {
      set.seed(2100 + i +rep)
      
      d_total = data.frame(Y_total, X_total, D_total)
      
      data0 = d_total[d_total$X5 == 0, !(names(d_total) %in% c("X5", "Y_total", "D_total"))]
      
      dummy = dummyVars(~ X2 + X3 + X4 + M1 + M2, data = data0)
      dummy = predict(dummy, data0)
      data0 = as.matrix(cbind(d_total$X1[d_total$X5 == 0], dummy))
      
      data = d_total[, !(names(d_total) %in% c("X5", "Y_total", "D_total"))]
      
      dummy = dummyVars(~ X2 + X3 + X4 + M1 + M2, data = data)
      dummy = predict(dummy, data)
      data = as.matrix(cbind(d_total$X1, dummy))
      
      ##### tuning the parameters with 5-fold CV
      sample.fraction.seq = 0.5 
      honesty.fraction.seq = 0.5 
      # min.node.size.seq = c(15, 25, 40)
      min.node.size.seq = c(15, 25, 40)
      mtry.seq = c(3, 5, 7)
      train.result = array(0, dim = c(1,1,3,3))
      
      for(h in 1:length(sample.fraction.seq)){
        for(j in 1:length(honesty.fraction.seq)){
          for(q in 1:length(min.node.size.seq)){
            for(r in 1:length(mtry.seq)){
              model.train = survival_forest(X = data0,
                                            Y = Y_total[d_total$X5 == 0],
                                            D = D_total[d_total$X5 == 0],
                                            sample.fraction = sample.fraction.seq[h],
                                            honesty.fraction = honesty.fraction.seq[j],
                                            min.node.size = min.node.size.seq[q],
                                            mtry = mtry.seq[r],
                                            num.trees = number_tree,
                                            honesty = honest_input,
                                            seed = 2000+h+j+q+r+rep,
                                            num.threads = 16)
              s.pred.nelson.aalen <- predict(model.train, prediction.type = "Nelson-Aalen")
              chf.score <- rowSums(-log(s.pred.nelson.aalen$predictions))
              
              train.result[h,j,q,r] = concordance(Surv(Y_total[d_total$X5 == 0], D_total[d_total$X5 == 0]) ~ chf.score, reverse = TRUE)$concordance
            }
          }
        }
      }
      tune.par.idx = which(train.result == max(train.result), arr.ind = TRUE)
      tune.black.male = tune.par.idx
      
      model0 = survival_forest(X = data0,
                               Y = Y_total[d_total$X5 == 0],
                               D = D_total[d_total$X5 == 0],
                               sample.fraction = sample.fraction.seq[tune.par.idx[1]],
                               honesty.fraction = honesty.fraction.seq[tune.par.idx[2]],
                               min.node.size = min.node.size.seq[tune.par.idx[3]],
                               mtry = mtry.seq[tune.par.idx[4]],
                               num.trees = number_tree,
                               honesty = honest_input,
                               seed = 2200+i+rep,
                               num.threads = 16)
      
      expectation0 = grf:::predict.survival_forest(model0, newdata = data,
                                                   failure.times = c(100, 200, 300))
      
      sub_res1 = w1 * expectation0$predictions
      omega_b_1 = colSums(sub_res1)
      
      sub_res2 = w2 * expectation0$predictions
      omega_b_2 = colSums(sub_res2)
      
      sub_res3 = w3 * expectation0$predictions
      omega_b_3 = colSums(sub_res3)
      
      sub_res4 = w4 * expectation0$predictions
      omega_b_4 = colSums(sub_res4)
      
      res1 = omega_b_1
      res2 = omega_b_2
      res3 = omega_b_3
      res4 = omega_b_4
    }
  }
  res[[rep]] = cbind(res1, res2, res3, res4)
}

res_final1 = list.rbind(lapply(res, function(x) x[,1] - x[,4]))
res_final2 = list.rbind(lapply(res, function(x) x[,2] - x[,4]))
res_final3 = list.rbind(lapply(res, function(x) x[,3] - x[,4]))

res_list = list()

res_list$hyp_mat = matrix(0, ncol = 3, nrow = 3)

res_list$hyp_mat[1,1] = (res_final1[1,1] - sd(res_final1[-1,1]) * qnorm(0.975) * scale) <= 0 & (res_final1[1,1] + sd(res_final1[-1,1]) * qnorm(0.975) * scale) >= 0  
res_list$hyp_mat[2,1] = (res_final1[1,2] - sd(res_final1[-1,2]) * qnorm(0.975) * scale) <= 0 & (res_final1[1,2] + sd(res_final1[-1,2]) * qnorm(0.975) * scale) >= 0  
res_list$hyp_mat[3,1] = (res_final1[1,3] - sd(res_final1[-1,3]) * qnorm(0.975) * scale) <= 0 & (res_final1[1,3] + sd(res_final1[-1,3]) * qnorm(0.975) * scale) >= 0  

res_list$hyp_mat[1,3] = (res_final2[1,1] - sd(res_final2[-1,1]) * qnorm(0.975) * scale) <= 0 & (res_final2[1,1] + sd(res_final2[-1,1]) * qnorm(0.975) * scale) >= 0  
res_list$hyp_mat[2,3] = (res_final2[1,2] - sd(res_final2[-1,2]) * qnorm(0.975) * scale) <= 0 & (res_final2[1,2] + sd(res_final2[-1,2]) * qnorm(0.975) * scale) >= 0  
res_list$hyp_mat[3,3] = (res_final2[1,3] - sd(res_final2[-1,3]) * qnorm(0.975) * scale) <= 0 & (res_final2[1,3] + sd(res_final2[-1,3]) * qnorm(0.975) * scale) >= 0 

res_list$hyp_mat[1,2] = (res_final3[1,1] - sd(res_final3[-1,1]) * qnorm(0.975) * scale) <= 0 & (res_final3[1,1] + sd(res_final3[-1,1]) * qnorm(0.975) * scale) >= 0  
res_list$hyp_mat[2,2] = (res_final3[1,2] - sd(res_final3[-1,2]) * qnorm(0.975) * scale) <= 0 & (res_final3[1,2] + sd(res_final3[-1,2]) * qnorm(0.975) * scale) >= 0  
res_list$hyp_mat[3,2] = (res_final3[1,3] - sd(res_final3[-1,3]) * qnorm(0.975) * scale) <= 0 & (res_final3[1,3] + sd(res_final3[-1,3]) * qnorm(0.975) * scale) >= 0 

colnames(res_list$hyp_mat) = c("m1 - no", "m2 - no", "m1m2 - no")
rownames(res_list$hyp_mat) = c("100 days", "200 days", "300 days")

print(res_list)

if(honest_input){
  saveRDS(res_list, paste0("./res_files_honest/res_", seed_num, ".rds"))
} else {
  saveRDS(res_list, paste0("./res_files_nohonest/res_", seed_num, ".rds"))
}

