##############################################################################################################

# Running the toy version of simulation

##############################################################################################################
rm(list=ls())
##############################################################################################################

seed_num_seq <- 1:100
for(seed_num in seed_num_seq){
  honest_input <- FALSE
  cat("\n########## Monte Carlo Iteration Number:", seed_num, "#################### \n")
  source("toy.R")
}

folder_path <- paste0("./res_files_nohonest")
file_list <- list.files(folder_path, full.names = TRUE)

na_nh <- 0 # since the number of samples is too small, there might be some cases with missing results, 
            # there was no missing in the original simulation
for(i in file_list){
  name = i
  temp = readRDS(name)
  if(i == file_list[1]){
    res = temp$hyp_mat
  } else {
    if(is.na(temp$hyp_mat[1,1])){
      na_nh = na_nh + 1
      next
    } else {
      res = res + temp$hyp_mat
    }
  }
}
res/(length(file_list) - na_nh)


seed_num_seq <- 1:100
for(seed_num in seed_num_seq){
  honest_input <- TRUE
  cat("########## Monte Carlo Iteration Number:", seed_num, "#################### \n")
  source("toy.R")
}

folder_path <- paste0("./res_files_honest")
file_list <- list.files(folder_path, full.names = TRUE)

na_h <- 0 # since the number of samples is too small, there might be some cases with missing results, 
          # there was no missing in the original simulation
for(i in file_list){
  name = i
  temp = readRDS(name)
  if(i == file_list[1]){
    res = temp$hyp_mat
  } else {
    if(is.na(temp$hyp_mat[1,1])){
      na_h = na_h + 1
      next
    } else {
      res = res + temp$hyp_mat
    }
  }
}
res/(length(file_list) - na_h)
