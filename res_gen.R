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

for(i in file_list){
  name = i
  temp = readRDS(name)
  if(i == file_list[1]){
    res = temp$hyp_mat
  } else {
    res = res + temp$hyp_mat
  }
}
res/length(file_list)


seed_num_seq <- 1:100
for(seed_num in seed_num_seq){
  honest_input <- TRUE
  cat("########## Monte Carlo Iteration Number:", seed_num, "#################### \n")
  source("toy.R")
}

folder_path <- paste0("./res_files_honest")
file_list <- list.files(folder_path, full.names = TRUE)

for(i in file_list){
  name = i
  temp = readRDS(name)
  if(i == file_list[1]){
    res = temp$hyp_mat
  } else {
    res = res + temp$hyp_mat
  }
}
res/length(file_list)
