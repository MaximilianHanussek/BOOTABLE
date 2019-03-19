#Parse ThreadedBioBenchSuite output



chooseCRANmirror(ind = 33)

old_warnvalue <- getOption("warn")
options(warn = -1)

load.fun <- function(x) {
  x <- as.character(substitute(x))
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text=paste("require(", x, ")", sep="")))
  } else {
    update.packages() # recommended before installing so that dependencies are the latest version
    eval(parse(text=paste("install.packages('", x, "')", sep="")))
    eval(parse(text=paste("require(", x, ")", sep="")))
  }
} 

load.fun(grid)
load.fun(gridExtra)
load.fun(RColorBrewer)
load.fun(stringr)

workingdir <- "./"


args <- commandArgs(trailingOnly = TRUE)
argslen <- length(args)
if (argslen == 1){ 
  scaling_flag <- args[1]
} else {
  scaling_flag <- FALSE
}

# Get host information

cpu_info <- system("lscpu", intern = TRUE)
ram_info <- system("lsmem", intern = TRUE)
compiler_info <- system2("gcc", args = c("-v"), stderr = TRUE)
hostname <- system("hostname", intern = TRUE)

cpu_architecture <- str_squish(unlist(strsplit(cpu_info[grep("Architecture:", cpu_info)], split = ":")))[2]
cpu_vendor <- str_squish(unlist(strsplit(cpu_info[grep("Vendor ID:", cpu_info)], split = ":")))[2]
cpu_family <- str_squish(unlist(strsplit(cpu_info[grep("CPU family:", cpu_info)], split = ":")))[2]
cpu_model <- str_squish(unlist(strsplit(cpu_info[grep("Model:", cpu_info)], split = ":")))[2]
cpu_model_name <- unlist(strsplit(str_squish(unlist(strsplit(cpu_info[grep("Model name:", cpu_info)], split = ":")))[2], split = "@"))[1]
cpu_cache3_size <- paste0(str_squish(unlist(strsplit(cpu_info[grep("L3 cache:", cpu_info)], split = ":")))[2], "B")
cpu_cache2_size <- paste0(str_squish(unlist(strsplit(cpu_info[grep("L2 cache:", cpu_info)], split = ":")))[2], "B")
cpu_proc_speed <- toString(unlist(strsplit(str_squish(unlist(strsplit(cpu_info[grep("Model name:", cpu_info)], split = ":")))[2], split = "@ "))[2])
cpu_num_cores <- str_squish(unlist(strsplit(cpu_info[grep("CPU\\(s\\):", cpu_info)], split = ":")))[2]
ram_max_mem <- paste0(str_squish(unlist(strsplit(ram_info[grep("Total online memory:", ram_info)], split = ":")))[2], "B")
compiler_version <- paste(unlist(strsplit(compiler_info[grep("gcc-Version", compiler_info)], split = " "))[1], unlist(strsplit(compiler_info[7], split = " "))[2], collapse = " ")

if (cpu_cache3_size == "NAB") {
  system_information_row_names <- c("Hostname",
                                    "Architecture", 
                                    "Vendor ID", 
                                    "CPU Family", 
                                    "CPU Model", 
                                    "CPU Model Name", 
                                    "L2 Cache", 
                                    "Processor Speed",
                                    "CPU cores",
                                    "RAM",
                                    "Compiler")
  cpu_cache_size <- cpu_cache2_size
} else{
  system_information_row_names <- c("Hostname",
                                    "Architecture", 
                                    "Vendor ID", 
                                    "CPU Family", 
                                    "CPU Model", 
                                    "CPU Model Name", 
                                    "L3 Cache", 
                                    "Processor Speed",
                                    "CPU cores",
                                    "RAM",
                                    "Compiler")
  cpu_cache_size <- cpu_cache3_size
}
  
df_system_information <- data.frame(matrix(ncol = 1, nrow = 11))
colnames(df_system_information) <- c("System Information")
rownames(df_system_information) <- system_information_row_names
df_system_information[1,1] <- hostname
df_system_information[2,1] <- cpu_architecture
df_system_information[3,1] <- cpu_vendor
df_system_information[4,1] <- cpu_family
df_system_information[5,1] <- cpu_model
df_system_information[6,1] <- cpu_model_name
df_system_information[7,1] <- cpu_cache_size
df_system_information[8,1] <- cpu_proc_speed
df_system_information[9,1] <- cpu_num_cores
df_system_information[10,1] <- ram_max_mem
df_system_information[11,1] <- compiler_version


# Get all summary files, containing the time output
summary_file_paths <- list.files(path = workingdir, pattern = "^benchmark_summary_.*\\.txt", full.names = TRUE)
summary_file_names <- list.files(path = workingdir, pattern = "^benchmark_summary_.*\\.txt", full.names = FALSE)
number_of_summary_files <- length(summary_file_paths)
names_time_vector <- c("real", "user", "sys")

# Initialize scaling capabilities vector
scaling_cores_vector <- c()
scaling_mean_real_times_vector <- c()
scaling_number_of_used_tools <- c()

for (summary_file in summary_file_paths){
  used_tools <- c()
  used_replica <- c()
  real_values_all_vector <- c()
  user_values_all_vector <- c()
  sys_values_all_vector <- c()
  summary_name <- summary_file_names[which(summary_file == summary_file_paths)]
  summary_name <- strsplit(summary_name, "\\.")[[1]][1]
  summary_file_name <- paste(summary_name, ".pdf", sep = "")
  
  summary_file_lines <- readLines(con = summary_file)
  replica_entries_line_numbers <- grep("Replica", summary_file_lines)
  real_entries_line_numbers <- grep ("real ", summary_file_lines)
  user_entries_line_numbers <- grep ("user ", summary_file_lines)
  sys_entries_line_numbers <- grep ("sys ", summary_file_lines)
  
  number_of_cores <- as.integer(strsplit(summary_file_lines[1], "\\s+")[[1]][4])
  scaling_cores_vector <- c(scaling_cores_vector, number_of_cores)
  
  for (replica_entry in replica_entries_line_numbers){
    used_replica <- c(used_replica, strsplit(summary_file_lines[replica_entry], "\\s+")[[1]][1])
    used_tools <- c(used_tools, strsplit(summary_file_lines[replica_entry], "\\s+")[[1]][2])
  }
  
  used_tools_unique <- unique(used_tools)
  number_of_used_tools <- length(used_tools_unique)
  scaling_number_of_used_tools <- number_of_used_tools
  used_replica_unique <- unique(used_replica)
  number_of_used_replica <- length(used_replica_unique)
  
  df_all <- data.frame(matrix(data = list(), ncol = number_of_used_tools, nrow = number_of_used_replica))
  df_all_printable <- data.frame(matrix(ncol = number_of_used_tools, nrow = number_of_used_replica))
  colnames(df_all) <- used_tools_unique
  rownames(df_all) <- used_replica_unique
  colnames(df_all_printable) <- used_tools_unique
  rownames(df_all_printable) <- used_replica_unique
  
  i <- 1
  j <- 1
  
  for (entry in 1:length(real_entries_line_numbers)){
    real_entry <- as.numeric(strsplit(summary_file_lines[real_entries_line_numbers[entry]], "\\s+")[[1]][2])
    user_entry <- as.numeric(strsplit(summary_file_lines[user_entries_line_numbers[entry]], "\\s+")[[1]][2])
    sys_entry  <- as.numeric(strsplit(summary_file_lines[sys_entries_line_numbers[entry]], "\\s+")[[1]][2])
    time_vector <- c(real_entry, user_entry, sys_entry)
    df_all[[i, j]] <- time_vector
    df_all_printable[[i, j]] <- paste(time_vector, collapse = " / ")
    
    if (j >= number_of_used_tools){
    i <- i + 1
    j <- 1
    }
    else {
      j <- j + 1
    }
  }
  
  used_time_values <- c("average real", "average user", "average sys")
  df_toolwise <- data.frame(matrix(ncol = number_of_used_tools, nrow = 3))
  colnames(df_toolwise) <- used_tools_unique
  rownames(df_toolwise) <- used_time_values
  real_mean_values_tools_vector <- c() 
  user_mean_values_tools_vector <- c()
  sys_mean_values_tools_vector <- c()
  
  for (tool in 1:number_of_used_tools) {
    real_values_vector <- c()
    user_values_vector <- c()
    sys_values_vector  <- c()
    
    for (value in 1:number_of_used_replica){
      real_values_vector <- c(real_values_vector, df_all[[tool]][[value]][1])
      
      user_values_vector <- c(user_values_vector, df_all[[tool]][[value]][2])
      
      sys_values_vector <- c(sys_values_vector, df_all[[tool]][[value]][3])
    }
    df_toolwise[[1, tool]] <- format(round(mean(real_values_vector), digits = 2), nsmall = 2)
    df_toolwise[[2, tool]] <- format(round(mean(user_values_vector), digits = 2), nsmall = 2)
    df_toolwise[[3, tool]] <- format(round(mean(sys_values_vector), digits = 2), nsmall = 2)
    
    real_mean_values_tools_vector <- c(real_mean_values_tools_vector, df_toolwise[[1, tool]])
    user_mean_values_tools_vector <- c(user_mean_values_tools_vector, df_toolwise[[2, tool]])
    sys_mean_values_tools_vector <- c(sys_mean_values_tools_vector, df_toolwise[[3, tool]])
    
    
    
    real_values_all_vector <- c(real_values_all_vector, real_values_vector)
    user_values_all_vector <- c(user_values_all_vector, user_values_vector)
    sys_values_all_vector  <- c(sys_values_all_vector, sys_values_vector) 
  }
  
  scaling_mean_real_times_vector <- c(scaling_mean_real_times_vector, real_mean_values_tools_vector)
  
  df_replicawise <- data.frame(matrix(ncol = number_of_used_replica, nrow = 3))
  colnames(df_replicawise) <- used_replica_unique
  rownames(df_replicawise) <- used_time_values
  
  real_mean_values_replica_vector <- c() 
  user_mean_values_replica_vector <- c()
  sys_mean_values_replica_vector <- c()
  
  for (replica in 1:number_of_used_replica) {
    real_values_vector <- c()
    user_values_vector <- c()
    sys_values_vector  <- c()
    
    for (value in 1:number_of_used_tools){
      real_values_vector <- c(real_values_vector, df_all[[value]][[replica]][1])
      
      user_values_vector <- c(user_values_vector, df_all[[value]][[replica]][2])
      
      sys_values_vector <- c(sys_values_vector, df_all[[value]][[replica]][3])
    }
    df_replicawise[[1, replica]] <- format(round(mean(real_values_vector), digits = 2), nsmall = 2)
    df_replicawise[[2, replica]] <- format(round(mean(user_values_vector), digits = 2), nsmall = 2)
    df_replicawise[[3, replica]] <- format(round(mean(sys_values_vector), digits = 2), nsmall = 2)
    
    real_mean_values_replica_vector <- c(real_mean_values_replica_vector, df_replicawise[[1, replica]])
    user_mean_values_replica_vector <- c(user_mean_values_replica_vector, df_replicawise[[2, replica]])
    sys_mean_values_replica_vector <- c(sys_mean_values_replica_vector, df_replicawise[[3, replica]])
  }
  total_time_real <- format(sum(real_values_all_vector), nsmall = 2)
  total_time_user <- format(sum(user_values_all_vector), nsmall = 2)
  total_time_sys  <- format(sum(sys_values_all_vector), nsmall = 2)

  t0 <- tableGrob(df_system_information, rows = rownames(df_system_information), cols = colnames(df_system_information))
  t1 <- tableGrob(t(df_all_printable), rows = colnames(df_all_printable), cols = rownames(df_all_printable))
  t2 <- tableGrob(t(df_toolwise), rows = colnames(df_toolwise), cols = rownames(df_toolwise))
  t3 <- tableGrob(t(df_replicawise), rows = colnames(df_replicawise), cols = rownames(df_replicawise))
  
  
  title <- textGrob("Statistical Report of BOOTABLE", gp = gpar(fontsize=20,font=2), just = "center")
  subtitle <- textGrob(date(), gp = gpar(fontsize = 14, font = 1), just = "center")
  margin <- unit(0.5, "line")
  
  real_tools_percentage_vector <- format(round(prop.table(as.numeric(real_mean_values_tools_vector))*100, digits = 2))
  legend_vector_tools_real <- unlist(strsplit(paste(used_tools_unique, real_tools_percentage_vector, sep = " ", collapse = "," ), split = ","))
  legend_vector_tools_real <- unlist(strsplit(paste(legend_vector_tools_real, "%", sep = "", collapse = ","), split = ","))
  
  user_tools_percentage_vector <- format(round(prop.table(as.numeric(user_mean_values_tools_vector))*100, digits = 2))
  legend_vector_tools_user <- unlist(strsplit(paste(used_tools_unique, user_tools_percentage_vector, sep = " ", collapse = "," ), split = ","))
  legend_vector_tools_user <- unlist(strsplit(paste(legend_vector_tools_user, "%", sep = "", collapse = ","), split = ","))
  
  sys_tools_percentage_vector <- format(round(prop.table(as.numeric(sys_mean_values_tools_vector))*100, digits = 2))
  legend_vector_tools_sys <- unlist(strsplit(paste(used_tools_unique, sys_tools_percentage_vector, sep = " ", collapse = "," ), split = ","))
  legend_vector_tools_sys <- unlist(strsplit(paste(legend_vector_tools_sys, "%", sep = "", collapse = ","), split = ","))
  
  real_replica_percentage_vector <- format(round(prop.table(as.numeric(real_mean_values_replica_vector))*100, digits = 2))
  legend_vector_replica_real <- unlist(strsplit(paste(used_replica_unique, real_replica_percentage_vector, sep = " ", collapse = "," ), split = ","))
  legend_vector_replica_real <- unlist(strsplit(paste(legend_vector_replica_real, "%", sep = "", collapse = ","), split = ","))
  
  user_replica_percentage_vector <- format(round(prop.table(as.numeric(user_mean_values_replica_vector))*100, digits = 2))
  legend_vector_replica_user <- unlist(strsplit(paste(used_replica_unique, user_replica_percentage_vector, sep = " ", collapse = "," ), split = ","))
  legend_vector_replica_user <- unlist(strsplit(paste(legend_vector_replica_user, "%", sep = "", collapse = ","), split = ","))
  
  sys_replica_percentage_vector <- format(round(prop.table(as.numeric(sys_mean_values_replica_vector))*100, digits = 2))
  legend_vector_replica_sys <- unlist(strsplit(paste(used_replica_unique, sys_replica_percentage_vector, sep = " ", collapse = "," ), split = ","))
  legend_vector_replica_sys <- unlist(strsplit(paste(legend_vector_replica_sys, "%", sep = "", collapse = ","), split = ","))
  
  
  pdf(summary_file_name, height = 12, width = 10)
  
  grob_list0 <- list(title, subtitle)
  grob_list1 <- list(t0, t1, t2, t3)
  
  grid_matrix0 <- rbind(c(1,1),
                        c(2,2))
  
  grid_matrix1 <- rbind(c(1, 1),
                        c(2, 2),
                        c(3, 4))
  
  g0 <- arrangeGrob(grobs = grob_list0, layout_matrix = grid_matrix0, heights = unit(c(-12, 14), c("cm", "cm")))
  g1 <- arrangeGrob(grobs = grob_list1, layout_matrix = grid_matrix1, heights = unit(c(0,17,0,14), c("cm", "cm", "cm", "cm")))

  grid.arrange(g0, g1, nrow = 2)
  
  
  par(mfrow = c(2,3))
  
  pie(as.numeric(real_mean_values_tools_vector), 
      labels=c("","","","","","","",""), font=2, 
      main="Averaged real runtimes of all tools",
      clockwise = FALSE,
      density = NULL,
      col = brewer.pal(length(used_tools_unique), "Set1"),
      border = NULL,
      lty = NULL,
      radius = 0.8)
  legend("bottom", legend=legend_vector_tools_real, cex=1.1, bty = "n", fill = brewer.pal(length(used_tools_unique), "Set1"))
  
  
  pie(as.numeric(user_mean_values_tools_vector), 
      labels=c("","","","","","","",""), font=2, 
      main="Averaged user runtimes of all tools",
      clockwise = FALSE,
      density = NULL,
      col = brewer.pal(length(used_tools_unique), "Set1"),
      border = NULL,
      lty = NULL,
      radius = 0.8)
  legend("bottom", legend=legend_vector_tools_user, cex=1.1, bty = "n", fill = brewer.pal(length(used_tools_unique), "Set1"))
  
 
  pie(as.numeric(sys_mean_values_tools_vector), 
      labels=c("","","","","","","",""), font=2, 
      main="Averaged system call times of all tools",
      clockwise = FALSE,
      density = NULL,
      col = brewer.pal(length(used_tools_unique), "Set1"),
      border = NULL,
      lty = NULL,
      radius = 0.8)
  legend("bottom", legend=legend_vector_tools_sys, cex=1.1, bty = "n", fill = brewer.pal(length(used_tools_unique), "Set1"))
  
  pie(as.numeric(real_mean_values_replica_vector), 
      labels=c("","",""), font=2, 
      main="Averaged real runtimes of all replica",
      clockwise = FALSE,
      density = NULL,
      col = brewer.pal(length(used_replica_unique), "Set1"),
      border = NULL,
      lty = NULL,
      radius = 0.8)
  legend("bottom", legend=legend_vector_replica_real, cex=1.1, bty = "n", fill = brewer.pal(length(used_replica_unique), "Set1"))
  
  pie(as.numeric(user_mean_values_replica_vector), 
      labels=c("","","","","","","",""), font=2, 
      main="Averaged user runtimes of all replica",
      clockwise = FALSE,
      density = NULL,
      col = brewer.pal(length(used_replica_unique), "Set1"),
      border = NULL,
      lty = NULL,
      radius = 0.8)
  legend("bottom", legend=legend_vector_replica_user, cex=1.1, bty = "n", fill = brewer.pal(length(used_replica_unique), "Set1"))
  
  pie(as.numeric(sys_mean_values_replica_vector), 
      labels=c("","","","","","","",""), font=2, 
      main="Averaged system call times of all replica",
      clockwise = FALSE,
      density = NULL,
      col = brewer.pal(length(used_replica_unique), "Set1"),
      border = NULL,
      lty = NULL,
      radius = 0.8)
  legend("bottom", legend=legend_vector_replica_sys, cex=1.1, bty = "n", fill = brewer.pal(length(used_replica_unique), "Set1"))
  invisible(dev.off())
}

if (scaling_flag == "scaling") {
  date <- format(Sys.time(), "%Y-%m-%d_%H:%M")
  scaling_output_filename <- paste("scaling_plot_", date, ".pdf", sep = "")
  
  pdf(scaling_output_filename)
  
  par(mfrow = c(2,2))
  
  for (i in 1:scaling_number_of_used_tools) {
  y_value_vector <- c()
    for (j in 1:length(scaling_cores_vector)) {
      y_value_vector <- c(y_value_vector, as.numeric(scaling_mean_real_times_vector[(i * j)]))
    }    

  plot(scaling_cores_vector, 
       y_value_vector,
       main = paste("Scaling behaviour of averaged real times \n for", used_tools_unique[i], sep = " "),
       col = brewer.pal(length(scaling_cores_vector), "Set1"), 
       xlab = "Number of used CPU cores",
       ylab = "Wall clock time in seconds",
       xaxt = "n",
       cex.main = 0.8)
    
  lines(scaling_cores_vector,
        y_value_vector)
  
 axis(side = 1, 
       at = scaling_cores_vector, 
       labels = scaling_cores_vector,
       tck=-.02)
  }
  invisible(dev.off())
}

options(warn = old_warnvalue)
