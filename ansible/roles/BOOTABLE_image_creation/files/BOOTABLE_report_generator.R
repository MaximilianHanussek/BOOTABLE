# Author: Maximilian Hanussek, maximilian.hanussek@uni-tuebingen.de
# Parse and generate plots of BOOTABLE

# Choose CRAN mirror in beforehand to run the script without interaction
chooseCRANmirror(ind = 33)

# Save the old warnvalue
old_warnvalue <- getOption("warn")

# Suppress warnings on stdout for a cleaner user experience
options(warn = -1)

# Load or install package if necessary
load.fun <- function(x) {  # x is the package name
  x <- as.character(substitute(x))
  if(isTRUE(x %in% .packages(all.available=TRUE))) {  # Check if package is already installed
    eval(parse(text=paste("require(", x, ")", sep="")))  # If yes, load it
  } else {  # Else install it and load it
    update.packages() # Recommended before installing so that dependencies are the latest version
    eval(parse(text=paste("install.packages('", x, "')", sep=""))) # Install package
    eval(parse(text=paste("require(", x, ")", sep="")))  # Load package
  }
} 

load.fun(grid)          # Load grid package, required for tables
load.fun(gridExtra)     # Load gridExtra package, required for tables
load.fun(RColorBrewer)  # Load RColorBrewer package, required for pie plots
load.fun(stringr)       # Load stringr package, required for host information parsing

# Set current working directory to the directory where the script is executed
workingdir <- "./"

# Get argument from the command line, whether scaling plot should be generated or not
args <- commandArgs(trailingOnly = TRUE)
argslen <- length(args)
if (argslen == 1){          # Handle error if flag is not set
  scaling_flag <- args[1]
} else {                    # Set default value FALSE
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

# Check if Layer3 cache is available otherwise use Layer2 cache
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

# Create dataframe with host information for printing it in tabular form
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


# Get all summary text files produced from BOOTABLE, containing the time output
summary_file_paths <- list.files(path = workingdir, pattern = "^benchmark_summary_.*\\.txt", full.names = TRUE)  # Get the full paths of all report
summary_file_names <- list.files(path = workingdir, pattern = "^benchmark_summary_.*\\.txt", full.names = FALSE) # Get only the report names
ordered_summary_file_paths <- summary_file_paths[order(as.numeric(gsub("[^\\d]+", "\\1", summary_file_paths, perl = TRUE)))] # Full paths ordered by integer of the used CPU cores
number_of_summary_files <- length(summary_file_paths)  # Get the number of created summary text files
names_time_vector <- c("real", "user", "sys")  # Create a name vector for the 3 different measured times

# Initialize scaling capabilities vectors
scaling_cores_vector <- c()             # Initialize vector for the used cores numbers (e.g. 1,7,14,28)
scaling_mean_real_times_vector <- c()   # Initialize vector for the mean values of the real times of all replica
scaling_number_of_used_tools <- c()     # Initialize vector for the number of tools used in the benchmark run


###########################################################
### START Iterate over all available summary text files ###
###########################################################
for (summary_file in ordered_summary_file_paths){                                # Iterate over the ordered filepaths vector
  used_tools <- c()                                                              # Initialize used tools name vector   
  used_replica <- c()                                                            # Initialize used replica value vector  
  real_values_all_vector <- c()                                                  # Initialize real time values vector
  user_values_all_vector <- c()                                                  # Initialize user time values vector
  sys_values_all_vector <- c()                                                   # Initialize sys time values vector
  summary_name <- summary_file_names[which(summary_file == summary_file_paths)]  # Get the index of the current summary file in the vector
  summary_name <- strsplit(summary_name, "\\.")[[1]][1]                          # Get rid of the .txt extension
  summary_file_name <- paste(summary_name, ".pdf", sep = "")                     # Create complete report filename 
  
  summary_file_lines <- readLines(con = summary_file)                            # Read in current summary file, linewise
  replica_entries_line_numbers <- grep("Replica", summary_file_lines)            # Get onyl the replica line numbers
  real_entries_line_numbers <- grep ("real ", summary_file_lines)                # Get onyl the real time values line numbers
  user_entries_line_numbers <- grep ("user ", summary_file_lines)                # Get only the user time values line numbers
  sys_entries_line_numbers <- grep ("sys ", summary_file_lines)                  # Get onyl the sys time values line numbers
  
  number_of_cores <- as.integer(strsplit(summary_file_lines[1], "\\s+")[[1]][4]) # Get the number of cores used in this summary file and this benchmark run
  scaling_cores_vector <- c(scaling_cores_vector, number_of_cores)               # Append the number of used cpu cores for the scalability plot later

  #######################################################################################
  ### START Iterate over all replica entries to get the tool names and replica values ###
  #######################################################################################
  for (replica_entry in replica_entries_line_numbers){
    used_replica <- c(used_replica, strsplit(summary_file_lines[replica_entry], "\\s+")[[1]][1]) # Get the used replica string "Replica_N"
    used_tools <- c(used_tools, strsplit(summary_file_lines[replica_entry], "\\s+")[[1]][2])     # Get used tool names
  }
  #####################################################################################
  ### END Iterate over all replica entries to get the tool names and replica values ###
  #####################################################################################
  
  used_tools_unique <- unique(used_tools)                                                                 # Remove duplciates of extracted tool names
  number_of_used_tools <- length(used_tools_unique)                                                       # Get the number of used tools
  scaling_number_of_used_tools <- number_of_used_tools                                                    # Duplicate vector above for scaling functionality
  used_replica_unique <- unique(used_replica)                                                             # Remove duplicate of extracted replicate strings
  number_of_used_replica <- length(used_replica_unique)                                                   # Get the number of used replicas
  
  df_all <- data.frame(matrix(data = list(), ncol = number_of_used_tools, nrow = number_of_used_replica)) # Initialize data frame for all colected data
  df_all_printable <- data.frame(matrix(ncol = number_of_used_tools, nrow = number_of_used_replica))      # Initialize data frame in printable version
  colnames(df_all) <- used_tools_unique                                                                   # Create column names (Tool)
  rownames(df_all) <- used_replica_unique                                                                 # Create row names (Replica)
  colnames(df_all_printable) <- used_tools_unique                                                         # Create column names (Tool) for printabel version
  rownames(df_all_printable) <- used_replica_unique                                                       # Create row names (Replica) for printabel version
  
  i <- 1  # Initialize inner-loop variable
  j <- 1  # Initialize inner-loop variable
  
  #############################################################################################
  ### START Iterate over all time values and add them to the empty dataframe for all values ###
  #############################################################################################
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
  ###########################################################################################
  ### END Iterate over all time values and add them to the empty dataframe for all values ###
  ###########################################################################################
  
  used_time_values <- c("average real", "average user", "average sys")      # Create average name vector      
  df_toolwise <- data.frame(matrix(ncol = number_of_used_tools, nrow = 3))  # Create empty data frame for table in the report
  colnames(df_toolwise) <- used_tools_unique                                # Column names are the tool names
  rownames(df_toolwise) <- used_time_values                                 # Row names are the time values (average)
  real_mean_values_tools_vector <- c()                                      # Initialize real average values vector
  user_mean_values_tools_vector <- c()                                      # Initialize user average values vector
  sys_mean_values_tools_vector <- c()                                       # Initialize sys average values vector
  
  ##############################################################
  ### START Iterate over all tools and fill empty data frame ###
  ##############################################################
  for (tool in 1:number_of_used_tools) {
    real_values_vector <- c()            # Initialize real time values vector
    user_values_vector <- c()            # Initialize real user values vector
    sys_values_vector  <- c()            # Initialize real sys values vector
    
    #########################################################################################
    ### START Iterate over all replicas for every tool to create the mean (average) value ###
    #########################################################################################
    for (value in 1:number_of_used_replica){
      real_values_vector <- c(real_values_vector, df_all[[tool]][[value]][1]) # Get all real values from every replica
      user_values_vector <- c(user_values_vector, df_all[[tool]][[value]][2]) # Get all user values from every replica
      sys_values_vector <- c(sys_values_vector, df_all[[tool]][[value]][3])   # Get all sys values from every replica
    }
    #######################################################################################
    ### END Iterate over all replicas for every tool to create the mean (average) value ###
    #######################################################################################
    
    df_toolwise[[1, tool]] <- format(round(mean(real_values_vector), digits = 2), nsmall = 2)  # Calculate mean of real values over replica with 2 digits after the comma
    df_toolwise[[2, tool]] <- format(round(mean(user_values_vector), digits = 2), nsmall = 2)  # Calculate mean of user values over replica with 2 digits after the comma
    df_toolwise[[3, tool]] <- format(round(mean(sys_values_vector), digits = 2), nsmall = 2)   # Calculate mean of sys values over replica with 2 digits after the comma
    
    real_mean_values_tools_vector <- c(real_mean_values_tools_vector, df_toolwise[[1, tool]])  # Collect all real mean time values for every tool
    user_mean_values_tools_vector <- c(user_mean_values_tools_vector, df_toolwise[[2, tool]])  # Collect all user mean time values for every tool
    sys_mean_values_tools_vector <- c(sys_mean_values_tools_vector, df_toolwise[[3, tool]])    # Collect all sys mean time values for every tool
    
    real_values_all_vector <- c(real_values_all_vector, real_values_vector)                    # Collect all real time values from every replica
    user_values_all_vector <- c(user_values_all_vector, user_values_vector)                    # Collect all user time values from every replica
    sys_values_all_vector  <- c(sys_values_all_vector, sys_values_vector)                      # Collect all sys time values from every replica
  }
  ############################################################
  ### END Iterate over all tools and fill empty data frame ###
  ############################################################
  
  scaling_mean_real_times_vector <- c(scaling_mean_real_times_vector, real_mean_values_tools_vector) # Collect only the real time mean values from every summary txt file for scaling functionality
  
  df_replicawise <- data.frame(matrix(ncol = number_of_used_replica, nrow = 3))  # Create empty dataframe for table in report (average all tool times for one replica)
  colnames(df_replicawise) <- used_replica_unique                                # Column names are the replica (Replica_1, Replica_2, ...)
  rownames(df_replicawise) <- used_time_values                                   # Row names are the time values (average real, average user, average sys)
  
  real_mean_values_replica_vector <- c()                                         # Initialize real time average per replica vector
  user_mean_values_replica_vector <- c()                                         # Initialize user time average per replica vector
  sys_mean_values_replica_vector <- c()                                          # Initialize sys time average per replica vector
  
  #################################################################
  ### START Iterate over all replicas and fill empty data frame ###
  #################################################################
  for (replica in 1:number_of_used_replica) {
    real_values_vector <- c()                                                    # Initialize real value vector
    user_values_vector <- c()                                                    # Initialize user value vector
    sys_values_vector  <- c()                                                    # Initialize sys value vector
    
    #########################################################################################
    ### START Iterate over all tools for every replica to create the mean (average) value ###
    #########################################################################################
    for (value in 1:number_of_used_tools){
      real_values_vector <- c(real_values_vector, df_all[[value]][[replica]][1]) # Collect real values from every replica
      user_values_vector <- c(user_values_vector, df_all[[value]][[replica]][2]) # Collect user values from every replica
      sys_values_vector <- c(sys_values_vector, df_all[[value]][[replica]][3])   # Collect sys values from every replica
    }
    #######################################################################################
    ### END Iterate over all tools for every replica to create the mean (average) value ###
    #######################################################################################
    
    df_replicawise[[1, replica]] <- format(round(mean(real_values_vector), digits = 2), nsmall = 2)  # Calculate mean of real values over tools with 2 digits after the comma
    df_replicawise[[2, replica]] <- format(round(mean(user_values_vector), digits = 2), nsmall = 2)  # Calculate mean of user values over replica with 2 digits after the comma
    df_replicawise[[3, replica]] <- format(round(mean(sys_values_vector), digits = 2), nsmall = 2)   # Calculate mean of sys values over replica with 2 digits after the comma
    
    real_mean_values_replica_vector <- c(real_mean_values_replica_vector, df_replicawise[[1, replica]])  # Collect all real mean time values for every replica
    user_mean_values_replica_vector <- c(user_mean_values_replica_vector, df_replicawise[[2, replica]])  # Collect all user mean time values for every replica
    sys_mean_values_replica_vector <- c(sys_mean_values_replica_vector, df_replicawise[[3, replica]])    # Collect all sys mean time values for every replica
  }
  ###############################################################
  ### END Iterate over all replicas and fill empty data frame ###
  ###############################################################
  
  total_time_real <- format(sum(real_values_all_vector), nsmall = 2)  # Sum of all real values from every tool for all replica in the loop
  total_time_user <- format(sum(user_values_all_vector), nsmall = 2)  # Sum of all user values from every tool for all replica in the loop    
  total_time_sys  <- format(sum(sys_values_all_vector), nsmall = 2)   # Sum of all sys values from every tool for all replica in the loop

  t0 <- tableGrob(df_system_information, rows = rownames(df_system_information), cols = colnames(df_system_information))  # Create graphical table for the systeminfos with row and col names
  t1 <- tableGrob(t(df_all_printable), rows = colnames(df_all_printable), cols = rownames(df_all_printable))              # Create graphical table for all measured times wit row and col names
  t2 <- tableGrob(t(df_toolwise), rows = colnames(df_toolwise), cols = rownames(df_toolwise))                             # Create graphical table with averaged values over all replicas per tool
  t3 <- tableGrob(t(df_replicawise), rows = colnames(df_replicawise), cols = rownames(df_replicawise))                    # Create graphical table with averaged values over all tools per replica
  
  
  title <- textGrob("Statistical Report of BOOTABLE", gp = gpar(fontsize=20,font=2), just = "center")   # Create over all title for the report
  subtitle <- textGrob(date(), gp = gpar(fontsize = 14, font = 1), just = "center")                     # Create subtitle with the date and time of creation
  #margin <- unit(0.5, "line")  # Not used so far
  
  # Collect pie chart values (percentage) and create legends
  # Toolwise over all replica
  real_tools_percentage_vector <- format(round(prop.table(as.numeric(real_mean_values_tools_vector))*100, digits = 2))                            # Get the percentage of every tool over all replica from the overall real time
  legend_vector_tools_real <- unlist(strsplit(paste(used_tools_unique, real_tools_percentage_vector, sep = " ", collapse = "," ), split = ","))   # Paste the percentage values and the tool names togethter in one vector
  legend_vector_tools_real <- unlist(strsplit(paste(legend_vector_tools_real, "%", sep = "", collapse = ","), split = ","))                       # Add '%' char to the percentage values
  
  user_tools_percentage_vector <- format(round(prop.table(as.numeric(user_mean_values_tools_vector))*100, digits = 2))                            # Get the percentage of every tool over all replica from the overall user time
  legend_vector_tools_user <- unlist(strsplit(paste(used_tools_unique, user_tools_percentage_vector, sep = " ", collapse = "," ), split = ","))   # Paste the percentage values and the tool names togethter in one vector
  legend_vector_tools_user <- unlist(strsplit(paste(legend_vector_tools_user, "%", sep = "", collapse = ","), split = ","))                       # Add '%' char to the percentage values
  
  sys_tools_percentage_vector <- format(round(prop.table(as.numeric(sys_mean_values_tools_vector))*100, digits = 2))                              # Get the percentage of every tool over all replica from the overall sys time
  legend_vector_tools_sys <- unlist(strsplit(paste(used_tools_unique, sys_tools_percentage_vector, sep = " ", collapse = "," ), split = ","))     # Paste the percentage values and the tool names togethter in one vector
  legend_vector_tools_sys <- unlist(strsplit(paste(legend_vector_tools_sys, "%", sep = "", collapse = ","), split = ","))                         # Add '%' char to the percentage values
  
  # Replicawise over all tools
  real_replica_percentage_vector <- format(round(prop.table(as.numeric(real_mean_values_replica_vector))*100, digits = 2))                              # Get the percentage of every replica over all tools from the overall real time
  legend_vector_replica_real <- unlist(strsplit(paste(used_replica_unique, real_replica_percentage_vector, sep = " ", collapse = "," ), split = ","))   # Paste the percentage values and the tool names togethter in one vector
  legend_vector_replica_real <- unlist(strsplit(paste(legend_vector_replica_real, "%", sep = "", collapse = ","), split = ","))                         # Add '%' char to the percentage values
  
  user_replica_percentage_vector <- format(round(prop.table(as.numeric(user_mean_values_replica_vector))*100, digits = 2))                              # Get the percentage of every replica over all tools from the overall user time
  legend_vector_replica_user <- unlist(strsplit(paste(used_replica_unique, user_replica_percentage_vector, sep = " ", collapse = "," ), split = ","))   # Paste the percentage values and the tool names togethter in one vector
  legend_vector_replica_user <- unlist(strsplit(paste(legend_vector_replica_user, "%", sep = "", collapse = ","), split = ","))                         # Add '%' char to the percentage values
  
  sys_replica_percentage_vector <- format(round(prop.table(as.numeric(sys_mean_values_replica_vector))*100, digits = 2))                                # Get the percentage of every replica over all tools from the overall real time
  legend_vector_replica_sys <- unlist(strsplit(paste(used_replica_unique, sys_replica_percentage_vector, sep = " ", collapse = "," ), split = ","))     # Paste the percentage values and the tool names togethter in one vector
  legend_vector_replica_sys <- unlist(strsplit(paste(legend_vector_replica_sys, "%", sep = "", collapse = ","), split = ","))                           # Add '%' char to the percentage values
  
  
  # Start creating the pdf document
  pdf(summary_file_name, height = 12, width = 10)
  
  grob_list0 <- list(title, subtitle) # Convert titles to list as grob expects lists of grob tables
  grob_list1 <- list(t0, t1, t2, t3)  # Convert tables to list as grob expects lists of grob tables
  
  
  # Define layout on the page with grid matrices
  # One column for the title and subtitle each
  grid_matrix0 <- rbind(c(1,1),
                        c(2,2))
  
  # One column for the system informations, 
  # one column for the overall values table,
  # two columns for the toolwise and replicawise table -> printed beneath each other
  grid_matrix1 <- rbind(c(1, 1),
                        c(2, 2),
                        c(3, 4))
  
  # Connect table lists with layout and define fitting heights 
  g0 <- arrangeGrob(grobs = grob_list0, layout_matrix = grid_matrix0, heights = unit(c(-12, 14), c("cm", "cm")))
  g1 <- arrangeGrob(grobs = grob_list1, layout_matrix = grid_matrix1, heights = unit(c(0,17,0,14), c("cm", "cm", "cm", "cm")))

  # Finalize table setting process on pdf page
  grid.arrange(g0, g1, nrow = 2)
  
  
  # Next page starting with the pie chart plots (6 per page)
  par(mfrow = c(2,3))
  
  # Percentage pie plots toolwise (real, user, sys)
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
  
  
  # Percentage pie plots replicawise (real, user, sys)
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
  
  # Shut down opened pdf device and hide print on stdout (invisible) for convenience reasons 
  invisible(dev.off()) 
}
#########################################################
### END Iterate over all available summary text files ###
#########################################################

# Starting scaling plot possibility
if (scaling_flag == "scaling") {                                              # Check if scaling is desired by scaling flag
  date <- format(Sys.time(), "%Y-%m-%d_%H:%M")                                # Save date in variable (2042-10-23_07:00) for output name
  scaling_output_filename <- paste("scaling_plot_", date, ".pdf", sep = "")   # Create final output name
  
  pdf(scaling_output_filename)                                                # Open pdf device with file output name
  
  par(mfrow = c(2,2))                                                         # Generate 4 plots per page
  
  ######################################################
  ### START Iterate over every tool used for scaling ###
  ######################################################
  for (i in 1:scaling_number_of_used_tools) {
  y_value_vector <- c()                                                                             # Initialize empty vector for y axis values
  index <- 0                                                                                        # Initialize index variable with 0
  
    #####################################################################################################################  
    ### START Iterate over mean real time values vector and pick out the ones corresponding to the correct used cores ###
    #####################################################################################################################
    for (core_number in scaling_cores_vector) {             
      y_value_vector <- c(y_value_vector, as.numeric(scaling_mean_real_times_vector[(i + index)]))
      index <- index + scaling_number_of_used_tools
    }
    ###################################################################################################################
    ### END Iterate over mean real time values vector and pick out the ones corresponding to the correct used cores ###
    ###################################################################################################################
  
  # Set up plot for every used number of cores
  plot(scaling_cores_vector, 
       y_value_vector,
       main = paste("Scaling behaviour of averaged real times \n for", used_tools_unique[i], sep = " "),
       col = brewer.pal(length(scaling_cores_vector), "Set1"), 
       xlab = "Number of used CPU cores",
       ylab = "Wall clock time in seconds",
       xaxt = "n",
       cex.main = 0.8,
       log = "yx")      # Logarithmic x and y axes to see a linear behavior
  
  # Plot line for tool in set up empty plot for the used cores    
  lines(scaling_cores_vector,
        y_value_vector)
  
  # Define own lables for x-axis
  axis(side = 1, 
       at = scaling_cores_vector, 
       labels = scaling_cores_vector,
       tck=-.02)
  }
  ####################################################
  ### END Iterate over every tool used for scaling ###
  ####################################################
  
  # Shut down opened pdf device and hide print on stdout (invisible) for convenience reasons
  invisible(dev.off())    
}

# Set the general warn level bvack to the old one in order to not disrupt other ongoing programming
options(warn = old_warnvalue)