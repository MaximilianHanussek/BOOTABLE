#!/bin/bash

# This script will collect all relevant output and results file from BOOTABLE in a .tar archive with the current time stamp 
dt=$(date '+%Y-%m-%d_%H-%M')
tar -cf BOOTABLE_results_collection_$dt.tar benchmark_summary_* bootable_system_info.txt scaling_plot* nmon_stats/*
