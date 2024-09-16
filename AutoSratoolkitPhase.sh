#!/bin/bash

# Define input file and directories
input_file="/home/gevazatorsky/Documents/AlonProjects/EranSegalAntbiotics/RowFilesFolder/err_WGS.csv"
output_dir="/home/gevazatorsky/Documents/AlonProjects/EranSegalAntbiotics/RowFilesFolder"
temp_dir="/media/gevazatorsky/storage2/Alon_Projects/Eran_Antib_ERRs_after_PhaseV/temp_RowFilesFolder"
results_dir="/media/gevazatorsky/storage2/Alon_Projects/Eran_Antib_ERRs_after_PhaseV"
err_done="/home/gevazatorsky/Documents/AlonProjects/EranSegalAntbiotics/RowFilesFolder/err_WGS_done.csv"
err_2read="/home/gevazatorsky/Documents/AlonProjects/EranSegalAntbiotics/RowFilesFolder/err_WGS_need_to_read.csv"

# Create temp_dir if it doesn't exist
#mkdir -p "$temp_dir"

# Read done errors into an array
declare -A err_done_set
while IFS= read -r accession; do
    err_done_set["$accession"]=1
done < "$err_done"

# Initialize loop counter
loop_count=0

# Loop through each accession in err_2read
while IFS= read -r accession; do
    # Check if accession is already processed
    if [[ -z ${err_done_set["$accession"]} ]]; then
        (( loop_count++ ))
        echo "Processing accession: $accession #$loop_count"

        # Download fastq files using fasterq-dump
        echo "Downloading $accession"
        fastq-dump --split-files --outdir "$temp_dir" "$accession"

        # Rename downloaded files to p1 and p2
        mv "${temp_dir}/${accession}_1.fastq" "${temp_dir}/p1.fastq"
        mv "${temp_dir}/${accession}_2.fastq" "${temp_dir}/p2.fastq"
        
        # Run PhaseFinder
        echo "Running PhaseFinder1"
        python /home/gevazatorsky/PhaseFinder/PhaseFinder.py ratio -i "${temp_dir}/bact1.ID.fasta" -1 "${temp_dir}/p1.fastq" -2 "${temp_dir}/p2.fastq" -p 16 -o "${results_dir}/out_${accession}_1"
        echo "Running PhaseFinder2"
        python /home/gevazatorsky/PhaseFinder/PhaseFinder.py ratio -i "${temp_dir}/bact2.ID.fasta" -1 "${temp_dir}/p1.fastq" -2 "${temp_dir}/p2.fastq" -p 16 -o "${results_dir}/out_${accession}_2"

        # Clean up temp_dir
        echo "Cleaning up temp directory"
        rm -f "${temp_dir}/p1.fastq" "${temp_dir}/p2.fastq"
    else
        echo "Skipping accession $accession as it is already processed."
    fi
done < "$err_2read"

