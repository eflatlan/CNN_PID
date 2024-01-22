#!/bin/bash

# The directory where you downloaded the Google Drive folder
directory="/root/AllH5_to_fix/"

# Check if gdown is installed
if ! command -v gdown &> /dev/null; then
    echo "gdown is not installed. Attempting to install..."
    pip install gdown
fi

# Check if parallel is installed
if ! command -v parallel &> /dev/null; then
    echo "parallel is not installed. Attempting to install..."
    sudo yum install parallel
fi

# Define a function to process a file
process_file() {
    file=$1
    # Get the file name
    file_name=$(basename "$file")

    # Check if file already exists in the current directory and rename if necessary
    counter=1
    new_file_name=$file_name
    while [ -f "$new_file_name" ]; do
        new_file_name="${file_name%.*}_$counter.${file_name##*.}"
        ((counter++))
    done
    if [ "$new_file_name" != "$file_name" ]; then
        cp "$file" "$new_file_name"
        echo "File $file_name already exists. Renamed to $new_file_name"
        file_name=$new_file_name
    fi

    # Print the file name
    echo "Processing file: $file_name"

    # Process the file with your C++ macro
    root -l -b -q 'UpdateH5.C("'$file_name'")'
}

# Export the function so it's available to parallel
export -f process_file

# Specify the number of parallel jobs
num_jobs=5

# Find all .h5 files in the directory and its subdirectories, and process them in parallel
find "$directory" -name "*.h5" -type f | parallel -j "$num_jobs" process_file

