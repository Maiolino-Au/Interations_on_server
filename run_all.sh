#!/bin/bash

# Automatically copies the file from the external drives and loades the docker image
# echo "Copying folder from Drive"
# cp -r /media/"$USER"/Elements/scrnaseq_on_server "$HOME"/Maiolino_scrna
# echo "File copied - Loading Docker image"
# docker load -i "$HOME"/Maiolino_scrna/scrnaseq_on_server/scrnaseq_on_server_doc.tar
# echo "Docker image loaded"

# Enshure that the Results folder exists
mkdir -p "$HOME"/Maiolino_scrna/Results

echo "Starting analyses"

# Start timing
start_time=$(date +%s)

# Define the 8 timepoints that need to be analyzed
timepoints=("23days" "1month" "1.5month" "2month" "3month" "4month" "5month" "6month")

# Starts 8 parallel docker containers from the same image, one for each timepoint
for tp in "${timepoints[@]}"; do
    docker run -it --rm \
        -v "$HOME"/Maiolino_scrna/Data:/Data \
        -v "$HOME"/Maiolino_scrna/Results:/Results \
        scrnaseq_on_server_doc Rscript scripts/RUN.R "$tp"
done
wait # Waits for every docker to finish

# End timing
end_time=$(date +%s)
runtime=$((end_time - start_time))
hours=$((runtime / 3600))
minutes=$((runtime / 60 - hours * 60))
seconds=$((runtime % 60))

echo "Analysis finished"

# Copy automatically the results on the external drive
# echo "Transferring results"
# cp -r "$HOME"/Maiolino_scrna/Results /media/"$USER"/Elements/Iterations_on_server_result
# echo "Results transferred - END"

# Print total runtime and saves it in /Results
echo "Total runtime: ${hours} hours, ${minutes} minutes and ${seconds} seconds" | tee "$HOME"/Maiolino_scrna/Results/runtime.txt

