#!/bin/bash

# echo "Copying folder from Drive"

# cp -r /media/"$USER"/Elements/scrnaseq_on_server "$HOME"/Maiolino_scrna

# echo "File copied - Loading Docker image"

# docker load -i "$HOME"/Maiolino_scrna/scrnaseq_on_server/scrnaseq_on_server_doc.tar

# echo "Docker image loaded"
echo "Starting analyses"

# Start timing
start_time=$(date +%s)

timepoints=("23days" "1month" "1.5month" "2month" "3month" "4month" "5month" "6month")

for tp in "${timepoints[@]}"; do
    docker run -it --rm \
        -v "$HOME"/Maiolino_scrna/Data:/Data \
        -v "$HOME"/Maiolino_scrna/Results:/Results \
        scrnaseq_on_server_doc Rscript scripts/RUN.R "$tp"
done

wait

# End timing
end_time=$(date +%s)
runtime=$((end_time - start_time))
minutes=$((runtime / 60))
seconds=$((runtime % 60))

echo "Analysis finished"

#echo "Transferring results"

# cp -r "$HOME"/Maiolino_scrna/Results /media/"$USER"/Elements/Iterations_on_server_result

#echo "Results transferred - END"

echo "Total runtime: ${minutes} minutes and ${seconds} seconds"
echo "Total runtime: ${minutes} minutes and ${seconds} seconds" > "$HOME"/Maiolino_scrna/Results/runtime.txt
