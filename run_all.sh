cp -r /media/user/Elements/scrnaseq_on_server /home/user/Documents/Maiolino_scrna

echo "Filed copied"

docker load -i /home/user/Documents/Maiolino_scrna/scrnaseq_on_server/scrnaseq_on_server_doc.tar

timepoints=("23days" "1month" "1.5month" "2month" "3month" "4month" "5month" "6month")

for tp in ${timepoints[@]}; do
    docker run -it --rm  
        -v /home/user/Documents/Maiolino_scrna/Data:/Data 
        -v /home/user/Documents/Maiolino_scrna/Results:/Results
        scrnaseq_on_server_doc Rscript RUN.R ${i}
done

wait

echo "Analysis finished"

cp -r /home/user/Documents/Maiolino_scrna/Results /media/user/Elements/Iterations_on_server_result

echo "Results copied"