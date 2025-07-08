docker build -t scrnaseq_on_server_doc .
docker save -o scrnaseq_on_server_doc.tar scrnaseq_on_server_doc
cp ./scrnaseq_on_server_doc.tar /media/user/Elements/scrnaseq_on_server/
cp ./run_all.sh /media/user/Elements/scrnaseq_on_server/