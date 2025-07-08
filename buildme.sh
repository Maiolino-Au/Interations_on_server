docker build -t scrnaseq_on_server_doc .
echo "image build"
docker save -o scrnaseq_on_server_doc.tar scrnaseq_on_server_doc
echo "image saved"
cp ./scrnaseq_on_server_doc.tar /media/user/Elements/scrnaseq_on_server/
echo "image transfered"
cp ./run_all.sh /media/user/Elements/scrnaseq_on_server/
echo "file transfered"