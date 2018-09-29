for file in `ls *.pdf`
do
    pdfcrop --margin 0 $file $file
done
