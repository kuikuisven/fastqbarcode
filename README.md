# fastqbarcode
Find most frequent barcodes in fastq files and dump records in separate fastq files. 

This tool is VERY fast, as it's written in C.

Note that this should be run on a computer that can contain the whole file in RAM.



Usage: fastqbarcode -f FASTQ.

   -- get data from pipe.

   -p prefix for output files

   -m number of barcode output (default 10)

   -h Displays this help message.

   -d Verbose



A version that doesn't require large amouts of RAM is in developpment.
