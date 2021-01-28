# Creates a set of shorter-named symbolic links from a set of KOMODO2 input files
# You can run this script from inside the input folder
# or copy the command to the command line and change the input path.
# If your file names are longer/shorter, change the string slice
# e.g. ${i:0:13}
for i in `ls *.txt`; do ln -s $i ${i:0:13}; done 