Muon lifetime
=======================

To compile, do the following:  
$ ./autogen.sh  
$ ./configure --prefix=\`pwd\`  
$ make  
$ make install  

To run everything in the assignment, run the muonlifetime script in
$PREFIX/bin:  
$ ./bin/muonlifetime  
To suppress the display of plots, pass the '-q' flag to muonlifetime.  
$ ./bin/muonlifetime -q  

Warning: Some of the calibration steps use d\_type from the dirent struct.
Because of that, this code can only be run on a filesystem that supports
d\_type. For more details, see the man page for readdir.
