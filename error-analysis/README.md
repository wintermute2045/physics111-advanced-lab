EAX homework
=======================

To compile, do the following:  
$ ./autogen.sh  
$ ./configure --prefix=\`pwd\`  
$ make  
$ make install  

To run everything in the assignment, run the eaxhomework executable in
$PREFIX/bin:  
$ ./bin/eaxhomework  
To suppress the display of plots, pass the '-q' flag to eaxhomework.  
$ ./bin/eaxhomework -q  
To run only a specific problem (3, 4, or 5), specify it with the '-p' flag:  
$ ./bin/eaxhomework -p 3
