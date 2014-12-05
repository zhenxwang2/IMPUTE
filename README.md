=======
///////    Project of Biostatistics 830  /////////////

C++ code to do imputation of hyplotypes with a HMM...

compile cpp code :

g++ -I ./libStatGen/include/ -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS main.cpp ../libStatGen/libStatGen.a -lm -lz -o impute

