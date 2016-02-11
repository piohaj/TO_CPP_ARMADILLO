INCLUDES=-I/usr/include -I/usr/local/include
LIBS=-L/usr/local/lib 
LDFLAGS=-larmadillo -ltbb -llapack

vf3: my_vectfit.cpp vf3_main.cpp
	icpc -g -O3 -xAVX -DICC -debug inline-debug-info  $(INCLUDES) $(LIBS) vf3_main.cpp my_vectfit.cpp -ovf3 $(LDFLAGS)

