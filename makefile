INCLUDES=-I/usr/include -I/usr/local/include
LIBS=-L/usr/local/lib 
LDFLAGS=-larmadillo -llapack 

vf3: my_vectfit.cpp vf3_main.cpp
	icpc -g0 -Ofast $(INCLUDES) $(LIBS) vf3_main.cpp my_vectfit.cpp -ovf3 $(LDFLAGS)

