INCLUDES=-I/usr/include -I/usr/local/include -I/opt/intel/mkl/include
LIBS=-L/usr/local/lib -L/opt/intel/lib/intel64 -L/opt/intel/mkl/lib/intel64 -L/opt/intel/compilers_and_libraries_2016.1.150/linux/compiler/lib/intel64_lin
LDFLAGS=-larmadillo  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -ltbb

vf3: my_vectfit.cpp vf3_main.cpp network_model.cpp
	icpc -g0 -Ofast $(INCLUDES) $(LIBS) vf3_main.cpp my_vectfit.cpp -ovf3 $(LDFLAGS)

