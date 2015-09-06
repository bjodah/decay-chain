CXXFLAGS=-std=c++11 -g -Wall -Wno-unused-local-typedefs -I.
ifeq ($(RELEASE),1)
  CXXFLAGS+=-O2 -DNDEBUG
else
  CXXFLAGS+=-O1
endif
CXXFLAGS+=$(EXTRA_CXXFLAGS)
        log10abstol = -12;
        log10reltol = -12;
        log10tend = 0;
        log10dx = -14;
        adaptive = 1;
        N = 27;
        p = 1;
        a = 27;

test: all
	for method in 0 1 2 3 ; do \
		./decay_chain_dp -12 -12 0 -14 1 27 1 27 $$method ; \
		./decay_chain_cpp -12 -12 0 -14 1 27 1 27 $$method ; \
		./decay_chain_gmp -12 -12 0 -14 1 27 1 27 $$method ; \
	done
	# lu decomposition fails with mpfr (proably a bug in either odeint, ublas or mpfr):
	#./decay_chain_mpfr

all: decay_chain_dp decay_chain_cpp decay_chain_mpfr decay_chain_gmp

.PHONY: all test

decay_chain_dp: decay_chain.cpp
	$(CXX) $(CXXFLAGS) -DVALUE_TYPE_IDX=0 -o $@ $^

decay_chain_cpp: decay_chain.cpp
	$(CXX) $(CXXFLAGS) -DVALUE_TYPE_IDX=1 -o $@ $^

decay_chain_mpfr: decay_chain.cpp
	$(CXX) $(CXXFLAGS) -DVALUE_TYPE_IDX=2 -o $@ $^ -lmpfr

decay_chain_gmp: decay_chain.cpp
	$(CXX) $(CXXFLAGS) -DVALUE_TYPE_IDX=3 -o $@ $^ -lgmp
