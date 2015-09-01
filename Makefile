CXXFLAGS=-std=c++11 -g -Wall -Wno-unused-local-typedefs -I.
ifeq ($(RELEASE),1)
  CXXFLAGS+=-O2 -DNDEBUG
else
  CXXFLAGS+=-O1
endif
CXXFLAGS+=$(EXTRA_CXXFLAGS)

test: all
	./decay_chain_dp
	./decay_chain_cpp
	./decay_chain_gmp
	# lu decomposition fails with mpfr:
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
