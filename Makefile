CXXFLAGS=-std=c++11 -g -Wall -Wno-unused-local-typedefs -I.
ifeq ($(RELEASE),1)
  CXXFLAGS+=-O2 -DNDEBUG
else
  CXXFLAGS+=-O1
endif
CXXFLAGS+=$(EXTRA_CXX_FLAGS)

all: cda_mp cda_dp

.PHONY: all

cda_mp: coupled_decay_arbitrary.cpp
	$(CXX) $(CXXFLAGS) -DWITH_MULTIPRECISION -o $@ $^

cda_dp: coupled_decay_arbitrary.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^
