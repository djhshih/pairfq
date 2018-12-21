CXX=g++
CPPFLAGS=-O3
DESTDIR ?= /usr/local

all: pairfq
	

pairfq: pairfq.cpp
	$(CXX) $(CPPFLAGS) $? -o $@

check: pairfq.cpp
	mkdir -p out
	$(CXX) -coverage -O0 $? -o check
	./check || true
	./check data/r1.fq data/r2.fq out/out.fq natural null
	diff data/out.fq out/out.fq
	./check data/r1.fq data/r2.fq out/out.fq natural out/unpaired.fq
	diff data/out.fq out/out.fq
	diff data/unpaired.fq out/unpaired.fq

coverage: check
	gcov pairfq.cpp

test: check
	

install: pairfq
	mkdir -p $(DESTDIR)/bin/
	install pairfq $(DESTDIR)/bin/

clean:
	rm -f pairfq check
	rm -f *.exe *.gcov *.gcno *.gcda
	rm -rf out

