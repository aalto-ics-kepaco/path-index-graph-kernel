CC = g++

# Use 64-bit integers in a 64-bit environment.
SIZE_FLAGS = -DMASSIVE_DATA_RLCSA

# Parallelism is supported by either libstdc++ Parallel Mode or MCSTL.
# PARALLEL_FLAGS = -DMULTITHREAD_SUPPORT -D_GLIBCXX_PARALLEL -fopenmp
# MCSTL_ROOT = /fs-3/d/jltsiren/suds/mcstl
# PARALLEL_FLAGS = -DMULTITHREAD_SUPPORT -I$(MCSTL_ROOT)/c++ -fopenmp

# Vectors using nibble codes instead of delta codes are faster, but they also
# take up more space.
# VECTOR_FLAGS = -DUSE_NIBBLE_VECTORS

#CPPFLAGS = -Wall -g $(SIZE_FLAGS) $(PARALLEL_FLAGS) $(VECTOR_FLAGS)
CPPFLAGS = -Wall -O3 $(SIZE_FLAGS) $(PARALLEL_FLAGS) $(VECTOR_FLAGS)
OBJS = rlcsa.o rlcsa_builder.o sasamples.o lcpsamples.o bits/array.o bits/bitvector.o bits/deltavector.o bits/rlevector.o bits/nibblevector.o misc/parameters.o misc/utils.o qsufsort/qsufsort.o

VPATH = bits:misc:utils


default: rlcsa.a


rlcsa.a: $(OBJS)
	ar rcs rlcsa.a $(OBJS)

depend:
	g++ -MM *.cpp bits/*.cpp misc/*.cpp qsufsort/*.c utils/*.cpp > dependencies.mk

rlcsa_test: rlcsa_test.o rlcsa.a
	$(CC) $(CPPFLAGS) -o rlcsa_test rlcsa_test.o rlcsa.a

lcp_test: lcp_test.o rlcsa.a
	$(CC) $(CPPFLAGS) -o lcp_test lcp_test.o rlcsa.a

parallel_build: parallel_build.o rlcsa.a
	$(CC) $(CPPFLAGS) -o parallel_build parallel_build.o rlcsa.a

build_rlcsa: build_rlcsa.o rlcsa.a
	$(CC) $(CPPFLAGS) -o build_rlcsa build_rlcsa.o rlcsa.a

locate_test: locate_test.o rlcsa.a
	$(CC) $(CPPFLAGS) -o locate_test locate_test.o rlcsa.a

display_test: display_test.o rlcsa.a
	$(CC) $(CPPFLAGS) -o display_test display_test.o rlcsa.a

read_bwt: read_bwt.o rlcsa.a
	$(CC) $(CPPFLAGS) -o read_bwt read_bwt.o rlcsa.a

extract_sequence: extract_sequence.o rlcsa.a
	$(CC) $(CPPFLAGS) -o extract_sequence extract_sequence.o rlcsa.a

rlcsa_grep: rlcsa_grep.o rlcsa.a
	$(CC) $(CPPFLAGS) -o rlcsa_grep rlcsa_grep.o rlcsa.a

build_plcp: build_plcp.o rlcsa.a
	$(CC) $(CPPFLAGS) -o build_plcp build_plcp.o rlcsa.a

sample_lcp: sample_lcp.o rlcsa.a
	$(CC) $(CPPFLAGS) -o sample_lcp sample_lcp.o rlcsa.a

extract_text: extract_text.o
	$(CC) $(CPPFLAGS) -o utils/extract_text extract_text.o

convert_patterns: convert_patterns.o
	$(CC) $(CPPFLAGS) -o utils/convert_patterns convert_patterns.o misc/utils.o

split_text: split_text.o
	$(CC) $(CPPFLAGS) -o utils/split_text split_text.o misc/utils.o

genpatterns: genpatterns.c
	gcc -O3 -Wall -o utils/genpatterns utils/genpatterns.c

clean:
	rm -f rlcsa.a rlcsa_test lcp_test parallel_build
	rm -f build_rlcsa locate_test display_test read_bwt extract_sequence rlcsa_grep build_plcp sample_lcp
	rm -f utils/extract_text utils/convert_patterns utils/split_text
	rm -f *.o bits/*.o misc/*.o qsufsort/*.o utils/*.o

package:
	mkdir rlcsa
	mkdir rlcsa/bits rlcsa/misc rlcsa/qsufsort rlcsa/utils
	cp LICENSE Makefile README dependencies.mk *.cpp *.h rlcsa
	cp bits/*.cpp bits/*.h rlcsa/bits
	cp misc/*.cpp misc/*.h rlcsa/misc
	cp qsufsort/*.c qsufsort/*.h rlcsa/qsufsort
	cp utils/*.cpp utils/*.py rlcsa/utils
	tar cfz rlcsa.tgz rlcsa
	rm -r rlcsa/*
	rmdir rlcsa

include dependencies.mk
