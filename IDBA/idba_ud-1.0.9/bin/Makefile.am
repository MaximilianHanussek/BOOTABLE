AM_LDFLAGS = -fopenmp -pthread
LIBS = $(top_srcdir)/lib/libassembly.a
AM_CXXFLAGS = -Wall -O3 -fopenmp -pthread
INCLUDES =  -I$(top_srcdir)/src -I$(top_srcdir)/gtest_src 

parallel_blat_SOURCES = $(top_srcdir)/src/release/parallel_blat.cpp
sim_reads_SOURCES = $(top_srcdir)/src/release/sim_reads.cpp
fq2fa_SOURCES = $(top_srcdir)/src/release/fq2fa.cpp
idba_SOURCES = $(top_srcdir)/src/release/idba.cpp
raw_n50_SOURCES = $(top_srcdir)/src/release/raw_n50.cpp
idba_ud_SOURCES = $(top_srcdir)/src/release/idba_ud.cpp
validate_contigs_blat_SOURCES = $(top_srcdir)/src/release/validate_contigs_blat.cpp
split_scaffold_SOURCES = $(top_srcdir)/src/release/split_scaffold.cpp
filterfa_SOURCES = $(top_srcdir)/src/release/filterfa.cpp 

bin_PROGRAMS = \
	idba 

noinst_PROGRAMS = \
	parallel_blat \
	sim_reads \
	fq2fa \
	raw_n50 \
	idba_ud \
	validate_contigs_blat \
	split_scaffold \
	filterfa 

