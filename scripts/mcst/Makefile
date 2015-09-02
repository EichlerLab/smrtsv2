INCLUDE=../../dist/htslib ../../dist/samtools
GCC_INCLUDE=$(addprefix -I,$(INCLUDE))
LIBRARY=../../dist/htslib ../../dist/samtools
GCC_LIBRARY=$(addprefix -L,$(LIBRARY))
BLASR=../../dist/blasr/common
LINKS=-lbam -lm -lz -lpthread
BINARIES=screenInversions maxq bedq bamacc hardstop splitread coverage

all: $(BINARIES)

screenInversions: ScreenInversions.cpp InversionAlign.h
	g++ -O3 $< -o $@ -I $(BLASR) -lpthread

coverage: PrintGenomeCoverage.cpp ../../dist/htslib/libhts.a
	g++ -g $^ $(GCC_INCLUDE) $(GCC_LIBRARY) $(LINKS) -o $@

maxq: MaxQ.cpp ../../dist/htslib/libhts.a
	g++ -O3 $^ $(GCC_INCLUDE) $(GCC_LIBRARY) $(LINKS) -o $@

bedq: MaxQBed.cpp ../../dist/htslib/libhts.a
	g++ -O3 $^ $(GCC_INCLUDE) $(GCC_LIBRARY) $(LINKS) -o $@

bamc: PrintCoverage.cpp ../../dist/htslib/libhts.a
	g++ -O3 $^ $(GCC_INCLUDE) $(GCC_LIBRARY) $(LINKS) -o $@

hardstop: PrintHardStop.cpp ../../dist/htslib/libhts.a
	g++ -g -O3 $^ $(GCC_INCLUDE) $(GCC_LIBRARY) $(LINKS) -o $@

splitread: PrintSplitReads.cpp ../../dist/htslib/libhts.a
	g++ -g -O3 $^ $(GCC_INCLUDE) $(GCC_LIBRARY) $(LINKS) -o $@

bamacc: ReadStats.cpp ../../dist/htslib/libhts.a
	g++ -O3 $^ $(GCC_INCLUDE) $(GCC_LIBRARY) $(LINKS) -o $@

clean:
	rm -f $(BINARIES)
