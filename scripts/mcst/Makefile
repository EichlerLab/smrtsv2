INCLUDE=/usr/local/include
LIBRARY=/usr/local/lib
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/common
LINKS=-lbam -lm -lz -lpthread

all: screenInversions maxq bedq bamacc hardstop splitread coverage

screenInversions: ScreenInversions.cpp InversionAlign.h
	g++ -O3 $< -o $@ -I $(BLASR) -lpthread

coverage: PrintGenomeCoverage.cpp $(LIBRARY)/libhts.a
	g++ -g $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

maxq: MaxQ.cpp $(LIBRARY)/libhts.a
	g++ -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

bedq: MaxQBed.cpp $(LIBRARY)/libhts.a
	g++ -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

bamc: PrintCoverage.cpp $(LIBRARY)/libhts.a
	g++ -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

hardstop: PrintHardStop.cpp $(LIBRARY)/libhts.a
	g++ -g -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

splitread: PrintSplitReads.cpp $(LIBRARY)/libhts.a
	g++ -g -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

bamacc: ReadStats.cpp $(LIBRARY)/libhts.a
	g++ -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@


