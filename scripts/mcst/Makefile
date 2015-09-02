INCLUDE=/usr/local/include
LIBRARY=/usr/local/lib
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/common
LINKS=-lbam -lm -lz -lpthread -lhts
BINARIES=screenInversions maxq bedq bamacc hardstop splitread coverage

all: $(BINARIES)

screenInversions: ScreenInversions.cpp InversionAlign.h
	g++ -O3 $< -o $@ -I $(BLASR) -lpthread

coverage: PrintGenomeCoverage.cpp
	g++ -g $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

maxq: MaxQ.cpp
	g++ -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

bedq: MaxQBed.cpp
	g++ -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

bamc: PrintCoverage.cpp
	g++ -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

hardstop: PrintHardStop.cpp
	g++ -g -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

splitread: PrintSplitReads.cpp
	g++ -g -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

bamacc: ReadStats.cpp
	g++ -O3 $^ -I $(INCLUDE) -L $(LIBRARY) $(LINKS) -o $@

clean:
	rm -f $(BINARIES)
