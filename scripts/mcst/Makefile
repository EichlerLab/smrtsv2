SAMTOOLS=/net/eichler/vol5/home/mchaisso/software/samtools-1.2
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/common

all: screenInversions maxq bedq bamacc hardstop splitread coverage

screenInversions: ScreenInversions.cpp InversionAlign.h
	g++ -O3 $< -o $@ -I $(BLASR) -lpthread

coverage: PrintGenomeCoverage.cpp
	g++ -g  $^ -I $(SAMTOOLS)  -L $(SAMTOOLS)  -lbam -lm -lz -o $@

maxq: MaxQ.cpp
	g++ -O3 $^ -I $(SAMTOOLS)  -L $(SAMTOOLS)  -lbam -lm -lz -o $@

bedq: MaxQBed.cpp
	g++ -O3 $^ -I $(SAMTOOLS)  -L $(SAMTOOLS)  -lbam -lm -lz -o $@

bamc: PrintCoverage.cpp
	g++ -O3 $^ -I $(SAMTOOLS)  -L $(SAMTOOLS)  -lbam -lm -lz -o $@

hardstop: PrintHardStop.cpp
	g++ -g -O3  $^ -I $(SAMTOOLS)  -L $(SAMTOOLS)  -lbam -lm -lz -o $@

splitread: PrintSplitReads.cpp
	g++ -g -O3  $^ -I $(SAMTOOLS)  -L $(SAMTOOLS)  -lbam -lm -lz -o $@

bamacc: ReadStats.cpp
	g++ -O3 $^ -I $(SAMTOOLS)  -L $(SAMTOOLS)  -lbam -lm -lz -o $@


