SAMTOOLS=$(HOME)/software/samtools/samtools-0.1.18/
all: maxq bedq bamacc hardstop splitread coverage

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


