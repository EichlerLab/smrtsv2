PREFIX = $(abspath ..)

# Note: Do not install the version of samblaster in bwakit. It will cause the
# genotyper to fail on some samples with:
# samblaster: New buffer length exceeds maximum while changing field value.


.PHONY: install
install: ${PREFIX}/bin/bwa-postalt.js ${PREFIX}/bin/k8 ${PREFIX}/bin/seqtk

${PREFIX}/bin/bwa-postalt.js ${PREFIX}/bin/k8 ${PREFIX}/bin/seqtk: bwa.kit/bwa
	mkdir -p ${PREFIX}/bin
	install bwa.kit/k8 ${PREFIX}/bin
	install bwa.kit/seqtk ${PREFIX}/bin
	install bwa.kit/bwa-postalt.js ${PREFIX}/bin
	install bwa.kit/samblaster ${PREFIX}/bin

bwa.kit/bwa: bwakit-0.7.15_x64-linux.tar.bz2
	tar -xjf bwakit-0.7.15_x64-linux.tar.bz2

bwakit-0.7.15_x64-linux.tar.bz2:
	wget https://downloads.sourceforge.net/project/bio-bwa/bwakit/bwakit-0.7.15_x64-linux.tar.bz2
