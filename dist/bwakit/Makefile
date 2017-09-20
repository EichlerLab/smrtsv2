PREFIX = $(realpath ..)

.PHONY: install
install: ${PREFIX}/bin/bwa-postalt.js ${PREFIX}/bin/k8 ${PREFIX}/bin/seqtk ${PREFIX}/bin/samblaster

${PREFIX}/bin/bwa-postalt.js ${PREFIX}/bin/k8 ${PREFIX}/bin/seqtk ${PREFIX}/bin/samblaster: bwa.kit/bwa
	mkdir -p ${PREFIX}/bin
	install bwa.kit/k8 ${PREFIX}/bin
	install bwa.kit/seqtk ${PREFIX}/bin
	install bwa.kit/bwa-postalt.js ${PREFIX}/bin
	install bwa.kit/samblaster ${PREFIX}/bin

bwa.kit/bwa: bwakit-0.7.15_x64-linux.tar.bz2
	tar -xjvf bwakit-0.7.15_x64-linux.tar.bz2

bwakit-0.7.15_x64-linux.tar.bz2:
	wget https://downloads.sourceforge.net/project/bio-bwa/bwakit/bwakit-0.7.15_x64-linux.tar.bz2
