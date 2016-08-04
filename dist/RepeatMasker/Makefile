PWD  = $(shell pwd)
PERL = $(shell which perl)

RepeatMasker/RepeatMaskerConfig.pm: RepeatMasker
	cd $< && sed 's#/usr/local/hmmer#$(PWD)/../../bin#;s#TRF_PRGM = ""#TRF_PRGM = "$(PWD)/../../bin/trf"#;s#DEFAULT_SEARCH_ENGINE = "crossmatch"#DEFAULT_SEARCH_ENGINE = "hmmer"#' RepeatMaskerConfig.tmpl > RepeatMaskerConfig.pm && sed -i 's%#!/u1/local/bin/perl%#!$(PERL)%' RepeatMasker DateRepeats DupMasker ProcessRepeats RepeatProteinMask

RepeatMasker: RepeatMasker-open-4-0-5.tar.gz
	tar zxvf $< && touch $@

RepeatMasker-open-4-0-5.tar.gz:
	wget http://www.repeatmasker.org/$@ -O $@
