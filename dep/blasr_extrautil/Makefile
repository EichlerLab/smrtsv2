# Install PrintTupleCountTable from "extrautils" in the blasr repository

ENV_PACBIO = $(abspath ../conda/build/envs/pacbio)
ENV_TOOLS = $(abspath ../conda/build/envs/tools)

# Compile
bin/PrintTupleCountTable: blasr/extrautils/PrintTupleCountTable.cpp
	mkdir -p bin
	$(ENV_TOOLS)/bin/g++ \
		-std=c++11 -I$(ENV_TOOLS)/include \
		-I$(ENV_PACBIO)/include \
		-I$(ENV_PACBIO)/include/libblasr \
		-L$(ENV_TOOLS)/lib \
		-Wl,-rpath=$(ENV_TOOLS)/lib \
		-L$(ENV_PACBIO)/lib \
		-Wl,-rpath=$(ENV_PACBIO)/lib \
		-lblasr \
		-o $@ \
		$^


# Clone
blasr/extrautils/PrintTupleCountTable.cpp:
	git clone -b 5.3.2 https://github.com/PacificBiosciences/blasr.git

# Clean
.PHONY: clean
clean:
	rm -rf bin blasr
