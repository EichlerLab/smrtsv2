ENV_TOOLS = $(abspath ../conda/build/envs/tools)

CXX = $(ENV_TOOLS)/bin/g++

CXXFLAGS = -I$(abspath include) -I$(abspath ${PREFIX}/include) -std=c++11 -g

LDFLAGS = -L$(ENV_TOOLS)/lib -Wl,-rpath=$(ENV_TOOLS)/lib -lboost_program_options -lrt -lpthread -lbz2 -lhts -llzma

BIN_OUT = $(addprefix bin/, alignfixup phasereads tlenadd)

.PHONY: all
all: ${BIN_OUT}

bin/hdrtest: build/hdrtest.o build/util/err.o
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

bin/tlenadd: build/tlenadd.o build/util/err.o
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

bin/phasereads: build/phasereads.o build/phase/PhaseTable.o build/phase/SNVEntry.o build/phase/ReadPhaser.o build/phase/PhaseWriter.o build/phase/StringUtil.o build/util/err.o
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) ${LDFLAGS} -o $@ $^

bin/alignfixup: build/alignfixup.o build/util/err.o
	mkdir -p $(dir $@)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${LDFLAGS} -o $@ $^

build/%.o: src/%.cpp
	mkdir -p $(dir $@)
	${CXX} -c ${CPPFLAGS} ${CXXFLAGS} -o $@ $<

build/phase/%.o: src/phase/%.cpp include/phase/%.h
	mkdir -p $(dir $@)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

build/util/%.o: src/util/%.cpp include/util/%.h
	mkdir -p $(dir $@)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -rf build bin
