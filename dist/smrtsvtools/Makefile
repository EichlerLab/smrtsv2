PREFIX = $(abspath ..)

CXXFLAGS = -I$(abspath include) -I$(abspath ../include) --std=c++14

LDFLAGS = -L$(abspath ../lib)
LDLIBS = -lboost_program_options -lrt -lpthread -lbz2 -lhts -llzma

BIN_OUT = $(addprefix bin/, alignfixup phasereads)

.PHONY: all
all: ${BIN_OUT}
	install bin/* ${PREFIX}/bin

bin/phasereads: build/phasereads.o build/phase/PhaseTable.o build/phase/SNVEntry.o build/phase/ReadPhaser.o build/phase/PhaseWriter.o build/phase/StringUtil.o build/util/err.o
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) ${LDFLAGS} ${LDLIBS} -o $@ $^

bin/alignfixup: build/alignfixup.o build/util/err.o
	mkdir -p $(dir $@)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${LDFLAGS} ${LDLIBS} -o $@ $^

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
