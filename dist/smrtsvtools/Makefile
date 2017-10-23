PREFIX = $(abspath ..)

CXXFLAGS = -Iinclude -I../include --std=c++14

LDFLAGS = -L../lib
LDLIBS = -lboost_program_options -lrt -lpthread -lbz2 -lhts

.PHONY: all
all: bin/alignfixup bin/phasereads

bin/phasereads: build/phasereads.o build/phase/PhaseTable.o build/phase/SNVEntry.o build/phase/ReadPhaser.o build/phase/PhaseWriter.o build/phase/StringUtil.o
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
