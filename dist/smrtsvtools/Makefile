PREFIX = $(abspath ..)

CXXFLAGS += -I$(abspath include) -I$(abspath ../include) --std=c++14 -g

#LDFLAGS = -L$(abspath ../lib)
LDFLAGS += -lboost_program_options -lrt -lpthread -lbz2 -lhts -llzma

ifdef LD_LIBRARY_PATH
    LDFLAGS += $(foreach lpath,$(subst :, ,$(LD_LIBRARY_PATH)),-L$(lpath) -Wl,-rpath=$(lpath))
endif

BIN_OUT = $(addprefix bin/, alignfixup phasereads)

.PHONY: all
all: ${BIN_OUT}
	install bin/* ${PREFIX}/bin

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
