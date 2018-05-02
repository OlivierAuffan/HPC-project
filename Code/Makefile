CXX=mpic++ -Wall -O3

CPPFLAGS=-Iinc

LDFLAGS=-lboost_program_options -lm

FILES=export forward init memory parse_args shalw
OBJECTS=$(addsuffix .o, $(FILES))
BIN=bin/shalw

all : $(BIN)

$(BIN) : $(addprefix obj/, $(OBJECTS))
	$(CXX) -g -o $@ $^ $(LDFLAGS)

obj/%.o : src/%.cpp
	$(CXX) -g -c -o $@ $^ $(CPPFLAGS)

obj/%.o : src/%.c
	$(CXX) -g -c -o $@ $^ $(CPPFLAGS)

clean :
	rm -f bin/* obj/*
