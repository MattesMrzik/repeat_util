# g++ -std=gnu++17 main.cpp argparser.cpp ConsecutiveKmers.cpp -o main -lhts -DSEQAN_HAS_ZLIB -lz
CXX = g++
CXXFLAGS = -std=gnu++17
LDFLAGS = -lhts -DSEQAN_HAS_ZLIB -lz
MAX_SEQ_LEN ?= 160
CXXFLAGS += -DMAX_SEQ_LEN=$(MAX_SEQ_LEN)

SRCS = main.cpp argparser.cpp ConsecutiveKmers.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = main

TEST_SRCS = test.cpp argparser.cpp ConsecutiveKmers.cpp
TEST_OBJS = $(TEST_SRCS:.cpp=.o)
TEST_EXEC = test

.PHONY: all main test clean

all: $(EXEC)

main: $(OBJS)
		$(CXX) $(CXXFLAGS) -DMAX_SEQ_LEN=$(MAX_SEQ_LEN) $^ -o $@ $(LDFLAGS)

test: $(TEST_OBJS)
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)
		./test

%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
		rm -f $(OBJS) $(EXEC) $(TEST_OBJS) $(TEST_EXEC)
