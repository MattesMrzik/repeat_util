# g++ -std=gnu++17 main.cpp argparser.cpp ConsecutiveKmers.cpp -o main -lhts -DSEQAN_HAS_ZLIB -lz
CXX = g++
CXXFLAGS = -std=gnu++17 -Wall
LDFLAGS = -lhts -DSEQAN_HAS_ZLIB -lz
MAX_SEQ_LEN ?= 160
CXXFLAGS += -DMAX_SEQ_LEN=$(MAX_SEQ_LEN)

SRCDIR = src
INCDIR = include
TESTDIR = test
BUILDDIR = build
BINDIR = bin

# SRCS and OBJS is all but main since i want to include them in the test
# and this avoids having two objects files with main functions
SRCS := $(filter-out $(SRCDIR)/main.cpp, $(wildcard $(SRCDIR)/*.cpp))
OBJS := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(SRCS))
EXEC = $(BINDIR)/main

TEST_SRCS := $(wildcard $(TESTDIR)/*.cpp)
TEST_OBJS := $(TEST_SRCS:$(TESTDIR)/%.cpp=$(BUILDDIR)/%.o)
TEST_EXEC = $(BINDIR)/test

.PHONY: all main test clean

all: $(EXEC) $(TEST_EXEC)

main: $(EXEC)

test: $(TEST_EXEC)
		$(TEST_EXEC)

$(EXEC): $(OBJS)  $(BUILDDIR)/main.o | $(BINDIR)
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TEST_EXEC): $(OBJS) $(TEST_OBJS) | $(BINDIR)
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
		$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ $(LDFLAGS)

$(BUILDDIR)/%.o: $(TESTDIR)/%.cpp | $(BUILDDIR)
		$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@ $(LDFLAGS)

$(BINDIR):
		mkdir -p $(BINDIR)

$(BUILDDIR):
		mkdir -p $(BUILDDIR)

clean:
		rm -rf $(OBJS) $(EXEC) $(TEST_OBJS) $(TEST_EXEC) $(BINDIR) $(BUILDDIR)
