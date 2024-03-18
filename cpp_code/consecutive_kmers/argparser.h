#ifndef ARGSPARSER_H
#define ARGSPARSER_H

#include <string>
#include <vector>

struct Args {
    int k = 3;               // Default value for k
    int threshold = 8;               // Default value for t
    std::vector<std::string> files; // Vector to store input files
    size_t max_read_len = 160; // Default value for m
};

void printHelp(const char* programName);
Args parseArgs(int argc, char* argv[]);

#endif // ARGSPARSER_H