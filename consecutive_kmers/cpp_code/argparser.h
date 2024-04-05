#ifndef ARGSPARSER_H
#define ARGSPARSER_H

#include <string>
#include <vector>

struct Args {
    int k = 3;               // Default value for k
    int threshold = 15;               // Default value for t
    std::vector<std::string> files; // Vector to store input files
    size_t max_read_len = 0; // Default value for m
    std::string score = "max"; // Default value for the score type [max, acc, <kmer>]
    std::string output_dir = "../out/"; // Default value for the output directory
};

void printHelp(const char* programName);
Args parseArgs(int argc, char* argv[]);

#endif // ARGSPARSER_H
