#ifndef ARGSPARSER_H
#define ARGSPARSER_H

#ifndef MAX_SEQ_LEN
#define MAX_SEQ_LEN 160 // Default value
#endif
#define RESULT_LEN (MAX_SEQ_LEN * 4 / 3)

#include <string>
#include <vector>

struct Args
{
    int k = 3;                          // Default value for k
    int threshold = 15;                 // Default value for t
    std::vector<std::string> files;     // Vector to store input files
    std::string score = "max";          // Default value for the score type [max, acc, <kmer>]
    std::string output_dir = "../out/"; // Default value for the output directory
    size_t max_seq_len = MAX_SEQ_LEN;  // Default value for the maximal read length
    size_t result_len = RESULT_LEN;     // Default value for the result length
    bool use_max_seq_len = false;      // Default value for the use of maximal read length
    bool verbose = false;               // Default value for the verbose mode
};

void printHelp(const char *programName);
Args parseArgs(int argc, char *argv[]);

#endif // ARGSPARSER_H
