#include "argparser.h"
#include <cstdlib> // For atoi function
#include <iostream>
#include <filesystem>

void printHelp(const char* programName) {
    std::cout << "Usage: " << programName << " [-k <kmer_size>] [-t <threshold>] [-i <file1> <file2> ...]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -k <kmer_size>: Set kmer size (default: 3)" << std::endl;
    std::cout << "  -t <threshold>: Set threshold (default: 8)" << std::endl;
    std::cout << "  -m <max_read_len>: Set the maximal read length (default: 160)" << std::endl;
    std::cout << "  -i <file1> <file2> ...: Input files" << std::endl;
}

Args parseArgs(int argc, char* argv[]) {
    Args args;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-k" && i + 1 < argc) {
            args.k = std::atoi(argv[++i]);
        } else if (arg == "-t" && i + 1 < argc) {
            args.threshold = std::atoi(argv[++i]);
        } else if (arg == "-m" && i + 1 < argc) {
            args.max_read_len = std::atoi(argv[++i]);
        } else if (arg == "-i") {
            while (i + 1 < argc && argv[i + 1][0] != '-') {
                args.files.push_back(argv[++i]);
            }
        } else {
            std::cerr << "Error: Unknown option or missing argument: " << arg << std::endl;
            printHelp(argv[0]);
            std::exit(1);
        }
    }
    if (args.files.empty()) {
        std::cerr << "Error: No input files" << std::endl;
        printHelp(argv[0]);
        std::exit(1);
    }

    return args;
}
