#include "argparser.h"
#include <cstdlib> // For atoi function
#include <iostream>
#include <filesystem>

void printHelp(const char *programName)
{
    std::cout << "Usage: " << programName << " [-k <kmer_size>] [-t <threshold>] [-m <max_read_len>]" << std::endl;
    std::cout << "           [-s <score_type>] [-i <file1> <file2> ...] [-o <output_dir>]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -h:                      To print this help" << std::endl;
    std::cout << "  -k <kmer_size>:          Set kmer size (default: 3)" << std::endl;
    std::cout << "  -t <threshold>:          Set threshold (default: 15)" << std::endl;
    std::cout << "  -m <max_read_len>:       Set the maximal read length (default: 160)" << std::endl;
    std::cout << "  -s <score_type>:         Set the score type [[max], acc, <kmer>]" << std::endl;
    std::cout << "                               max:    the largest repeat defines the score" << std::endl;
    std::cout << "                               acc:    repeats add to the score, interruptions" << std::endl;
    std::cout << "                                       decrease the score" << std::endl;
    std::cout << "                               <kmer>: repeats of this kmer add to the score," << std::endl;
    std::cout << "                                       interruptions and other kmers decrease the score" << std::endl;
    std::cout << "  -i <file1> <file2> ...:  Input files" << std::endl;
    std::cout << "  -o <output_dir>:         Output directory (default: ../out)" << std::endl;
}

Args parseArgs(int argc, char *argv[])
{
    Args args;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-k" && i + 1 < argc)
        {
            args.k = std::atoi(argv[++i]);
        }
        else if (arg == "-t" && i + 1 < argc)
        {
            args.threshold = std::atoi(argv[++i]);
        }
        else if (arg == "-m" && i + 1 < argc)
        {
            args.max_read_len = std::atoi(argv[++i]);
        }
        else if (arg == "-s" && i + 1 < argc)
        {
            args.score = argv[++i];
        }
        else if (arg == "-o" && i + 1 < argc)
        {
            args.output_dir = argv[++i];
        }
        else if (arg == "-h")
        {
            printHelp(argv[0]);
            std::exit(0);
        }
        else if (arg == "-i")
        {
            while (i + 1 < argc && argv[i + 1][0] != '-')
            {
                args.files.push_back(argv[++i]);
            }
        }
        else
        {
            std::cerr << "Error: Unknown option or missing argument: " << arg << std::endl;
            printHelp(argv[0]);
            std::exit(1);
        }
    }
    if (args.files.empty())
    {
        std::cerr << "Error: No input files" << std::endl;
        printHelp(argv[0]);
        std::exit(1);
    }

    if (args.score != "max" && args.score != "acc")
    {
        for (auto c : args.score)
        {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            {
                std::cerr << "Error: Unknown score type: " << args.score << std::endl;
                printHelp(argv[0]);
                std::exit(1);
            }
        }
    }

    return args;
}
