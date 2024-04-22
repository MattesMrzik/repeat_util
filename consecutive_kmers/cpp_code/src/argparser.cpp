#include "argparser.h"
#include <cstdlib> // For atoi function
#include <iostream>
#include <filesystem>

void printHelp(const char *programName)
{
  std::cout << R"(Usage: )" << programName << R"( [-k <kmer_size>] [-t <threshold>] [-m <max_seq_len>]
           [-s <score_type>] [-i <file1> <file2> ...] [-o <output_dir>]
Options:
  -h:                      To print this help
  -k <kmer_size>:          Set kmer size (default: 3)
  -t <threshold>:          Set threshold (default: 15)
  -m:                      If your sequences are all about the same length,
                           you may set this flag at runtime as well as
                           -DMAX_SEQ_LEN=[int] at compile time. This
                           program will then only consider the first
                           MAX_SEQ_LEN many bases of the sequence. This
                           may increase the performance of the program.
                           The current program was compiled with
                           MAX_SEQ_LEN set to )"
            << MAX_SEQ_LEN << R"(.
  -s <score_type>:         Set the score type [[max], acc, <kmer>]
                               max:    the largest repeat defines the score
                               acc:    repeats add to the score, interruptions
                                       decrease the score
                               <kmer>: repeats of this kmer add to the score,
                                       interruptions and other kmers decrease
                                       the score
  -v                       print verbose to stdout
  -i <file1> <file2> ...:  Input files
  -c:                      Output coordinates if input is a BAM file
  -r:                      Use reverse complement kmers for the repeats
                           when producing a BAM file
  -o <output_dir>:         Output directory (default: ../out)
)";
}

Args parseArgs(int argc, char *argv[])
{
  return parseArgs(argc, argv, false);
}

Args parseArgs(int argc, char *argv[], bool quiet)
{
  Args args;
  args.root_path = std::filesystem::canonical(std::filesystem::absolute(argv[0])).parent_path().parent_path();
  for (int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if (arg == "-k" && i + 1 < argc)
    {
      args.k = std::stoul(argv[++i]);
      if (args.k < 2)
      {
        std::cout << "Error: k must be greater than 1" << std::endl;
        std::exit(1);
      }
    }
    else if (arg == "-t" && i + 1 < argc)
    {
      args.threshold = std::atoi(argv[++i]);
    }
    else if (arg == "-s" && i + 1 < argc)
    {
      args.score = argv[++i];
    }
    else if (arg == "-m")
    {
      args.use_max_seq_len = true;
    }
    else if (arg == "-v")
    {
      args.verbose = true;
    }
    else if (arg == "-c")
    {
      args.bam_output_as_coords = true;
    }
    else if (arg == "-r")
    {
      args.use_reverse_complement_kmer_for_bam = true;
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
      if (!quiet)
      {
        std::cerr << "Error: Unknown option or missing argument: " << arg << std::endl;
        printHelp(argv[0]);
        std::exit(1);
      }
    }
  }
  if (args.files.empty())
  {
    if (!quiet)
    {
      std::cerr << "Error: No input files" << std::endl;
      printHelp(argv[0]);
      std::exit(1);
    }
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
