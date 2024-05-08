#include "argparser.h"
#include <cstdlib> // For atoi function
#include <iostream>
#include <fstream>
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
  -c:                      Output coordinates if input is a BAM file. If a repeat
                           is found, its atomic pattern is written to the output.
                           See also option -r to use reverse complement kmers in
                           the determination of the the atomic pattern.
  -r:                      Also use the reverse complement of the kmer in the
                           determination of the atomic pattern.
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
  if (!quiet)
  {
    to_file(args, args.output_dir);
  }
  return args;
}

std::string to_string(Args args)
{
  std::ostringstream string_stream;
  string_stream << "Time and date: " << get_time_and_date() << "\n";
  string_stream << "k: " << args.k << "\n";
  string_stream << "threshold: " << args.threshold << "\n";
  string_stream << "score: " << args.score << "\n";
  string_stream << "max_seq_len: " << MAX_SEQ_LEN << "\n";
  string_stream << "use_max_seq_len: " << args.use_max_seq_len << "\n";
  string_stream << "verbose: " << args.verbose << "\n";
  string_stream << "bam_output_as_coords: " << args.bam_output_as_coords << "\n";
  string_stream << "use_reverse_complement_kmer_for_bam: " << args.use_reverse_complement_kmer_for_bam << "\n";
  string_stream << "output_dir: " << args.output_dir << "\n";
  string_stream << "root_path: " << args.root_path << "\n";
  string_stream << "files: ";
  for (auto file : args.files)
  {
    string_stream << file << ", ";
  }
  string_stream << "\n";

  return string_stream.str();
}

std::string get_time_and_date()
{
  std::time_t now = std::time(nullptr);
  std::tm *time_info = std::localtime(&now);
  std::stringstream ss;
  ss << std::put_time(time_info, "%Y-%m-%d_%H.%M.%S");
  return ss.str();
}

void to_file(Args args, std::string output_dir)
{
  std::ofstream file;
  std::filesystem::path output_file_path(output_dir);
  output_file_path /= "log";
  std::filesystem::create_directories(output_file_path);
  output_file_path /= "args_log_" + get_time_and_date() + ".txt";
  file.open(output_file_path);
  file << to_string(args);
  file.close();
}
