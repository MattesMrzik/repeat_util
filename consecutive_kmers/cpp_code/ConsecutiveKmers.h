#ifndef CONSECUTIVE_KMERS_H
#define CONSECUTIVE_KMERS_H

#include <string>
#include "argparser.h"

class ConsecutiveKmers
{
private:
  Args args;
public:

  ConsecutiveKmers(Args args) {
    this->args = args;
  }

  /// @brief calculates a string that is the input sequence, but repeating k-mers are replaced by (kmer)_n.
  /// @param seq the input sequence.
  /// @param frame the frame to start from (0, 1, 2, ... , k-1).
  /// @param result the result string. Must be long enough to hold the result.
  /// @tparam ResultType the type of the result. Must be either char* or std::string.
  /// @return the score of the sequence, i.e. the maximum number of consecutive k-mers.
  template <typename ResultType>
  size_t get_repeats(const std::string &seq, int frame, ResultType &result);

  /// @brief iterates over the frames of a sequence and writes the result to a file.
  /// @param seq_name the name of the sequence.
  /// @param seq the sequence.
  /// @param outfile the file to write the results to.
  void iterate_over_frames(const std::string &seq_name,
                           const std::string &seq,
                           std::ofstream &outfile);

  /// @brief appends the repeat (kmer)_repeatSize to the result an updates the scores.
  /// @param result the result string.
  /// @param current_kmer the current kmer.
  /// @param n_found the number of consecutive repeats found.
  /// @param acc_score the accumulated score.
  /// @param max_acc_score the maximum accumulated score.
  /// @param kmer_score the kmer score.
  /// @param max_kmer_score the maximum kmer score.
  template <typename ResultType>
  void append_repeat(ResultType &result,
                     char *current_kmer,
                     size_t &n_found,
                     size_t &acc_score,
                     size_t &max_acc_score,
                     size_t &kmer_score,
                     size_t &max_kmer_score);

  /// @brief for /path/to/inputfile get /path/to/output_dir/inputfile.out
  /// @param inputPath the path to the input file.
  /// @return the path to the output file.
  std::string getOutputFileName(const std::string &inputPath);

  /// @brief scans a file for sequences and calls iterate_over_frames() to find repeats for each sequence.
  /// @param file_name the name of the file to scan.
  /// @param outfile the file to write the results to.
  void scan_fasta_and_fastq(std::string file_name);

  void scan_bam(std::string filename);

  void scan_fasta_fastq_gz(std::string filename);

  void scan_file(std::string file_name);
};

#endif // CONSECUTIVE_KMERS_H