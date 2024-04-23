#ifndef CONSECUTIVE_KMERS_H
#define CONSECUTIVE_KMERS_H

#include "argparser.h"

#include <string>
#include <unordered_map>

#include <htslib/sam.h>

class ConsecutiveKmers
{
private:
  Args args;
  std::unordered_map<std::string, std::string> atomic_patterns;

public:
  ConsecutiveKmers(Args args)
  {
    this->args = args;
  }

  size_t get_atomic_pattern_size()
  {
    return atomic_patterns.size();
  }

  void init_atomic_patterns();

  static std::vector<uint32_t> cigar_str_to_array(const std::string &cigarString);

  static std::string cigar_array_to_str(const std::vector<uint32_t> &cigarArray);

  static std::string cigarToString(const uint32_t *cigar, int numCigarOps);

  static std::vector<uint32_t> get_aligned_reference_positions(bam1_t *read,
                                                               const std::vector<uint32_t> &readPositions);

  static std::vector<uint32_t> get_aligned_reference_positions(uint32_t const ref_start,
                                                               uint32_t *cigar,
                                                               size_t const n_cigar,
                                                               const std::vector<uint32_t> &readPositions);

  std::string get_atomic_pattern(const std::string &kmer, bool reverse_complement = false);

  template <typename OutputStream>
  void write_repeat_coordinates(const std::string &seq_name,
                                const std::string &seq,
                                OutputStream &outfile,
                                const std::string &chrom,
                                int position,
                                const std::string &cigar,
                                bool reverse_complement = false);

  /// @brief expands a string that contains collapsed repeats. For example, (GAT)_3 becomes GATGATGAT.
  /// @param str the input string.
  /// @return the expanded string.
  std::string expand_collapsed_repeats(const std::string &str);

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

  /// @brief iterates over the frames of a sequence and writes the result to a file.
  /// @param seq_name the name of the sequence.
  /// @param seq the sequence.
  /// @param outfile the file to write the results to.
  /// @param bam_aux the auxiliary information to write if the input is a BAM file.
  void iterate_over_frames(const std::string &seq_name,
                           const std::string &seq,
                           std::ofstream &outfile,
                           const std::string &bam_aux);

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
  std::string get_output_file_name(const std::string &inputPath);

  /// @brief scans a file for sequences and calls iterate_over_frames() to find repeats for each sequence.
  /// @param file_name the name of the file to scan.
  /// @param outfile the file to write the results to.
  void scan_fasta_and_fastq(std::string file_name);

  void scan_bam(std::string filename);

  void scan_fasta_fastq_gz(std::string filename);

  void scan_file(std::string file_name);
};

#endif // CONSECUTIVE_KMERS_H