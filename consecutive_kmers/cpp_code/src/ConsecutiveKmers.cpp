#include <iostream>
#include <string>
#include <filesystem>
#include <algorithm> // for std::sort
#include <stdexcept>

#include <htslib/hts.h>
#include <htslib/sam.h>
// #include <htslib/faidx.h>
#include <seqan/seq_io.h>

#include "ConsecutiveKmers.h"

namespace fs = std::filesystem;

#include "argparser.h"

inline void init_string_type(char *str, const char *init)
{
  std::strcpy(str, init);
}

inline void init_string_type(std::string &str, const char *init)
{
  str = init;
}

inline void append_string_type(char *str, const char *append)
{
  std::strcat(str, append);
}

inline void append_string_type(std::string &str, const char *append)
{
  str.append(append);
}

std::string ConsecutiveKmers::expand_collapsed_repeats(const std::string &str)
{
  std::string result;
  size_t pos = 0;

  while (pos < str.size())
  {
    if (str[pos] == '(')
    {
      size_t closingPos = str.find(')', pos);
      if (closingPos == std::string::npos)
      {
        return "Error: Closing parenthesis not found.";
      }
      std::string repeat = str.substr(pos + 1, closingPos - pos - 1);
      size_t space_pos = str.find(' ', pos);
      if (space_pos == std::string::npos)
      {
        return "Error: Space not found after opening parenthesis.";
      }
      std::string repeatCountStr = str.substr(closingPos + 2, space_pos - closingPos - 2);
      int repeatCount = std::stoi(repeatCountStr);
      std::string toRepeat = repeat;
      for (int i = 0; i < repeatCount; ++i)
      {
        result += toRepeat;
      }
      pos = space_pos + 1;
    }
    else
    {
      result += str[pos];
      ++pos;
    }
  }
  return result;
}

template <typename ResultType>
void ConsecutiveKmers::append_repeat(ResultType &result,
                                     char *current_kmer,
                                     size_t &n_found,
                                     size_t &acc_score,
                                     size_t &max_acc_score,
                                     size_t &kmer_score,
                                     size_t &max_kmer_score)
{

  append_string_type(result, "(");
  append_string_type(result, current_kmer);
  append_string_type(result, ")_");
  append_string_type(result, std::to_string(n_found + 1).c_str());
  append_string_type(result, " ");
  // dont forget to decrease in get_repeats() if no repeat was found
  acc_score += n_found + 1;
  if (acc_score > max_acc_score)
  {
    max_acc_score = acc_score;
  }
  if (current_kmer == args.score)
  {
    kmer_score += n_found + 1;
    if (kmer_score > max_kmer_score)
    {
      max_kmer_score = kmer_score;
    }
  }
  else
  {
    if (kmer_score > n_found + 1)
    {
      kmer_score -= n_found + 1;
    }
  }
  n_found = 0;
}

template <typename ResultType>
size_t ConsecutiveKmers::get_repeats(const std::string &seq,
                                     int frame,
                                     ResultType &result)
{
  size_t l = seq.length();
  size_t k = args.k;
  if (seq.length() > args.max_seq_len && args.use_max_seq_len)
  {
    l = args.max_seq_len;
  }
  size_t n_found = 0;
  size_t max_n_found = 0;
  size_t acc_score = 0;
  size_t max_acc_score = 0;
  size_t kmer_score = 0;
  size_t max_kmer_score = 0;
  char current_kmer[k + 1];
  current_kmer[k] = '\0';
  init_string_type(result, seq.substr(0, frame).c_str());
  size_t last_position_to_check = l - ((l - frame) % k) - k - k;

  for (size_t i = frame; i <= last_position_to_check; i += k)
  {
    bool is_repeat = true;
    for (size_t j = 0; j < k; j++) // this is me trying to make it lightweight
    {
      current_kmer[j] = seq[i + j];
      if (seq[i + j] != seq[i + k + j])
      {
        is_repeat = false;
      }
    }
    if (is_repeat)
    {
      n_found++;
      if (n_found > max_n_found)
      {
        max_n_found = n_found;
      }
    }
    else
    {
      if (n_found > 0)
      {
        append_repeat(result, current_kmer, n_found, acc_score, max_acc_score, kmer_score, max_kmer_score);
      }
      else
      {
        append_string_type(result, current_kmer);
        if (acc_score > 0)
        {
          acc_score -= 1;
        }
        if (kmer_score > 0)
        {
          kmer_score -= 1;
        }
      }
    }
  }
  size_t remainder_start;
  size_t remainder_length;
  if (n_found > 0)
  {
    append_repeat(result, current_kmer, n_found, acc_score, max_acc_score, kmer_score, max_kmer_score);
    remainder_start = last_position_to_check + k + k;
    remainder_length = l - remainder_start;
  }
  else
  {
    remainder_start = last_position_to_check + k;
    remainder_length = l - remainder_start;
  }
  append_string_type(result, seq.substr(remainder_start, remainder_length).c_str());

  if (args.score == "acc")
  {
    return max_acc_score;
  }
  else if (args.score == "max")
  {
    if (max_n_found > 0)
    {
      return max_n_found + 1;
    }
    return 0;
  }
  return max_kmer_score;
}

std::string reverseComplement(const std::string &kmer)
{
  std::string reverseComp;
  for (int i = kmer.length() - 1; i >= 0; --i)
  {
    char base = kmer[i];
    // Complement each nucleotide base
    switch (base)
    {
    case 'A':
      reverseComp += 'T';
      break;
    case 'T':
      reverseComp += 'A';
      break;
    case 'C':
      reverseComp += 'G';
      break;
    case 'G':
      reverseComp += 'C';
      break;
    default:
      // If the base is not A, T, C, or G, handle it as needed
      break;
    }
  }
  return reverseComp;
}

void ConsecutiveKmers::init_atomic_patterns()
{
  // TODO this may be faster than calculating the atomic patterns on the fly
  // since with this method i dont need to do the lookup
}

// lets see if this can be called with char* since this may produce an r value which perhaps cant be passed by reference
std::string ConsecutiveKmers::get_atomic_pattern(const std::string &kmer, bool reverse_complement)
{
  // check if the kmer is in the atomic pattern map
  if (atomic_patterns.find(kmer) != atomic_patterns.end())
  {
    return atomic_patterns[kmer];
  }
  else
  {
    if (args.verbose)
    {
      std::cout << "Calculating atomic pattern for " << kmer << std::endl;
    }
    int size;
    if (reverse_complement)
    {
      size = args.k * 2;
    }
    else
    {
      size = args.k;
    }
    std::string all[size];
    for (size_t i = 0; i < args.k; i++)
    {
      if (reverse_complement)
      {
        all[2 * i] = kmer.substr(i) + kmer.substr(0, i);
        all[2 * i + 1] = reverseComplement(kmer.substr(i) + kmer.substr(0, i));
      }
      else
      {
        all[i] = kmer.substr(i) + kmer.substr(0, i);
      }
    }
    std::sort(all, all + sizeof(all) / sizeof(all[0]));
    atomic_patterns[kmer] = all[0];
    return atomic_patterns[kmer];
  }
}

std::string ConsecutiveKmers::cigar_array_to_str(const uint32_t *cigar, int numCigarOps)
{
  std::stringstream ss;
  for (int i = 0; i < numCigarOps; ++i)
  {
    int opLen = bam_cigar_oplen(cigar[i]);
    char opChar = bam_cigar_opchr(cigar[i]);
    ss << opLen << opChar;
  }
  return ss.str();
}

uint32_t ConsecutiveKmers::n_cigar_from_str(const std::string &cigar_string)
{
  uint32_t n_cigar = 0;
  for (auto c : cigar_string)
  {
    if (c == 'M' || c == 'I' || c == 'D' || c == 'N' || c == 'S' || c == 'H' || c == 'P' || c == '=' || c == 'X')
    {
      n_cigar++;
    }
  }
  return n_cigar;
}

std::string ConsecutiveKmers::cigar_vector_to_str(const std::vector<uint32_t> &cigarArray)
{
  std::string cigarString;
  for (uint32_t cigarElem : cigarArray)
  {
    int cigarLen = cigarElem >> 4; // Extract length from upper 28 bits
    char cigarOp;
    switch (cigarElem & 0x0F)
    { // Extract operation from lower 4 bits
    case 0:
      cigarOp = 'M';
      break;
    case 1:
      cigarOp = 'I';
      break;
    case 2:
      cigarOp = 'D';
      break;
    case 3:
      cigarOp = 'N';
      break;
    case 4:
      cigarOp = 'S';
      break;
    case 5:
      cigarOp = 'H';
      break;
    case 6:
      cigarOp = 'P';
      break;
    case 7:
      cigarOp = '=';
      break;
    case 8:
      cigarOp = 'X';
      break;
    default:
      cigarOp = '?';
      break;
    }
    cigarString += std::to_string(cigarLen) + cigarOp;
  }
  return cigarString;
}

std::vector<uint32_t> ConsecutiveKmers::cigar_str_to_vector(const std::string &cigarString)
{
  std::vector<uint32_t> cigarArray;
  size_t pos = 0;
  while (pos < cigarString.length())
  {
    size_t nextPos = cigarString.find_first_of("MIDNSHP=X", pos);
    if (nextPos == std::string::npos)
    {
      break;
    }
    int cigarLen = std::stoi(cigarString.substr(pos, nextPos - pos));
    char cigarOp = cigarString[nextPos];
    uint32_t cigarElem = cigarLen << 4; // Shift the length into the upper 28 bits
    switch (cigarOp)
    {
    case 'M':
      cigarElem |= 0;
      break;
    case 'I':
      cigarElem |= 1;
      break;
    case 'D':
      cigarElem |= 2;
      break;
    case 'N':
      cigarElem |= 3;
      break;
    case 'S':
      cigarElem |= 4;
      break;
    case 'H':
      cigarElem |= 5;
      break;
    case 'P':
      cigarElem |= 6;
      break;
    case '=':
      cigarElem |= 7;
      break;
    case 'X':
      cigarElem |= 8;
      break;
    }
    cigarArray.push_back(cigarElem);
    pos = nextPos + 1;
  }
  return cigarArray;
}

uint32_t *ConsecutiveKmers::cigar_str_to_array(const std::string &cigarString)
{
  std::vector<uint32_t> cigarArray = cigar_str_to_vector(cigarString);
  uint32_t *cigar = new uint32_t[cigarArray.size()];
  std::copy(cigarArray.begin(), cigarArray.end(), cigar);
  return cigar;
}
std::vector<uint32_t> ConsecutiveKmers::get_aligned_reference_positions(bam1_t *read,
                                                                        const std::vector<uint32_t> &readPositions)
{
  uint32_t ref_start = read->core.pos;
  uint32_t *cigar = bam_get_cigar(read);
  auto n_cigar = read->core.n_cigar;
  return get_aligned_reference_positions(ref_start, cigar, n_cigar, readPositions); // TOD maybe return a shared pointer
}

// https://samtools.github.io/hts-specs/SAMv1.pdf
std::vector<uint32_t> ConsecutiveKmers::get_aligned_reference_positions(uint32_t const ref_start,
                                                                        uint32_t *cigar_array,
                                                                        size_t const n_cigar,
                                                                        const std::vector<uint32_t> &read_positions)
{
  std::vector<uint32_t> aligned_ref_positions;
  if (read_positions.empty())
  {
    return aligned_ref_positions;
  }

  for (size_t i = 1; i < read_positions.size(); i++)
  {
    if (read_positions[i] < read_positions[i - 1])
    {
      throw std::runtime_error("[Error] read_positions are not sorted!");
    }
  }
  uint32_t read_pos_idx = 0;
  uint32_t read_pos = read_positions[0];
  uint32_t previous_read_pos = 0;
  size_t pos_in_ref = 0;
  uint32_t pos_in_query = 0;

  for (size_t i = 0; i < n_cigar; ++i)
  {
    uint32_t op = bam_cigar_op(cigar_array[i]);
    uint32_t len = bam_cigar_oplen(cigar_array[i]);
    bool ref_moved = false;
    // bool query_moved = false;
    switch (op)
    {
    case 0: // M
    case 7: // =
    case 8: // X
      pos_in_query += len;
      pos_in_ref += len;
      // query_moved = true;
      ref_moved = true;
      break;
    case 1: // I
    case 4: // S
      pos_in_query += len;
      // query_moved = true;
      break;
    case 2: // D
    case 3: // N
      pos_in_ref += len;
      ref_moved = true;
      break;
    case 5: // H: H can only be present as the first and/or last operation
    case 6: // P
      break;
    }

    // Check if the current read position is covered by this CIGAR operation
    while (read_pos < pos_in_query)
    {
      uint32_t pos;
      if (!ref_moved)
      {
        pos = ref_start + pos_in_ref;
      }
      else
      {
        pos = ref_start + (pos_in_ref + read_pos - pos_in_query);
      }
      aligned_ref_positions.push_back(pos);
      read_pos_idx++;
      if (read_pos_idx < read_positions.size())
      {
        previous_read_pos = read_pos;
        read_pos = read_positions[read_pos_idx];
      }
      else
      {
        return aligned_ref_positions;
      }
      if (previous_read_pos > read_pos)
      {
        // maybe this is never reached
        std::cerr << "This is not implemented: reach out to mattes@mrzik.de with code z54j09hrijtsbv08q";
        exit(1);
      }
    }
  }
  if (read_pos_idx < read_positions.size())
  {
    std::string msg = "[Error] query index is not covered by CIGAR string!\n";
    msg += "Query index: " + std::to_string(read_positions[read_pos_idx]) + "\n";
    if (n_cigar > 0)
    {
      msg += "CIGAR: " + cigar_array_to_str(cigar_array, n_cigar) + "\n";
    }
    else
    {
      msg += "CIGAR: empty\n";
    }
    throw std::runtime_error(msg);
  }
  return aligned_ref_positions;
}

template <typename OutputStream>
void ConsecutiveKmers::write_repeat_coordinate(const std::string &seq_name,
                                               const std::string &seq,
                                               OutputStream &outfile,
                                               const std::string &chrom,
                                               uint32_t found_repeat_start,
                                               uint32_t end_pos,
                                               const std::string &repeat,
                                               uint32_t seq_position,
                                               uint32_t *cigar_array,
                                               size_t n_cigar,
                                               bool reverse_complement)
{
  if (args.verbose)
  {
    std::cout << "[Verbose] "
              << "Found repeat: " << repeat << " at position " << found_repeat_start << " to " << end_pos << std::endl;
  }
  std::vector<uint32_t> aligned_ref_positions;
  try
  {
    aligned_ref_positions = get_aligned_reference_positions(seq_position,
                                                            cigar_array,
                                                            n_cigar,
                                                            std::vector<uint32_t>{found_repeat_start, end_pos});
  }
  catch (const std::exception &e)
  {
    std::string msg(e.what());
    msg += "for sequence " + seq_name + ": " + seq + ".\n";
    throw std::runtime_error( msg);
  }

  outfile << seq_name << '\t'
          << chrom << '\t'
          << found_repeat_start << '\t'
          << end_pos << '\t'
          << aligned_ref_positions[0] << '\t'
          << aligned_ref_positions[1] << '\t'
          << get_atomic_pattern(repeat, reverse_complement) << std::endl;
}

template void ConsecutiveKmers::write_repeat_coordinate(const std::string &seq_name,
                                                        const std::string &seq,
                                                        std::ostringstream &outfile,
                                                        const std::string &chrom,
                                                        uint32_t found_repeat_start,
                                                        uint32_t end_pos,
                                                        const std::string &repeat,
                                                        uint32_t seq_position,
                                                        uint32_t *cigar_array,
                                                        size_t n_cigar,
                                                        bool reverse_complement);

template <typename OutputStream>
void ConsecutiveKmers::write_repeat_coordinates(const std::string &seq_name,
                                                const std::string &seq,
                                                OutputStream &outfile,
                                                const std::string &chrom,
                                                uint32_t seq_position,
                                                uint32_t *cigar_array,
                                                const uint32_t n_cigar,
                                                bool reverse_complement)

{
  // TODO: dont pass seq as string but as byte array

  // std::cout << "Processing " << seq_name << std::endl;
  // std::cout << "seq: " << seq << std::endl;
  // std::cout << "chrom: " << chrom << std::endl;
  // std::cout << "position: " << position << std::endl;
  // std::cout << "cigar: " << cigar << std::endl;
  if (seq.length() < args.k * 2)
  {
    if (args.verbose)
    {
      std::cerr << "[Verbose] "
                << "Processing " << seq_name << " failed because the sequence is too short to contain a kmer of length " << args.k << std::endl;
    }
    return;
  }
  if (args.verbose)
  {
    std::cout << "[Verbose] "
              << "Processing " << seq_name
              << ", seq: " << seq << "." << std::endl;
  }
  char repeat[args.k];
  bool found_repeat = false;
  uint32_t found_repeat_start;
  for (size_t i = 0; i < seq.length() - args.k - args.k + 1; i++)
  {
    if (found_repeat)
    {
      if (seq[i + args.k - 1] != seq[i + args.k + args.k - 1])
      {
        // repeat ends
        // TODO do I want to create a hash table to check if this repeat was already found by another read?
        // or simply write it to a file,
        // TODO instead of calling get_aligned_reference_positions for every repeat separately
        // accumulate all repeat positions and call get_aligned_reference_positions only once

        uint32_t end_pos = i + args.k + args.k - 2;
        write_repeat_coordinate(seq_name,
                                seq,
                                outfile,
                                chrom,
                                found_repeat_start,
                                end_pos,
                                std::string(repeat, args.k),
                                seq_position,
                                cigar_array,
                                n_cigar,
                                reverse_complement);
        found_repeat = false;
        i = i + args.k + args.k - 2;
      }
    }
    else
    {
      bool continue_outer_loop = false;
      for (size_t j = i + args.k; j > i; j--)
      {
        size_t j_prime = j - 1;
        if (seq[j_prime] != seq[j_prime + args.k])
        {
          i = j_prime;
          found_repeat = false;
          continue_outer_loop = true;
          continue;
        }
        else
        {
          repeat[j_prime - i] = seq[j_prime];
        }
      }
      if (continue_outer_loop)
      {
        continue;
      }
      found_repeat_start = i;
      found_repeat = true;
    }
  }
  if (found_repeat)
  {
    write_repeat_coordinate(seq_name,
                            seq,
                            outfile,
                            chrom,
                            found_repeat_start,
                            static_cast<uint32_t>(seq.length() - 1),
                            std::string(repeat, args.k),
                            seq_position,
                            cigar_array,
                            n_cigar,
                            reverse_complement);
  }
}

// since I do not call this method in this class with type ostringstream
// I need to explicitly instantiate the template
template void ConsecutiveKmers::write_repeat_coordinates(const std::string &seq_name,
                                                         const std::string &seq,
                                                         std::ostringstream &outfile,
                                                         const std::string &chrom,
                                                         uint32_t seq_position,
                                                         uint32_t *cigar_array,
                                                         const uint32_t n_cigar,
                                                         bool reverse_complement = false);

void ConsecutiveKmers::iterate_over_frames(const std::string &seq_name,
                                           const std::string &seq,
                                           std::ofstream &outfile)
{
  iterate_over_frames(seq_name, seq, outfile, "");
}

void ConsecutiveKmers::iterate_over_frames(const std::string &seq_name,
                                           const std::string &seq,
                                           std::ofstream &outfile,
                                           const std::string &bam_aux)
{
  if (seq.length() < args.k * 2)
  {
    if (args.verbose)
    {
      std::cerr << "[Verbose] "
                << "Processing " << seq_name << " failed because the sequence is too short to contain a kmer of length " << args.k << std::endl;
    }
    return;
  }

  // TODO put this array in the loop, and make 2 separate loops, one for the max_seq_len and one for the normal case
  char result[RESULT_LEN]; // this gets overwritten in when processing the next sequence

  for (size_t frame = 0; frame < args.k; ++frame)
  {
    if (!args.use_max_seq_len)
    {
      if (args.verbose)
      {
        std::cerr << "[Verbose] "
                  << "Processing " << seq_name << " frame " << frame << " without a maximal sequence length" << std::endl;
      }
      std::string result_as_string;
      int score = get_repeats(seq, frame, result_as_string);
      if (score >= args.threshold)
      {
        outfile << seq_name
                << ", frame: " << frame
                << ", " << result_as_string
                << ", score_type: " << args.score
                << ", score: " << score
                << ", seqlen too long: " << (seq.length() > args.max_seq_len)
                << ", bam aux: " << bam_aux << std::endl;
      }
    }
    else
    {
      if (args.verbose)
      {
        std::cerr << "[Verbose] "
                  << "Processing " << seq_name << " frame " << frame << " with a maximal sequence length " << args.max_seq_len << std::endl;
      }
      int score = get_repeats(seq, frame, result);
      if (score >= args.threshold)
      {
        outfile << seq_name
                << ", frame: " << frame
                << ", " << result
                << ", score_type: " << args.score
                << ", score: " << score
                << ", seqlen too long: " << (seq.length() > args.max_seq_len)
                << ", bam aux: " << bam_aux << std::endl;
      }
    }
  }
}

std::string ConsecutiveKmers::get_output_file_name(const std::string &inputPath)
{
  fs::path inputFilePath(inputPath);
  std::string filename = inputFilePath.filename().string();
  fs::path outputPath(args.output_dir);
  outputPath /= filename + ".out";
  return outputPath.string();
}

void ConsecutiveKmers::scan_fasta_and_fastq(std::string file_name)
{
  if (args.verbose)
  {
    std::cerr << "[Verbose] "
              << "Processing file " << file_name << " as fasta or fastq" << std::endl;
  }
  std::string output_file_name = get_output_file_name(file_name);
  fs::create_directories(args.output_dir);
  std::ofstream outfile(output_file_name);

  if (!outfile.is_open())
  {
    std::cerr << "Error: Unable to open output file " << output_file_name + ".out" << std::endl;
    return;
  }
  // make the kmers parallel or is this too much overhead for to short seqs?
  std::ifstream file(file_name);
  if (!file.is_open())
  {
    std::cerr << "Error: Unable to open file " << file_name << std::endl;
    return;
  }

  std::string line; // maybe not necessary if sequences are always on one line in the input files
  std::string seq;  // maybe also convert to char array?
  std::string seq_name;

  bool in_quality = false; // in case it is a fastq file

  size_t count = 0;
  while (std::getline(file, line))
  {
    if (line.empty())
      continue;
    if (line[0] == '>' || line[0] == '@')
    {
      if (!seq.empty())
      {
        count++;
        if (count % 1000000 == 0)
        {
          std::cout << "Scanned " << count / 1000000 << " million sequences in file " << file_name << std::endl;
        }
        iterate_over_frames(seq_name.substr(1), seq, outfile);
        seq.clear();
      }
      seq_name = line;
      in_quality = false;
    }
    else if (line[0] == '+')
    {
      in_quality = true;
    }
    else if (!in_quality)
    {
      seq += line;
    }
  }
  if (!seq.empty())
  {
    iterate_over_frames(seq_name.substr(1), seq, outfile);
  }
  outfile.close();
}

void ConsecutiveKmers::scan_bam(std::string filename)
{
  if (args.verbose)
  {
    std::cerr << "[Verbose] "
              << "Processing file " << filename << " as BAM" << std::endl;
  }
  htsFile *bamFile = hts_open(filename.c_str(), "r");
  if (!bamFile)
  {
    std::cerr << "Error: Failed to open the BAM file " << filename << std::endl;
    return;
  }
  std::string output_file_name = get_output_file_name(filename);
  fs::create_directories(args.output_dir);
  std::ofstream outfile(output_file_name);
  if (!outfile.is_open())
  {
    std::cerr << "Error: Unable to open output file " << output_file_name + ".out" << std::endl;
    return;
  }

  bam_hdr_t *header = sam_hdr_read(bamFile);
  bam1_t *record = bam_init1();

  size_t count = 0;
  while (sam_read1(bamFile, header, record) >= 0)
  {
    count++;
    if (count % 1000000 == 0)
    {
      std::cout << "Scanned " << count / 1000000 << " million sequences in file " << filename << std::endl;
    }
    if (args.bam_output_as_coords)
    {
      // read is not aligned
      if (record->core.flag & 0x4)
      {
        continue;
      }
      // if read is secondary or supplementary
      if (record->core.flag & 0x100 || record->core.flag & 0x800)
      {
        continue;
      }
    }
    std::string seq_name = bam_get_qname(record);
    if (args.verbose)
    {
      std::cerr << "[Verbose] "
                << "Processing " << seq_name << std::endl;
    }
    uint8_t *succinct_seq = bam_get_seq(record);
    std::string seq;
    for (int i = 0; i < record->core.l_qseq; ++i)
    {
      seq += seq_nt16_str[bam_seqi(succinct_seq, i)];
    }
    uint32_t tid = record->core.tid; // Target ID (chromosome ID)
    const char *chrom;
    if (tid >= static_cast<uint32_t>(header->n_targets))
    {
      chrom = "unknown";
    }
    else
    {
      chrom = header->target_name[tid];
    }

    // TODO also allow option to not calculate the exact positions of repeats,
    // but just the starting coordinates of the read

    // TODO Bit 0x100 marks the alignment not to be used in certain analyses when the tools in use are aware
    // of this bit. It is typically used to flag alternative mappings when multiple mappings are presented
    // in a SAM
    if (args.bam_output_as_coords)
    {
      write_repeat_coordinates(seq_name,
                               seq,
                               outfile,
                               chrom,
                               record->core.pos,
                               bam_get_cigar(record),
                               record->core.n_cigar,
                               args.use_reverse_complement_kmer_for_bam);
    }
    else
    {
      std::string bam_aux = std::string(chrom) + "; " + std::to_string(record->core.pos);
      iterate_over_frames(seq_name, seq, outfile, bam_aux); // overload method and add parameters for seq, pos
    }
    // Process each alignment record here
    // std::cout << "Read name: " << bam_get_qname(record) << std::endl;
    // std::cout << "Read cigar: " << bam_get_cigar(record) << std::endl;
    // std::cout << "postion: " << record->core.pos << std::endl;
    // std::cout << "end position: " << bam_endpos(record) << std::endl;
    if (args.verbose)
    {
      std::cerr << "[Verbose] "
                << "Processed " << seq_name << std::endl;
    }
  }

  // Clean up
  bam_destroy1(record);
  bam_hdr_destroy(header);
  hts_close(bamFile);
}

void ConsecutiveKmers::scan_fasta_fastq_gz(std::string filename)
{
  if (args.verbose)
  {
    std::cerr << "[Verbose] "
              << "Processing file " << filename << " as gzipped fasta or fastq" << std::endl;
  }
  seqan::SeqFileIn seqFileIn;
  if (!open(seqFileIn, filename.c_str()))
  {
    std::cerr << "Error: Failed to open the compressed file." << std::endl;
  }

  std::string output_file_name = get_output_file_name(filename);
  fs::create_directories(args.output_dir);
  std::ofstream outfile(output_file_name);
  if (!outfile.is_open())
  {
    std::cerr << "Error: Unable to open output file " << output_file_name + ".out" << std::endl;
    return;
  }

  seqan::CharString seq_name, seq, qual;
  size_t count = 0;
  while (!atEnd(seqFileIn))
  {
    count++;
    if (count % 1000000 == 0)
    {
      std::cout << "Scanned " << count / 1000000 << " million sequences in file " << filename << std::endl;
    }
    readRecord(seq_name, seq, qual, seqFileIn);

    // Process the record (e.g., print ID and sequence)
    if (args.verbose)
    {
      std::cout << "sequence name: " << seq_name << std::endl;
      std::cout << "seq: " << seq << std::endl;
    }
    iterate_over_frames(toCString(seq_name), toCString(seq), outfile);
  }

  close(seqFileIn);
}

void ConsecutiveKmers::scan_file(std::string file_name)
{
  std::string extension = fs::path(file_name).extension().string();
  if (extension == ".bam")
  {
    scan_bam(file_name);
  }
  else if (extension == ".fasta" || extension == ".fa" || extension == ".fastq" || extension == ".fq")
  {
    scan_fasta_and_fastq(file_name);
  }
  else if (extension == ".gz")
  {
    scan_fasta_fastq_gz(file_name);
  }
  else
  {
    std::cerr << "Error: Unsupported file format " << extension << std::endl;
  }
}
