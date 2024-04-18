#include <iostream>
#include <string>
#include <filesystem>
#include <algorithm> // for std::sort

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

template <typename OutputStream>
void ConsecutiveKmers::get_repeat_coordinates(const std::string &seq_name,
                                              const std::string &seq,
                                              OutputStream &outfile,
                                              bool reverse_complement)

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
  char repeat[args.k];
  bool found_repeat = false;
  size_t found_repeat_start;
  for (size_t i = 0; i < seq.length(); i++)
  {
    if (found_repeat)
    {
      if (seq[i + args.k - 1] != seq[i + args.k + args.k - 1])
      {

        // repeat ends
        // TODO do I want to create a hash table to check if this repeat was already found by another read?
        // or simply write it to a file,
        outfile << seq_name << "\t" << found_repeat_start << "\t" << i + args.k + args.k - 2
                << "\t" << get_atomic_pattern(std::string(repeat, args.k), reverse_complement) << std::endl;
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
  if (found_repeat) // this is perhaps never called since std::string is always null terminated
  {
    std::cout << "Found repeat at the end of the sequence" << std::endl;
    outfile << seq_name << found_repeat_start << "\t" << seq.length() + args.k + args.k - 2
            << "\t" << get_atomic_pattern(std::string(repeat, args.k), reverse_complement) << std::endl;
  }
}

// since I do not call this method in this class with type ostringstream
// I need to explicitly instantiate the template
template void ConsecutiveKmers::get_repeat_coordinates(const std::string &seq_name,
                                                       const std::string &seq,
                                                       std::ostringstream &outfile,
                                                       bool reverse_complement);

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

std::string ConsecutiveKmers::getOutputFileName(const std::string &inputPath)
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
  std::string output_file_name = getOutputFileName(file_name);
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
  std::string output_file_name = getOutputFileName(filename);
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

    std::string bam_aux = std::string(chrom) + "; " + std::to_string(record->core.pos);
    iterate_over_frames(seq_name, seq, outfile, bam_aux); // overload method and add parameters for seq, pos
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

  std::string output_file_name = getOutputFileName(filename);
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
