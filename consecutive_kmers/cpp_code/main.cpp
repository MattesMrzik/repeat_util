#ifndef MAX_READ_LEN
#define MAX_READ_LEN 160 // Default value
#endif
#define RESULT_LEN (MAX_READ_LEN * 4 / 3)

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <omp.h>
#include <map>
#include <utility> // for std::make_pair
#include <chrono>
#include <thread>
#include <cstring> // for std::strcpy
#include <memory>
#include <filesystem>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include <seqan/seq_io.h>

namespace fs = std::filesystem;

#include "argparser.h"

Args args;

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
                   size_t &max_kmer_score)
{

    append_string_type(result, "(");
    append_string_type(result, current_kmer);
    append_string_type(result, ")_");
    append_string_type(result, std::to_string(n_found + 1).c_str());
    append_string_type(result, " ");
    // dont forget to decrease if no repeat was found in get_repeats()
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

/// @brief calculates a string that is the input sequence, but repeating k-mers are replaced by (kmer)_n.
/// @param seq the input sequence.
/// @param frame the frame to start from (0, 1, 2, ... , k-1).
/// @param result the result string. Must be long enough to hold the result.
/// @return the score of the sequence, i.e. the maximum number of consecutive k-mers.
template <typename ResultType>
size_t get_repeats(const std::string &seq,
                   int frame,
                   ResultType &result)
{
    size_t l = seq.length();
    size_t k = args.k;
    if (seq.length() > MAX_READ_LEN)
    {
        l = MAX_READ_LEN;
    }
    size_t n_found = 0;
    size_t max_n_found = 0;
    size_t acc_score = 0;
    size_t max_acc_score = 0;
    size_t kmer_score = 0;
    size_t max_kmer_score = 0;
    size_t current_result_len = 0;
    char current_kmer[k + 1];
    current_kmer[k] = '\0';
    init_string_type(result, seq.substr(0, frame).c_str());
    size_t last_position_to_check = l - ((l - frame) % k) - k;
    for (int i = frame; i < last_position_to_check; i += k)
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
    size_t rest_start = last_position_to_check;
    size_t rest_length = l - last_position_to_check;
    if (n_found > 0)
    {
        append_repeat(result, current_kmer, n_found, acc_score, max_acc_score, kmer_score, max_kmer_score);
        rest_length -= k;
    }

    append_string_type(result, seq.substr(rest_start, rest_length).c_str());

    if (args.score == "acc")
    {
        return max_acc_score;
    }
    else if (args.score == "max")
    {
        return n_found + 1;
    }
    return max_kmer_score;
}

/// @brief iterates over the frames of a sequence and writes the result to a file.
/// @param seq_name the name of the sequence.
/// @param seq the sequence.
/// @param outfile the file to write the results to.
void iterate_over_frames(const std::string &seq_name,
                         const std::string &seq,
                         std::ofstream &outfile)
{
    if (seq.length() < args.k * 2)
    {
        // std::cerr << "The sequence " << seq_name << " is too short to contain a kmer of length " << k << std::endl;
        return;
    }

    char result[RESULT_LEN]; // this gets overwritten in when processing the next sequence

    for (int frame = 0; frame < args.k; ++frame)
    {
        if (MAX_READ_LEN == 0)
        {
            std::string result_as_string;
            int score = get_repeats(seq, frame, result_as_string);
            if (score >= args.threshold)
            {
                outfile << seq_name
                        << ", frame: " << frame
                        << ", " << result_as_string
                        << ", score_type: " << args.score
                        << ", score: " << score
                        << ", seqlen too long: " << (seq.length() > MAX_READ_LEN) << std::endl;
            }
        }
        else
        {
            int score = get_repeats(seq, frame, result);
            if (score >= args.threshold)
            {
                outfile << seq_name
                        << ", frame: " << frame
                        << ", " << result
                        << ", score_type: " << args.score
                        << ", score: " << score
                        << ", seqlen too long: " << (seq.length() > args.max_read_len) << std::endl;
            }
        }
    }
}

/// @brief for /path/to/inputfile get /path/to/output_dir/inputfile.out
/// @param inputPath the path to the input file.
/// @return the path to the output file.
std::string getOutputFileName(const std::string &inputPath)
{
    fs::path inputFilePath(inputPath);
    std::string filename = inputFilePath.filename().string();
    fs::path outputPath(args.output_dir);
    outputPath /= filename + ".out";
    return outputPath.string();
}

/// @brief scans a file for sequences and calls iterate_over_frames() to find repeats for each sequence.
/// @param file_name the name of the file to scan.
/// @param outfile the file to write the results to.
void scan_fasta_and_fastq(std::string file_name)
{
    std::cout << "scanning file: " << file_name << std::endl;
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

    std::string line; // maybe not necessary if seqs are always on one line in the input files
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
                iterate_over_frames(seq_name, seq, outfile);
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
        iterate_over_frames(seq_name, seq, outfile);
    }
    outfile.close();
}

void scan_bam(std::string filename)
{
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

    // Initialize the htslib BAM file handler
    bam_hdr_t *header = sam_hdr_read(bamFile);
    bam1_t *record = bam_init1();

    size_t max_size = 0;
    // Read alignments from the compressed BAM file
    size_t count = 0;
    while (sam_read1(bamFile, header, record) >= 0)
    {
        count++;
        if (count % 1000000 == 0)
        {
            std::cout << "Scanned " << count / 1000000 << " million sequences in file " << filename << std::endl;
        }
        std::string seq_name = bam_get_qname(record);
        uint8_t *succinct_seq = bam_get_seq(record);
        std::string seq;
        for (int i = 0; i < record->core.l_qseq; ++i)
        {
            seq += seq_nt16_str[bam_seqi(succinct_seq, i)];
        }

        iterate_over_frames(seq_name, seq, outfile);
        // Process each alignment record here
        // std::cout << "Read name: " << bam_get_qname(record) << std::endl;
        // std::cout << "Read cigar: " << bam_get_cigar(record) << std::endl;
        // std::cout << "postion: " << record->core.pos << std::endl;
        // std::cout << "end position: " << bam_endpos(record) << std::endl;

        // You can access other fields of the alignment record as needed
    }

    // Clean up
    bam_destroy1(record);
    bam_hdr_destroy(header);
    hts_close(bamFile);
}

// void read_fasta(std::string filename)
// {
//     faidx_t *fastaIndex = fai_load(filename.c_str());
//     if (!fastaIndex)
//     {
//         std::cerr << "Error: Failed to open the compressed FASTA file." << std::endl;
//         return;
//     }

//     // Read sequences from the compressed FASTA file
//     int seq_len;
//     char *sequence = fai_fetch(fastaIndex, "sequence_name", &seq_len);
//     if (sequence)
//     {
//         std::cout << "Sequence: " << sequence << std::endl;
//         free(sequence); // Release memory allocated by fai_fetch
//     }

//     // Clean up
//     fai_destroy(fastaIndex);
// }

void read_fasta_fastq_gz(std::string filename)
{
    std::cout << "Reading file: " << filename << std::endl;
    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, filename.c_str()))
    {
        std::cerr << "Error: Failed to open the compressed file." << std::endl;
    }

    // Iterate over the records in the compressed file
    seqan::CharString id, seq, qual;
    while (!atEnd(seqFileIn))
    {
        readRecord(id, seq, qual, seqFileIn);

        // Process the record (e.g., print ID and sequence)
        std::cout << "ID: " << id << std::endl;
        std::cout << "Sequence: " << seq << std::endl;
    }

    // Close the compressed file
    close(seqFileIn);
}

void scan_file(std::string file_name)
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
        read_fasta_fastq_gz(file_name);
    }
    else
    {
        std::cerr << "Error: Unsupported file format " << extension << std::endl;
    }
}

int main(int argc, char *argv[])
{
    args = parseArgs(argc, argv);

    // #pragma omp parallel for
    for (auto &file_name : args.files)
    {
        scan_file(file_name);
    }
    return 0;
}