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

#include "argparser.h"

/// @brief calculates a string that is the input sequence, but repeating k-mers are replaced by (kmer)_n.
/// @param k the length of the k-mer.
/// @param seq the input sequence.
/// @param frame the frame to start from (0, 1, 2, ... , k-1).
/// @param result the result string. Must be long enough to hold the result.
/// @return the score of the sequence, i.e. the maximum number of consecutive k-mers.
size_t get_repeats(int k,
                   const std::string &seq,
                   int frame,
                   char *result,
                   size_t max_read_len,
                   bool only_used_prefix = false)
{
    size_t l = seq.length();
    if (only_used_prefix)
    {
        l = max_read_len;
    }
    size_t n_found = 0;
    size_t max_n_found = 0;
    size_t current_result_len = 0;
    char current_kmer[k + 1];
    current_kmer[k] = '\0';
    std::strcpy(result, seq.substr(0, frame).c_str());
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
                std::strcat(result, "(");
                std::strcat(result, current_kmer);
                std::strcat(result, ")_");
                std::strcat(result, std::to_string(n_found + 1).c_str());
                std::strcat(result, " ");
                n_found = 0;
            }
            else
            {
                std::strcat(result, current_kmer);
            }
        }
    }
    size_t rest_start = last_position_to_check;
    size_t rest_length = l - last_position_to_check;
    if (n_found > 0)
    {
        std::strcat(result, "(");
        std::strcat(result, current_kmer);
        std::strcat(result, ")_");
        std::strcat(result, std::to_string(n_found + 1).c_str());
        std::strcat(result, " ");
        rest_length -= k;
    }

    std::strcat(result, seq.substr(rest_start, rest_length).c_str());
    return max_n_found + 1;
}

/// @brief iterates over the frames of a sequence and writes the result to a file.
/// @param seq_name the name of the sequence.
/// @param seq the sequence.
/// @param outfile the file to write the results to.
void interate_over_frames(const std::string &seq_name,
                          const std::string &seq,
                          std::ofstream &outfile,
                          const size_t max_read_len,
                          const int k,
                          const int threshold)
{

    bool only_used_prefix = false;
    if (seq.length() > max_read_len)
    {
        // std::cerr << "The sequence " << seq_name << " is too longer than the specified max_len. "
        //           << "Only using the first max_len bases of the sequence." << std::endl;
        only_used_prefix = true;
    }
    if (seq.length() < k * 2)
    {
        // std::cerr << "The sequence " << seq_name << " is too short to contain a kmer of length " << k << std::endl;
        return;
    }

    size_t result_len = max_read_len * 4 / 3;
    char result[result_len];

    for (int frame = 0; frame < k; ++frame)
    {
        int score = get_repeats(k, seq, frame, result, max_read_len, only_used_prefix);
        if (score >= threshold)
        {
            outfile << seq_name << ", frame: " << frame << " " << result << ", max repeat l:"
                    << score << ", seqlen to long: " << only_used_prefix << std::endl;
        }
    }
}

/// @brief scans a file for sequences and calls interate_over_frames to find repeats for each sequence.
/// @param file_name the name of the file to scan.
/// @param outfile the file to write the results to.
void scan_file(std::string file_name,
               std::ofstream &outfile,
               const size_t max_read_len,
               const int k,
               const int threshold)
{
    // make the kmers parallel or is this too much overhead for to short seqs?
    std::ifstream file(file_name);
    if (!file.is_open())
    {
        std::cerr << "Error: Unable to open file " << file_name << std::endl;
        return;
    }

    std::string line; // maybe not necessary if seqs are always on one line in the input files
    std::string seq; // maybe also convert to char array?
    std::string seq_name;

    bool in_quality = false; // in case it is a fastq file

    while (std::getline(file, line))
    {
        if (line.empty())
            continue;
        if (line[0] == '>' || line[0] == '@')
        {
            if (!seq.empty())
            {
                interate_over_frames(seq_name, seq, outfile, max_read_len, k, threshold);
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
        interate_over_frames(seq_name, seq, outfile, max_read_len, k, threshold);
    }
}

int main(int argc, char *argv[])
{
    Args args = parseArgs(argc, argv);

    // #pragma omp parallel for
    for (auto &file_name : args.files)
    {
        std::ofstream outfile(file_name + ".out");
        if (!outfile.is_open())
        {
            std::cerr << "Error: Unable to open output file " << file_name + ".out" << std::endl;
            continue;
        }
        scan_file(file_name, outfile, args.max_read_len, args.k, args.threshold);
        outfile.close();
    }
    return 0;
}