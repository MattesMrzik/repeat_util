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

Args args;


/// @brief appends the repeat (kmer)_repeatSize to the result an updates the scores.
/// @param result the result string.
/// @param current_kmer the current kmer.
/// @param n_found the number of consecutive repeats found.
/// @param acc_score the accumulated score.
/// @param max_acc_score the maximum accumulated score.
/// @param kmer_score the kmer score.
/// @param max_kmer_score the maximum kmer score.
void append_repeat(char *result,
                   char *current_kmer,
                   size_t &n_found,
                   size_t &acc_score,
                   size_t &max_acc_score,
                   size_t &kmer_score,
                   size_t &max_kmer_score)
{
    std::strcat(result, "(");
    std::strcat(result, current_kmer);
    std::strcat(result, ")_");
    std::strcat(result, std::to_string(n_found + 1).c_str());
    std::strcat(result, " ");
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
size_t get_repeats(const std::string &seq,
                   int frame,
                   char *result,
                   bool only_used_prefix = false)
{
    size_t l = seq.length();
    size_t k = args.k;
    if (only_used_prefix)
    {
        l = args.max_read_len;
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
                append_repeat(result, current_kmer, n_found, acc_score, max_acc_score, kmer_score, max_kmer_score);
            }
            else
            {
                std::strcat(result, current_kmer);
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

    std::strcat(result, seq.substr(rest_start, rest_length).c_str());

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
void interate_over_frames(const std::string &seq_name,
                          const std::string &seq,
                          std::ofstream &outfile)
{

    bool only_used_prefix = false;
    if (seq.length() > args.max_read_len)
    {
        // std::cerr << "The sequence " << seq_name << " is too longer than the specified max_len. "
        //           << "Only using the first max_len bases of the sequence." << std::endl;
        only_used_prefix = true;
    }
    if (seq.length() < args.k * 2)
    {
        // std::cerr << "The sequence " << seq_name << " is too short to contain a kmer of length " << k << std::endl;
        return;
    }

    size_t result_len = args.max_read_len * 4 / 3;
    char result[result_len];

    for (int frame = 0; frame < args.k; ++frame)
    {
        int score = get_repeats(seq, frame, result, only_used_prefix);
        if (score >= args.threshold)
        {
            outfile << seq_name << ", frame: " << frame << " " << result << ", score " << args.score << ": "
                    << score << ", seqlen to long: " << only_used_prefix << std::endl;
        }
    }
}

/// @brief scans a file for sequences and calls interate_over_frames to find repeats for each sequence.
/// @param file_name the name of the file to scan.
/// @param outfile the file to write the results to.
void scan_file(std::string file_name)
{
    std::ofstream outfile(file_name + ".out");
    if (!outfile.is_open())
    {
        std::cerr << "Error: Unable to open output file " << file_name + ".out" << std::endl;
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
                if (count % 10000 == 0)
                {
                    std::cout << "Scanning sequence " << count << " of 32025228, " << count * 1.0f / 32025228 << std::endl;
                }
                interate_over_frames(seq_name, seq, outfile);
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
        interate_over_frames(seq_name, seq, outfile);
    }
    outfile.close();
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