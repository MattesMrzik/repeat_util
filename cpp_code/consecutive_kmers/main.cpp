#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <omp.h>

std::string get_repeat(int k,
                       const std::string & seq,
                       int frame)
{
    std::string result = seq.substr(0, frame);
    int n_found = 0;
    for (int i = frame; i < seq.length(); i += k)
    {
        std::string current_kmer = seq.substr(i, std::min(k, static_cast<int>(seq.length()) - i));
        std::string next_kmer;
        if (i + k < seq.length()) {
            next_kmer = seq.substr(i + k, std::min(k, static_cast<int>(seq.length()) - i - k));
        } else {
            next_kmer = "";
        }
        if (current_kmer == next_kmer)
        {
            n_found++;
        }
        else
        {
            if (n_found > 0)
            {
                result += "(" + current_kmer + ")_" + std::to_string(n_found + 1) + " ";
                n_found = 0;
            }
            else
            {
                result += current_kmer;
            }
        }
    }
    if (n_found > 0)
    {
        result += "(" + seq.substr(seq.length() - k, k) + ")_" + std::to_string(n_found + 1) + " ";
    }
    return result;
}

std::vector<std::string> read_sequences(const std::string &file_name, const std::string &file_format = "fasta")
{
    std::vector<std::string> seqs;
    std::ifstream file(file_name);
    if (!file.is_open())
    {
        std::cerr << "Error: Unable to open file " << file_name << std::endl;
        return seqs;
    }

    std::string line;
    std::string seq;
    bool read_sequence = false;

    while (std::getline(file, line))
    {
        if (line.empty())
            continue;
        if (line[0] == '>')
        {
            if (!seq.empty())
            {
                seqs.push_back(seq);
                seq.clear();
            }
            if (file_format == "fastq")
                std::getline(file, line); // Skip quality line
            read_sequence = true;
        }
        else if (read_sequence)
        {
            seq += line;
        }
    }
    if (!seq.empty())
    {
        seqs.push_back(seq);
    }
    return seqs;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <file(s)>" << std::endl;
        return 1;
    }

    // #pragma omp parallel for
    for (int i = 1; i < argc; ++i)
    {
        std::string file_name(argv[i]);
        std::vector<std::string> seqs;

        if (file_name.substr(file_name.find_last_of(".") + 1) == "fasta")
        {
            seqs = read_sequences(file_name);
        }
        else if (file_name.substr(file_name.find_last_of(".") + 1) == "fastq")
        {
            seqs = read_sequences(file_name, "fastq");
        }
        else
        {
            std::cerr << "Error: Unsupported file format" << std::endl;
            continue;
        }

        int k = 3; // default k-mer size

        // #pragma omp parallel for
        for (size_t j = 0; j < seqs.size(); ++j)
        {
            for (int k_val = 3; k_val <= k; ++k_val)
            {
                for (int frame = 0; frame < k_val; ++frame)
                {
                    std::cout << "File: " << file_name << ", Sequence: " << j << ", k: " << k_val << ", Frame: " << frame << std::endl;
                    std::cout << "Sequence: " << seqs[j] << std::endl;
                    std::cout << "Found " << k_val << "-mer repeats in frame " << frame << " are:" << std::endl;
                    std::cout << get_repeat(k_val, seqs[j], frame) << std::endl
                              << std::endl;
                }
            }
        }
    }

    return 0;
}
