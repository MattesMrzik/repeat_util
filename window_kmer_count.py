import argparse


def get_repeat(k, seq, frame):
    result = seq[:frame]
    n_found = 0
    for i in range(frame, len(seq), k):
        current_kmer = seq[i: i + k]
        next_kmer = seq[i + k: i + k * 2]
        if current_kmer == next_kmer:
            n_found += 1
        else:
            if n_found > 0:
                result += "(" + current_kmer + ")_" + str(n_found + 1) + " "
                n_found = 0
            else:
                result += current_kmer
    if n_found > 0:
        result += "(" + current_kmer + ")_" + str(n_found + 1) + " "
    return result


def read_sequences(fasta_file, file_format = "fasta"):
    from Bio import SeqIO
    seqs = []
    for record in SeqIO.parse(fasta_file, file_format):
        seqs.append(str(record.seq))
    return seqs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "shows kmer repeats in a sequence. Use only one sequence option at a time")
    parser.add_argument("--sequence", help = "sequence to be analyzed")
    # --seq AGTGACGACGACGTGGAGGGAGTAGTAGTAGCTACGT
    parser.add_argument("--fasta", help = "file with sequence to be analyzed in fasta format")
    parser.add_argument("--fastq", help = "file with sequence to be analyzed in fastq format")
    parser.add_argument("-k", help = "kmer size. if none provided k = 3,4,5 will be used.", type = int, default = -1)

    args = parser.parse_args()

    if not any([args.sequence, args.fasta, args.fastq]):
        print("you must provide a sequence")
        exit()
    if sum([bool(args.sequence), bool(args.fasta), bool(args.fastq)]) > 1:
        print("only one sequence option is allowed")
        exit()

    if args.sequence:
        seqs = [args.sequence]
    elif args.fasta:
        seqs = read_sequences(args.fasta)
    elif args.fastq:
        seqs = read_sequences(args.fastq, "fastq")

    if args.k == -1:
        ks = list(range(3, 6))
    else:
        ks = [args.k]
    print("running for k =", ks)
    for seq in seqs:
        for k in ks:
            for frame in range(k):
                print("the input sequence is:")
                print(seq);
                print(f"the found {k}-mer repeats in frame {frame} are:")
                print(get_repeat(k, seq, frame))
                print()
