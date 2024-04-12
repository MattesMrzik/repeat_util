#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <filesystem>

#include "ConsecutiveKmers.h"

// TODO also CHECK for the score to be correct
// TODO write integration tests for all the file types

TEST_CASE("get_repeats: repeats at start, middle and end, frame 0")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string result;
  std::string seq = "GATGATCGTGTTGTTGTTGATCCCCCC";
  int frame = 0;

  // act
  ck.get_repeats(seq, frame, result);

  // assert
  CHECK(result == "(GAT)_2 CGT(GTT)_3 GAT(CCC)_2 ");
  CHECK(seq == ck.expand_collapsed_repeats(result));
}
TEST_CASE("get_repeats: repeats at start, middle and end, frame 1")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string result;
  std::string seq = "AGATGATCGTGTTGTTGTTGATCCCCCC";
  int frame = 1;

  // act
  ck.get_repeats(seq, frame, result);

  // assert
  CHECK(result == "A(GAT)_2 CGT(GTT)_3 GAT(CCC)_2 ");
  CHECK(seq == ck.expand_collapsed_repeats(result));
}
TEST_CASE("get_repeats: repeats at start, middle and end, frame 1")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string result;
  std::string seq = "AAGATGATCGTGTTGTTGTTGATCCCCCC";
  int frame = 2;

  // act
  ck.get_repeats(seq, frame, result);

  // assert
  CHECK(result == "AA(GAT)_2 CGT(GTT)_3 GAT(CCC)_2 ");
  CHECK(seq == ck.expand_collapsed_repeats(result));
}

TEST_CASE("get_repeats: repeats in middle, one base at end")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string result;
  std::string seq = "ATAGATGATCGTGTTGTTGTTGATCCCCCCT";
  int frame = 0;

  // act
  ck.get_repeats(seq, frame, result);

  // assert
  CHECK(result == "ATA(GAT)_2 CGT(GTT)_3 GAT(CCC)_2 T");
  CHECK(seq == ck.expand_collapsed_repeats(result));
}
TEST_CASE("get_repeats: repeats in middle, two bases at end")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string result;
  std::string seq = "ATAGATGATCGTGTTGTTGTTGATCCCCCCTT";
  int frame = 0;

  // act
  ck.get_repeats(seq, frame, result);

  // assert
  CHECK(result == "ATA(GAT)_2 CGT(GTT)_3 GAT(CCC)_2 TT");
  CHECK(seq == ck.expand_collapsed_repeats(result));
}

TEST_CASE("get_repeats: only repeats")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string result;
  std::string seq = "GATGATGAT";
  int frame = 0;

  // act
  ck.get_repeats(seq, frame, result);

  // assert
  CHECK(result == "(GAT)_3 ");
  CHECK(seq == ck.expand_collapsed_repeats(result));
}

TEST_CASE("get_repeats: maximal length")
{
  // arrange
  Args args;
  args.use_max_seq_len = true;
  args.max_seq_len = 11;
  ConsecutiveKmers ck(args);
  std::string result;
  std::string seq = "AGATGATGATGATCC";
  int frame = 1;

  // act
  ck.get_repeats(seq, frame, result);

  // assert
  CHECK(result == "A(GAT)_3 G");
  CHECK(seq == ck.expand_collapsed_repeats(result.append("ATCC")));
}

void write_fasta_file(const std::filesystem::path &path)
{
  std::ofstream file(path);
  file << ">seq1\n";
  file << "GTGATGATGGCG\n";
  file << ">seq2\n";
  file << "GATGATCGTGTTGTTGTTGATCCCCCCG\n";
  file << "ATGATTGTT\n";
  file << ">seq3\n";
  file << "GATTCTGAA\n";
  file.close();
}

void write_fastq_file(const std::filesystem::path &path)
{
  std::ofstream file(path);
  file << "@seq1\n";
  file << "GTGATGATGGCG\n";
  file << "+seq1\n";
  file << "#FFFFFFFFF,F\n";
  file << "@seq2\n";
  file << "GATGATCGTGTTGTTGTTGATCCCCCCG\n";
  file << "ATGATTGTT\n";
  file << "+seq2\n";
  file << "FFFFFF:FFF\n";
  file << "@seq3\n";
  file << "GATTCTGAA\n";
  file << "+seq3\n";
  file << "FFFFFF:FF\n";
  file.close();
}

void check_fasta_fastq_output(const std::filesystem::path &path)
{
  std::ifstream file(path);
  std::string line;
  std::getline(file, line);
  CHECK(line == "seq1, frame: 0, GTG(ATG)_2 GCG, score_type: max, score: 2, seqlen too long: 0");
  std::getline(file, line);
  CHECK(line == "seq1, frame: 1, G(TGA)_2 TGGCG, score_type: max, score: 2, seqlen too long: 0");
  std::getline(file, line);
  CHECK(line == "seq1, frame: 2, GT(GAT)_2 GGCG, score_type: max, score: 2, seqlen too long: 0");
  std::getline(file, line);


  // TODO check this seq
  CHECK(line == "seq2, frame: 0, (GAT)_2 CGT(GTT)_3 GAT(CCC)_2 (GAT)_2 TGTT, score_type: max, score: 3, seqlen too long: 0");
  std::getline(file, line);
  CHECK(line == "seq2, frame: 1, GATGATCGTG(TTG)_3 ATCCCCCCGATGATTGTT, score_type: max, score: 3, seqlen too long: 0");
  std::getline(file, line);
  CHECK(line == "seq2, frame: 2, GATGATCG(TGT)_3 TGATCCCCCCGATGATTGTT, score_type: max, score: 3, seqlen too long: 0");
  file.close();
}

TEST_CASE("Integration test: fasta and fastq file")
{
  // arrange
  Args args;
  std::filesystem::path tempDir = std::filesystem::temp_directory_path() / "my_temp_dir_for_integration_test";
  std::filesystem::remove_all(tempDir); // TODO this is now safe!
  std::filesystem::create_directory(tempDir);
  std::filesystem::path fastaFile = tempDir / "test.fasta";
  std::filesystem::path fastqFile = tempDir / "test.fastq";
  write_fasta_file(fastaFile);
  write_fastq_file(fastqFile);
  args.output_dir = tempDir;
  args.threshold = 1;
  ConsecutiveKmers ck(args);

  // act
  ck.scan_fasta_and_fastq(fastaFile.string());
  ck.scan_fasta_and_fastq(fastqFile.string());

  // assert
  check_fasta_fastq_output(tempDir / "test.fasta.out");
  check_fasta_fastq_output(tempDir / "test.fastq.out");

  // clean
  std::filesystem::remove_all(tempDir);
}
