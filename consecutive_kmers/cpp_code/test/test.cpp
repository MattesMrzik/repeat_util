#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <filesystem>

#include "ConsecutiveKmers.h"

// TODO also CHECK for the score to be correct
// TODO write integration tests for all the file types

std::filesystem::path root_path;

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

TEST_CASE("get_repeat_coordinates: with reverse complement")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string seq = "GATGATCGTGTTGTTGTTGATCCCCCC";
  std::string seq_name = "seq1";
  std::ostringstream string_stream;

  // act
  ck.get_repeat_coordinates(seq_name, seq, string_stream, true);
  std::string output_string = string_stream.str();

  // assert
  CHECK(output_string == "seq1\t0\t5\tATC\nseq1\t8\t18\tAAC\nseq1\t21\t26\tCCC\n");
}

TEST_CASE("get_repeat_coordinates: with reverse complement, no repeats at start or stop")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string seq = "CGATGATCGTGTTGTTGTTGATCCCCCCAA";
  std::string seq_name = "seq1";
  std::ostringstream string_stream;

  // act
  ck.get_repeat_coordinates(seq_name, seq, string_stream, true);
  std::string output_string = string_stream.str();

  // assert
  CHECK(output_string == "seq1\t1\t6\tATC\nseq1\t9\t19\tAAC\nseq1\t22\t27\tCCC\n");
}

TEST_CASE("get_repeat_coordinates: without reverse complement")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string seq = "GATGATCGTGTTGTTGTTGATCCCCCC";
  std::string seq_name = "seq1";
  std::ostringstream string_stream;

  // act
  ck.get_repeat_coordinates(seq_name, seq, string_stream, false);
  std::string output_string = string_stream.str();

  // assert
  CHECK(output_string == "seq1\t0\t5\tATG\nseq1\t8\t18\tGTT\nseq1\t21\t26\tCCC\n");
}

TEST_CASE("get_repeat_coordinates: without reverse complement, no repeats at start or stop")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string seq = "CGATGATCGTGTTGTTGTTGATCCCCCCAA";
  std::string seq_name = "seq1";
  std::ostringstream string_stream;

  // act
  ck.get_repeat_coordinates(seq_name, seq, string_stream, false);
  std::string output_string = string_stream.str();

  // assert
  CHECK(output_string == "seq1\t1\t6\tATG\nseq1\t9\t19\tGTT\nseq1\t22\t27\tCCC\n");
}

TEST_CASE("get_atomic_pattern: with reverse complement")
{
  // arrange
  Args args;
  // args.verbose = true;
  ConsecutiveKmers ck(args);

  // act & assert
  CHECK(ck.get_atomic_pattern("GAT", true) == "ATC");
  CHECK(ck.get_atomic_pattern("AAA", true) == "AAA");
  CHECK(ck.get_atomic_pattern("CTG", true) == "AGC");
  CHECK(ck.get_atomic_pattern("GAA", true) == "AAG");
  CHECK(ck.get_atomic_pattern("TTC", true) == "AAG");
  ck.get_atomic_pattern("TTC", true);
  CHECK(ck.get_atomic_pattern_size() == 5);
}

TEST_CASE("get_atomic_pattern: without reverse complement")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);

  // act & assert
  CHECK(ck.get_atomic_pattern("GAT", false) == "ATG");
  CHECK(ck.get_atomic_pattern("AAA", false) == "AAA");
  CHECK(ck.get_atomic_pattern("GCG", false) == "CGG");
  CHECK(ck.get_atomic_pattern("GCA", false) == "AGC");
  CHECK(ck.get_atomic_pattern_size() == 4);
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

void _write_fastq_file(const std::filesystem::path &path)
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

void _check_fasta_fastq_output(const std::filesystem::path &path)
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
  // TODO maybe also test too long again here, not only in "get_repeats: maximal length"
  CHECK(line == "seq2, frame: 0, (GAT)_2 CGT(GTT)_3 GAT(CCC)_2 (GAT)_2 TGTT, score_type: max, score: 3, seqlen too long: 0");
  std::getline(file, line);
  CHECK(line == "seq2, frame: 1, GATGATCGTG(TTG)_3 ATCCCCCCGATGATTGTT, score_type: max, score: 3, seqlen too long: 0");
  std::getline(file, line);
  CHECK(line == "seq2, frame: 2, GATGATCG(TGT)_3 TGATCCCCCCGATGATTGTT, score_type: max, score: 3, seqlen too long: 0");
  file.close();
}

void compare_files(const std::filesystem::path &correct_file, const std::filesystem::path &out_file)
{
  std::ifstream correct_file_stream(correct_file);
  std::ifstream out_file_stream(out_file);
  std::string line1;
  std::string line2;
  bool correct_file_stream_has_next = bool(std::getline(correct_file_stream, line1));
  bool out_file_stream_has_next = bool(std::getline(out_file_stream, line2));
  while (correct_file_stream_has_next && out_file_stream_has_next)
  {
    CHECK(line1 == line2);
    correct_file_stream_has_next = bool(std::getline(correct_file_stream, line1));
    out_file_stream_has_next = bool(std::getline(out_file_stream, line2));
  }
  CHECK(correct_file_stream.eof() == out_file_stream.eof());
}

TEST_CASE("Integration test: fasta and fastq file")
{
  // arrange
  Args args;
  // tmp is used to not write anything to the git repo
  std::filesystem::path tempDir = std::filesystem::temp_directory_path() / "my_temp_dir_for_integration_test";
  std::filesystem::remove_all(tempDir); // TODO this is now safe!
  std::filesystem::create_directory(tempDir);
  std::filesystem::path resources = root_path / "test" / "resources";
  std::filesystem::path fastaFile = resources / "test.fasta";
  std::filesystem::path out_fastaFile = tempDir / "test.fasta.out";
  std::filesystem::path correct_fastaFile_out = resources / "correct_test.fasta.out";
  std::filesystem::path fastqFile = resources / "test.fastq";
  std::filesystem::path out_fastqFile = tempDir / "test.fastq.out";
  std::filesystem::path correct_fastqFile_out = resources / "correct_test.fastq.out";
  args.output_dir = tempDir;
  args.threshold = 1;
  ConsecutiveKmers ck(args);

  // act
  ck.scan_fasta_and_fastq(fastaFile.string());
  ck.scan_fasta_and_fastq(fastqFile.string());

  // assert
  CHECK(std::filesystem::exists(out_fastaFile));
  CHECK(std::filesystem::exists(out_fastqFile));
  compare_files(correct_fastaFile_out, out_fastaFile);
  compare_files(correct_fastqFile_out, out_fastqFile);

  // clean
  std::filesystem::remove_all(tempDir);
}

TEST_CASE("Integration test: gzipped fasta and fastq file")
{
  // arrange
  Args args;
  // tmp is used to not write anything to the git repo
  std::filesystem::path tempDir = std::filesystem::temp_directory_path() / "my_temp_dir_for_integration_test";
  std::filesystem::remove_all(tempDir); // TODO this is now safe!
  std::filesystem::create_directory(tempDir);
  std::filesystem::path resources = root_path / "test" / "resources";
  std::filesystem::path fastaFile = resources / "test.fasta.gz";
  std::filesystem::path out_fastaFile = tempDir / "test.fasta.gz.out";
  std::filesystem::path correct_fastaFile_out = resources / "correct_test.fasta.out";
  std::filesystem::path fastqFile = resources / "test.fastq.gz";
  std::filesystem::path out_fastqFile = tempDir / "test.fastq.gz.out";
  std::filesystem::path correct_fastqFile_out = resources / "correct_test.fastq.out";
  args.output_dir = tempDir;
  args.threshold = 1;
  ConsecutiveKmers ck(args);

  // act
  ck.scan_fasta_fastq_gz(fastaFile.string());
  ck.scan_fasta_fastq_gz(fastqFile.string());

  // assert
  CHECK(std::filesystem::exists(out_fastaFile));
  CHECK(std::filesystem::exists(out_fastqFile));
  compare_files(correct_fastaFile_out, out_fastaFile);
  compare_files(correct_fastqFile_out, out_fastqFile);

  // clean
  std::filesystem::remove_all(tempDir);
}

TEST_CASE("Integration test: bam file")
{
  // arrange
  Args args;
  // tmp is used to not write anything to the git repo
  std::filesystem::path tempDir = std::filesystem::temp_directory_path() / "my_temp_dir_for_integration_test";
  std::filesystem::remove_all(tempDir); // TODO this is now safe!
  std::filesystem::create_directory(tempDir);
  std::filesystem::path resources = root_path / "test" / "resources";
  std::filesystem::path bamFile = resources / "small.bam";
  std::filesystem::path out_bamFile = tempDir / "small.bam.out";
  std::filesystem::path correct_bamFile_out = resources / "correct_small.bam.out";
  args.output_dir = tempDir;
  args.threshold = 1;
  ConsecutiveKmers ck(args);

  // act
  ck.scan_bam(bamFile.string());

  // assert
  CHECK(std::filesystem::exists(out_bamFile));
  compare_files(correct_bamFile_out, out_bamFile);

  // clean
  std::filesystem::remove_all(tempDir);
}

int main(int argc, char *argv[])
{
  doctest::Context context(argc, argv);
  Args args = parseArgs(argc, argv, true);
  root_path = args.root_path;
  // TODO maybe instead make option to pass path to resources
  auto result = context.run();
  return result;
}