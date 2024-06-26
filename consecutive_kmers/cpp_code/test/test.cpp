#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <filesystem>
#include <stdexcept>

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

TEST_CASE("cigar_str_to_array")
{
  // arrange
  std::string cigar_str = "12M2I3D19N7S100000H9P19=5X";
  std::string cigar_str_02 = "3M3I3M4D2M";
  std::vector<uint32_t> cigar_array = {12 << 4 | 0,
                                       2 << 4 | 1,
                                       3 << 4 | 2,
                                       19 << 4 | 3,
                                       7 << 4 | 4,
                                       100000 << 4 | 5,
                                       9 << 4 | 6,
                                       19 << 4 | 7,
                                       5 << 4 | 8};
  std::vector<uint32_t> cigar_array_02 = {3 << 4 | 0,
                                          3 << 4 | 1,
                                          3 << 4 | 0,
                                          4 << 4 | 2,
                                          2 << 4 | 0};

  // act & assert
  CHECK(ConsecutiveKmers::cigar_str_to_vector(cigar_str) == cigar_array);
  CHECK(ConsecutiveKmers::cigar_str_to_vector(cigar_str) == cigar_array);
  CHECK(ConsecutiveKmers::cigar_vector_to_str(cigar_array) == cigar_str);
}

TEST_CASE("cigar_array_to_str")
{
  std::string cigar_str = "12M2I3D19N7S100000H9P19=5X";
  std::vector<uint32_t> cigar_array = {12 << 4 | 0,
                                       2 << 4 | 1,
                                       3 << 4 | 2,
                                       19 << 4 | 3,
                                       7 << 4 | 4,
                                       100000 << 4 | 5,
                                       9 << 4 | 6,
                                       19 << 4 | 7,
                                       5 << 4 | 8};

  // act & assert
  CHECK(ConsecutiveKmers::cigar_vector_to_str(cigar_array) == cigar_str);
  CHECK(ConsecutiveKmers::cigar_array_to_str(cigar_array.data(), 9) == cigar_str);
  CHECK(ConsecutiveKmers::cigar_vector_to_str(cigar_array) == cigar_str);
}

void check_cigar_output(std::vector<uint32_t> result, std::vector<uint32_t> correct)
{
  CHECK(result.size() == correct.size());
  for (size_t i = 0; i < result.size(); i++)
    if (result[i] != correct[i])
    {
      CHECK(result[i] == correct[i]);
      CHECK("the previous check failed: returned vector at position" == std::to_string(i));
    }
}

TEST_CASE("get_aligned_reference_positions")
{
  // arrange
  uint32_t ref_start = 1000;
  // reference   0 1 2       3 4 5 6 7 8 9 10 11 12    13    14 15
  //             | | |       | | |          |  |        |        |
  // read        0 1 2 3 4 5 6 7 8          9 10    11 12 13    14
  uint32_t *cigar = ConsecutiveKmers::cigar_str_to_array("3M3I3M4D2M1D1I1M1I1D1M");
  size_t n_cigar = 11;
  std::vector<uint32_t> true_ref_positions = {0, 1, 2, 3, 3, 3, 3, 4, 5, 10, 11, 13, 13, 14, 15};
  std::vector<uint32_t> read_positions____ = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  for (size_t i = 0; i < true_ref_positions.size(); i++)
  {
    true_ref_positions[i] = true_ref_positions[i] + ref_start;
  }

  // act
  auto result = ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions____);

  // assert
  check_cigar_output(result, true_ref_positions);
}

TEST_CASE("get_aligned_reference_positions with I at beginning and end")
{
  // arrange
  uint32_t ref_start = 1000;
  // reference   0 1 2       3 4 5 6 7 8 9 10 11 12    13    14 15
  //             | | |       | | |          |  |        |        |
  // read      0 1 2 3 4 5 6 7 8 9         10 11    12 13 14    15 16
  uint32_t *cigar = ConsecutiveKmers::cigar_str_to_array("1I3M3I3M4D2M1D1I1M1I1D1M1I");
  size_t n_cigar = 13;
  std::vector<uint32_t> true_ref_positions = {0, 0, 1, 2, 3, 3, 3, 3, 4, 5, 10, 11, 13, 13, 14, 15};
  std::vector<uint32_t> read_positions____ = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  for (size_t i = 0; i < true_ref_positions.size(); i++)
  {
    true_ref_positions[i] = true_ref_positions[i] + ref_start;
  }

  // act
  auto result = ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions____);

  // assert
  check_cigar_output(result, true_ref_positions);
}

TEST_CASE("get_aligned_reference_positions with D at beginning and end")
{
  // arrange
  uint32_t ref_start = 1000;
  // reference   0 1 2 3       4 5 6 7 8 9 10 11 12 13    14    15 16 17
  //               | | |       | | |          |  |         |        |
  // read          0 1 2 3 4 5 6 7 8          9 10     11 12 13    14
  uint32_t *cigar = ConsecutiveKmers::cigar_str_to_array("1D3M3I3M4D2M1D1I1M1I1D1M1D");
  size_t n_cigar = 13;
  std::vector<uint32_t> true_ref_positions = {1, 2, 3, 4, 4, 4, 4, 5, 6, 11, 12, 14, 14, 15, 16};
  std::vector<uint32_t> read_positions____ = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  for (size_t i = 0; i < true_ref_positions.size(); i++)
  {
    true_ref_positions[i] = true_ref_positions[i] + ref_start;
  }

  // act
  auto result = ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions____);

  // assert
  check_cigar_output(result, true_ref_positions);
}

TEST_CASE("get_aligned_reference_positions with gaps")
{
  // arrange
  uint32_t ref_start = 1000;
  // reference   0 1 2 3       4 5 6 7 8 9 10 11 12 13    14    15 16 17
  //               | | |       | | |          |  |         |        |
  // read          0 1 2 3 4 5 6 7 8          9 10     11 12 13    14
  uint32_t *cigar = ConsecutiveKmers::cigar_str_to_array("1D3M3I3M4D2M1D1I1M1I1D1M1D");
  size_t n_cigar = 13;
  std::vector<uint32_t> true_ref_positions = {1, 2, 3, 4, 11, 12, 14, 14, 15, 16};
  std::vector<uint32_t> read_positions____ = {0, 1, 2, 3, 9, 10, 11, 12, 13, 14};
  for (size_t i = 0; i < true_ref_positions.size(); i++)
  {
    true_ref_positions[i] = true_ref_positions[i] + ref_start;
  }

  // act
  auto result = ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions____);

  // assert
  check_cigar_output(result, true_ref_positions);
}

TEST_CASE("get_aligned_reference_positions only 2 positions")
{
  // arrange
  uint32_t ref_start = 1000;
  uint32_t *cigar = ConsecutiveKmers::cigar_str_to_array("4M");
  size_t n_cigar = 13;
  std::vector<uint32_t> true_ref_positions = {1, 3};
  std::vector<uint32_t> read_positions____ = {1, 3};
  for (size_t i = 0; i < true_ref_positions.size(); i++)
  {
    true_ref_positions[i] = true_ref_positions[i] + ref_start;
  }

  // act
  auto result = ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions____);

  // assert
  check_cigar_output(result, true_ref_positions);
}

TEST_CASE("get_aligned_reference_positions not in order in middle")
{
  // arrange
  uint32_t ref_start = 1000;
  uint32_t *cigar = ConsecutiveKmers::cigar_str_to_array("1D3M3I3M4D2M1D1I1M1I1D1M1D");
  size_t n_cigar = 13;
  std::vector<uint32_t> read_positions = {0, 1, 2, 12, 9, 10, 11, 12, 15, 14};

  // act & assert
  CHECK_THROWS(ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions));
  try
  {
    auto result = ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions);
  }
  catch (const std::runtime_error &e)
  {
    CHECK(std::string(e.what()).find("[Error] read_positions are not sorted!") == 0);
  }
}
TEST_CASE("get_aligned_reference_positions not in order at very end")
{
  // arrange
  uint32_t ref_start = 1000;
  uint32_t *cigar = ConsecutiveKmers::cigar_str_to_array("1D3M3I3M4D2M1D1I1M1I1D1M1D");
  size_t n_cigar = 13;
  std::vector<uint32_t> read_positions = {0, 1, 2, 3, 9, 10, 11, 12, 15, 14};

  // act & assert
  CHECK_THROWS(ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions));
  try
  {
    auto result = ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions);
  }
  catch (const std::runtime_error &e)
  {
    CHECK(std::string(e.what()).find("[Error] read_positions are not sorted!") == 0);
  }
}

// TODO what is returned if read_pos is not in the cigar string?

TEST_CASE("get_aligned_reference_positions requested position not in cigar string")
{
  // arrange
  uint32_t ref_start = 1000;
  uint32_t *cigar = ConsecutiveKmers::cigar_str_to_array("10M");
  size_t n_cigar = 1;
  std::vector<uint32_t> read_positions = {0, 1, 2, 11};

  // act & assert
  CHECK_THROWS(ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions));
  try
  {
    auto result = ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions);
  }
  catch (const std::runtime_error &e)
  {
    CHECK(std::string(e.what()).find("[Error] query index is not covered by CIGAR string!") == 0);
  }
}

TEST_CASE("get_aligned_reference_positions not in order at beginning")
{
  // arrange
  uint32_t ref_start = 1000;
  uint32_t *cigar = ConsecutiveKmers::cigar_str_to_array("1D3M3I3M4D2M1D1I1M1I1D1M1D");
  size_t n_cigar = 13;
  std::vector<uint32_t> read_positions = {1, 0, 2, 15, 9, 10, 11, 12, 15, 14};

  // act & assert
  CHECK_THROWS(ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions));
  try
  {
    auto result = ConsecutiveKmers::get_aligned_reference_positions(ref_start, cigar, n_cigar, read_positions);
  }
  catch (const std::runtime_error &e)
  {
    CHECK(std::string(e.what()).find("[Error] read_positions are not sorted!") == 0);
  }
}

TEST_CASE("get_aligned_reference_positions with 2 consecutive I, is this possible?")
{
}

// TODO, these need to be updated to the new output format, which must still be decided
TEST_CASE("write_repeat_coordinates: with reverse complement")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string seq_name = "seq1";
  std::string seq = "GATGATCGTGTTGTTGTTGATCCCCCC";
  std::ostringstream string_stream;
  std::string chrom = "chr1";
  uint32_t start = 100;

  uint32_t *cigar_array = ConsecutiveKmers::cigar_str_to_array("27M");
  uint32_t n_cigar = 1;
  bool reverse_complement = true;

  // act
  ck.write_repeat_coordinates(seq_name,
                              seq,
                              string_stream,
                              chrom,
                              start,
                              cigar_array,
                              n_cigar,
                              reverse_complement);
  std::string output_string = string_stream.str();

  // assert
  std::string result = "seq1\tchr1\t0\t5\t100\t105\tATC\n";
  result += "seq1\tchr1\t8\t18\t108\t118\tAAC\n";
  result += "seq1\tchr1\t21\t26\t121\t126\tCCC\n";
  CHECK(output_string == result);
}

TEST_CASE("write_repeat_coordinates: with reverse complement, no repeats at start or stop")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string seq_name = "seq1";
  std::string seq = "CGATGATCGTGTTGTTGTTGATCCCCCCAA";
  std::ostringstream string_stream;
  std::string chrom = "chr1";
  uint32_t start = 100;
  uint32_t *cigar_array = ConsecutiveKmers::cigar_str_to_array("22M3I5M");
  uint32_t n_cigar = 3;
  bool reverse_complement = true;

  // act
  ck.write_repeat_coordinates(seq_name,
                              seq,
                              string_stream,
                              chrom,
                              start,
                              cigar_array,
                              n_cigar,
                              reverse_complement);
  std::string output_string = string_stream.str();

  // assert
  std::string result = "seq1\tchr1\t1\t6\t101\t106\tATC\n";
  result += "seq1\tchr1\t9\t19\t109\t119\tAAC\n";
  result += "seq1\tchr1\t22\t27\t122\t124\tCCC\n";
  CHECK(output_string == result);
}

TEST_CASE("write_repeat_coordinates: without reverse complement")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string seq_name = "seq1";
  std::string seq = "GATGATCGTGTTGTTGTTGATCCCCCC";
  std::ostringstream string_stream;
  std::string chrom = "chr1";
  uint32_t start = 100;
  uint32_t *cigar_array = ConsecutiveKmers::cigar_str_to_array("27M");
  uint32_t n_cigar = 1;
  bool reverse_complement = false;

  // act
  ck.write_repeat_coordinates(seq_name,
                              seq,
                              string_stream,
                              chrom,
                              start,
                              cigar_array,
                              n_cigar,
                              reverse_complement);
  std::string output_string = string_stream.str();

  // assert
  std::string result = "seq1\tchr1\t0\t5\t100\t105\tATG\n";
  result += "seq1\tchr1\t8\t18\t108\t118\tGTT\n";
  result += "seq1\tchr1\t21\t26\t121\t126\tCCC\n";
  CHECK(output_string == result);
}

TEST_CASE("write_repeat_coordinates: without reverse complement, no repeats at start or stop")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string seq_name = "seq1";
  std::string seq = "CGATGATCGTGTTGTTGTTGATCCCCCCAA";
  std::ostringstream string_stream;
  std::string chrom = "chr1";
  uint32_t start = 100;
  uint32_t *cigar_array = ConsecutiveKmers::cigar_str_to_array("30M");
  uint32_t n_cigar = 1;
  bool reverse_complement = false;

  // act
  ck.write_repeat_coordinates(seq_name,
                              seq,
                              string_stream,
                              chrom,
                              start,
                              cigar_array,
                              n_cigar,
                              reverse_complement);
  std::string output_string = string_stream.str();

  // assert
  std::string result = "seq1\tchr1\t1\t6\t101\t106\tATG\n";
  result += "seq1\tchr1\t9\t19\t109\t119\tGTT\n";
  result += "seq1\tchr1\t22\t27\t122\t127\tCCC\n";
  CHECK(output_string == result);
}

TEST_CASE("write_repeat_coordinates, manual seq")
{
  // arrange
  Args args;
  ConsecutiveKmers ck(args);
  std::string seq_name = "seq1";
  std::string seq = "ATGTTATGTAACGCCTGACACAACTCATCGTATTGCTACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCAAGCGCGGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCATCTAATTTGTCCGGGAAGCTGAGTTCAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTGCTTGCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCAAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCGGTTGCGATCAGTTCACTCGTGCACCCAACTGATCTTTCAGCATCTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCATCAAGCGGATACATATTTGAATGTATTAGAAAAATAAACAAATAGGGGTTTCACGTTCCGAAAAGTGCCACCTTAAATTGTAAGCGTTAATATTTTGTTAAAATTCGCGTTAAATTTGTCTAAATCAGCTCATTTTTTAACCGATGAGCAAATCGGCAAAATCCCTTATAAATCAAAAGAATAGACCGAGATAGGGTTGAGTGTTGTTCCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGGCCCACTACGTGAACCATCACCCTAATCAGGTTTCAGTGAGGTAGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAAAGAAAGAAGCGAAAGGAGCGGGCGCTAGGGCGCTGGCAAGTGTACGGTCACGCTGCGCGTAACCACCACACCCGCCGCGCTTAATGCGCCGCTACAGGGCGCGTCCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCAATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAGCGCAAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAGCGCGCGTAATACGACTCACTATAGGGCGAATTGGGTACCCGGCAAAATCGTTAACTAGACGTGGATAACGCAGGAAAGGCGAAACTCGGACAAGACCCTCAGCAATGTAGACAAGTGAATAAGGGCGACCCTCGAGGTCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCCCGGGGCCGGGGCCGGGCCGGGGCCGCGCTCTCGAGGACTACAAAGACCACGACGGAGATTACAAATACACGACATCGATAAGCTTGATATCGAATTCCTGCAGCCCGGGGGATCCACTAGTTCTAGAGCGGCCGCCACCGCGTTCCAGCTTTTTGTTCCCTTTTTAAGTTTGAGGGTTAATTGCGCGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACAGAGCCGGAAGCATAAAGTGTAAACGGCTGGGGTGCCTAATGAGTGAGCTAAACCAGCAGTTTAATTCAGTGCGCTCACTGCCGCTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAGAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTTCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCACCTTCTCCCTTCGGGAACCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACGAGTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTTGGTAGCTCTTGATCCGGCAAACGAACCGCCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAATTTTTGTCAGATCCTTTGATATCTTTTCTACGGGGTCTGACGCTCAGTGGAACAAAACTCACGTTAAGGATTGGTCATGAGATTATCAAAAAGGATCTTCACCTCAATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTCATGCTGACGGTTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCAATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACAATTCAGGATTACCATCTGGCCCCAGTGCTGCAATGATACCGAGCAATACG";
  std::ostringstream string_stream;
  std::string chrom = "chr1";
  uint32_t start = 0;
  std::string cigar_str = "872S4M1I41M2D5M1I29M2D172M4D4M1D5M2D129M1D220M3I41M82I67M4I19M34I45M1D74M1D2M5D6M1I13M3I2M2I88M1I23M1I31M1I2M2I6M1D14M1D4M1D421M1D4M1D235M1D3M1D90M1I79M2D12M4D10M2I33M1D14M1D3M2D93M2I2M2D7M1I102M4D33M9S";
  uint32_t *cigar_array = ConsecutiveKmers::cigar_str_to_array(cigar_str);
  uint32_t n_cigar = ConsecutiveKmers::n_cigar_from_str(cigar_str);
  bool reverse_complement = false;

  // act
  ck.write_repeat_coordinates(seq_name,
                              seq,
                              string_stream,
                              chrom,
                              start,
                              cigar_array,
                              n_cigar,
                              reverse_complement);
  std::string output_string = string_stream.str();

  // out
  // std::cout << string_stream.str() << std::endl;


}

TEST_CASE("n_cigar_from_str")
{
  // arrange & act & assert
  CHECK(ConsecutiveKmers::n_cigar_from_str("1N1M1H1P1S1D1H1=1X") == 9);
  CHECK(ConsecutiveKmers::n_cigar_from_str("11") == 0);
  CHECK(ConsecutiveKmers::n_cigar_from_str("1N1M1H1P1S1D1H") == 7);
  CHECK(ConsecutiveKmers::n_cigar_from_str("144N1M121H14121P41S1D1H") == 7);
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

TEST_CASE("Integration test: bam file, not as coords")
{
  // arrange
  Args args;
  args.bam_output_as_coords = false;
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

TEST_CASE("Integration test: bam file with coords")
{
  // arrange
  Args args;
  args.bam_output_as_coords = true;
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
  // TODO
  // TODO
  // TODO
  // TODO
  //  compare_files(correct_bamFile_out, out_bamFile);

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