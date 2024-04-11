#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "ConsecutiveKmers.h"

TEST_CASE("get_repeats: repeats at start, middle and end")
{
  Args args;
  ConsecutiveKmers ck(args);
  std::string result;

  std::string seq = "GATGATCGTGTTGTTGTTGATCCCCCC";
  int frame = 0;
  ck.get_repeats(seq, frame, result);
  CHECK(result == "(GAT)_2 CGT(GTT)_3 GAT(CCC)_2 ");

  seq = "AGATGATCGTGTTGTTGTTGATCCCCCC";
  frame = 1;
  ck.get_repeats(seq, frame, result);
  CHECK(result == "A(GAT)_2 CGT(GTT)_3 GAT(CCC)_2 ");

  seq = "AAGATGATCGTGTTGTTGTTGATCCCCCC";
  frame = 2;
  ck.get_repeats(seq, frame, result);
  CHECK(result == "AA(GAT)_2 CGT(GTT)_3 GAT(CCC)_2 ");
}

TEST_CASE("get_repeats: only repeats")
{
  Args args;
  ConsecutiveKmers ck(args);
  std::string result;

  std::string seq = "GATGATGAT";
  int frame = 0;
  ck.get_repeats(seq, frame, result);
  CHECK(result == "(GAT)_3 ");
}

TEST_CASE("get_repeats: maximal length")
{
  Args args;
  args.use_max_seq_len = true;
  args.max_seq_len = 10;

  ConsecutiveKmers ck(args);
  std::string result;

  std::string seq = "AGATGATGATGATCC";
  int frame = 1;
  ck.get_repeats(seq, frame, result);
  CHECK(result == "A(GAT)_3 ");
}
