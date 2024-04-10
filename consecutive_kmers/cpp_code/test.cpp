#include <iostream>
#include <string>
#include "ConsecutiveKmers.h"

int main(int argc, char **argv)
{
  Args args;
  ConsecutiveKmers ck(args);
  std::string seq = "GTAGTAGGG";
  int frame = 0;
  std::string result;
  ck.get_repeats(seq, frame, result);
}
