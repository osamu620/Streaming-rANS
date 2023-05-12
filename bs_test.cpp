#include <cstdio>
#include "BitStream.hpp"

enum class ANS { MaxDepth = 10 };
int main() {
  BitstreamEncoder bse;
  bse.write_bit(true);
  bse.write_bit(false);
  bse.terminate_stream();
  printf("%d\n", ANS::MaxDepth);
  return 0;
}