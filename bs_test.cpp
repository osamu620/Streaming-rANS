#include <cstdio>
#include "BitStream.hpp"

int main() {
  BitstreamEncoder bse;
  bse.write_bit(true);
  bse.write_bit(false);
  bse.terminate_stream();
  return 0;
}