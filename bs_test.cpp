#include <cstdio>
#include "BitStream.hpp"

int main() {
  BitStream bse;

  bse.put_bit(0);
  bse.put_bit(0);
  bse.put_bit(0);
  bse.put_bit(0);

  bse.put_bit(1);
  bse.put_bit(1);
  bse.put_bit(0);
  bse.put_bit(0);

  bse.put_bit(1);
  bse.put_bit(0);
  bse.put_bit(0);
  bse.put_bit(1);

  bse.put_bit(1);
  bse.put_bit(0);
  bse.put_bit(1);
  bse.put_bit(1);

  bse.put_bit(0);
  bse.put_bit(0);
  bse.put_bit(0);
  bse.put_bit(1);

  bse.put_bit(1);
  bse.put_bit(1);
  bse.put_bit(1);
  bse.put_bit(1);

  bse.put_bit(1);
  bse.put_bit(1);
  bse.put_bit(0);
  bse.put_bit(0);

  std::vector<uint8_t> output;
  bse.flush(output);

  BitStream bsd(&output[0], 4, output[4]);
  int i;
  for (i = 0; i < 4 * 8 - output[4]; ++i) {
    printf("%d", bsd.get_bit());
  }
  printf("\n");
  return 0;
}