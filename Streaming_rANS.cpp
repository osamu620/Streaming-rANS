#include <cstdio>
#include "rANS_codec.hpp"

template <class T>
void show_data(std::vector<T> &data) {
  for (auto &e : data) {
    printf("%d ", e);
  }
  printf("\n");
}

void show_bitstream(std::vector<uint8_t> c) {
  for (int i = 0; i < c.size() - 2; ++i) {
    uint8_t tmp = c[i];
    for (int j = 7; j >= 0; --j) {
      printf("%d", (tmp & (1 << j)) >> j);
    }
  }
  if (c[c.size() - 1]) {
    uint8_t tmp = c[c.size() - 2];
    for (int j = 7; j >= c[c.size() - 1]; --j) {
      printf("%d", (tmp & (1 << j)) >> j);
    }
  }
  printf("\n");
}

int main() {
  // std::vector<uint32_t> F     = {2, 8, 2, 4};
  // std::vector<uint16_t> input = {2, 1, 1, 1, 3, 1, 1, 3, 0, 3, 3, 2, 1, 0, 1, 1};
  // std::vector<uint32_t> F = {95389, 246718, 207649, 83034, 42560, 21342, 6907, 2256,
  //                            57,    53,     19,     6,     1,     2,     0,    7};
  // FILE *fp                = fopen("real-low-nonminus1.raw", "rb");
  // std::vector<uint32_t> F = {191170, 493375, 415324, 165407, 85284, 42871, 13727, 4544,
  //                            128,    96,     43,     18,     2,     2,     0,     9};
  FILE *fp = fopen("ri-interleaved-low-nonminus1.raw", "rb");
  std::vector<uint16_t> input;

  int d;
  std::vector<uint32_t> freqs(256, 0);
  while ((d = fgetc(fp)) != EOF) {
    input.push_back(d);
    freqs[d]++;
  }

  St_rANSencoder enc(freqs);
  for (auto &s : input) {
    enc.encode_symbol(s);
  }
  std::vector<uint8_t> codestream;
  enc.bse.flush(codestream);
  const size_t cs_size = codestream.size() - 1;  // exclude last byte (length info)

  size_t numbits = 8 * (cs_size)-codestream[cs_size];
  printf("\nCompressed bitstream (%lu bits, %f bpp):\n", numbits, (double)numbits / (4096 * 4096));

  St_rANSdecoder dec(enc.get_state(), freqs, &codestream[0], cs_size, codestream[cs_size]);
  std::vector<uint16_t> rec(input.size());
  for (int i = 0; i < input.size(); ++i) {
    rec[rec.size() - 1 - i] = dec.decode_symbol();
  }

  // Verification
  for (int i = 0; i < rec.size(); ++i) {
    if (input[i] != rec[i]) {
      printf("ERROR!\n");
      return EXIT_FAILURE;
    }
  }
  printf("\nInput data are perfectly decoded.\n");
  return EXIT_SUCCESS;
}
