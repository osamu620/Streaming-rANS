#include <cstdio>
#include "rANS_codec.hpp"

int main() {
  FILE *fp = fopen("../../ri-interleaved-low-nonminus1.raw", "rb");
  std::vector<uint16_t> input;

  int d;
  std::vector<uint32_t> freqs(256, 0);
  while ((d = fgetc(fp)) != EOF) {
    input.push_back(d);
    freqs[d]++;
  }

  symbol_stats stats;
  stats.count_freqs(input);
  stats.calc_cum_freqs();
  stats.normalize_freqs();
  // cumlative->symbol table
  // this is super brute force
  uint8_t cum2sym[prob_scale];
  for (int s = 0; s < 256; s++)
    for (uint32_t i = stats.cum_freqs[s]; i < stats.cum_freqs[s + 1]; i++) cum2sym[i] = s;

  rANSencoder enc(stats);
  for (int i = 0; i < input.size(); ++i) {
    enc.encode_symbol(input[input.size() - i - 1]);
  }
  size_t num_bytes = enc.flush();

  std::vector<uint16_t> rec;
  rANSdecoder dec(enc.outbytes, stats);
  for (int i = 0; i < input.size(); ++i) {
    uint32_t s = dec.decode_symbol(cum2sym);
    rec.push_back((std::uint16_t)s);
  }
  // Verification
  for (size_t i = 0; i < rec.size(); ++i) {
    if (input[i] != rec[i]) {
      printf("ERROR!\n");
      return EXIT_FAILURE;
    }
  }
  printf("Codestream bytes = %zu\nbpp = %f\nInput data are perfectly decoded.\n", num_bytes,
         (float)num_bytes * 8 / (4096 * 4096));
  return EXIT_SUCCESS;
}
