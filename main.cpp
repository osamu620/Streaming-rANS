#include <cstdio>
#include <chrono>
#include "rANS_codec.hpp"

double get_elapsed_time(std::chrono::steady_clock::time_point st) {
  auto fduration = std::chrono::high_resolution_clock::now() - st;
  auto fcount    = std::chrono::duration_cast<std::chrono::microseconds>(fduration).count();
  double ftime   = static_cast<double>(fcount) / 1000.0;
  return ftime;
}

int main() {
  FILE *fp = fopen("../../ri-interleaved-low-nonminus1.raw", "rb");
  //  FILE *fp = fopen("../../xmax.raw", "rb");
  fseek(fp, 0, SEEK_END);
  const size_t length = ftell(fp) / 20000;
  fseek(fp, 0, SEEK_SET);

  auto input_data   = std::make_unique<int16_t[]>(length);
  auto decoded_data = std::make_unique<int16_t[]>(length);
  auto *const input = input_data.get();
  for (auto i = 0; i < length; ++i) {
    input_data[i] = fgetc(fp);
  }

  double enc_time, dec_time;

  // create stats for rANS enc/dec
  SymbolStats stats;
  stats.countFreqs(input, length);
  stats.calcCumFreqs();
  stats.normalizeFreqs();
  // cumlative->symbol table
  // this is super brute force
  uint8_t cum2sym[static_cast<uint32_t>(ANS::ProbScale)];
  for (auto s = 0; s < static_cast<size_t>(ANS::MaxDepth); s++)
    for (uint32_t i = stats.cumFreqs[s]; i < stats.cumFreqs[s + 1]; i++) cum2sym[i] = s;

  // Encoding
  auto fstart = std::chrono::high_resolution_clock::now();
  BitstreamEncoder bse;
  rANSEncode enc(stats);
  for (int i = 0; i < length; ++i) {
    enc.EncodeSymbol(bse, input[length - i - 1]);  //
  }
  size_t num_bytes = enc.DoneEncoding(bse);
  bse.terminate_stream();
  enc_time = get_elapsed_time(fstart);

  // Reverse compressed stream
  std::reverse(bse.stream.begin(), bse.stream.end());

  // Decoding
  fstart = std::chrono::high_resolution_clock::now();
  BitstreamDecoder bsd(&bse.stream[0], num_bytes);
  rANSDecode dec(bsd, stats);
  for (int i = 0; i < length; ++i) {
    uint32_t s      = dec.DecodeSymbol(bsd, cum2sym);
    decoded_data[i] = static_cast<int16_t>(s);
  }
  dec_time = get_elapsed_time(fstart);

  // Verification
  for (size_t i = 0; i < length; ++i) {
    if (input_data[i] != decoded_data[i]) {
      printf("ERROR at [%d]\n", i);
      return EXIT_FAILURE;
    }
  }
  printf("\nInput data are perfectly decoded!\n\n");
  printf("Codestream bytes = %zu\nbpp = %f\n", num_bytes, (float)num_bytes * 8 / (4096 * 4096));
  printf("Encoding time %-6.4lf[ms], %f [MB/s]\n", enc_time, (double)num_bytes / enc_time / 1000);
  printf("Decoding time %-6.4lf[ms], %f [MB/s]\n", dec_time, (double)num_bytes / dec_time / 1000);
  return EXIT_SUCCESS;
}
