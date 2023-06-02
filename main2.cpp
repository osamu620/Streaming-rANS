//
// Created by OSAMU WATANABE on 2023/05/31.
//
#include <cstdio>
#include <chrono>
#include "rANSbd_codec.hpp"

double get_elapsed_time(std::chrono::steady_clock::time_point st) {
  auto fduration = std::chrono::high_resolution_clock::now() - st;
  auto fcount    = std::chrono::duration_cast<std::chrono::microseconds>(fduration).count();
  double ftime   = static_cast<double>(fcount) / 1000.0;
  return ftime;
}

int main() {
  FILE *fp = fopen("../../dice4K_HIGH.dat", "rb");
  if (fp == nullptr) {
    printf("File not found.\n");
    exit(EXIT_FAILURE);
  }
  const size_t maxbd = 7;  // LOW = 4, MEDIUM = 5, HIGH = 7

  fseek(fp, 0, SEEK_END);
  const size_t length = ftell(fp) / sizeof(int16_t);
  fseek(fp, 0, SEEK_SET);

  auto input_data   = std::make_unique<int16_t[]>(length);
  auto decoded_data = std::make_unique<int16_t[]>(length);
  auto *const input = input_data.get();
  fread(input, sizeof(int16_t), length, fp);

  SymbolStats stats(maxbd, input, length);

  double enc_time, dec_time;

  // create stats for rANS enc/dec
  //  SymbolStats stats;
  //  stats.countFreqs(input, length);
  stats.calcCumFreqs();
  stats.normalizeFreqs();
  // cumlative->symbol table
  // this is super brute force
  uint8_t cum2sym[static_cast<uint32_t>(ANS::MaxDepth)];
  //  for (auto s = 0; s < static_cast<size_t>(ANS::MaxDepth); s++)
  //    for (uint32_t i = stats.cumFreqs[s]; i < stats.cumFreqs[s + 1]; i++) cum2sym[i] = s;

  // Encoding
  auto fstart = std::chrono::high_resolution_clock::now();
  std::vector<BitstreamEncoder> bse(maxbd + 1);
  std::vector<rANSEncode> enc;
  for (int i = 0; i <= maxbd; ++i) {
    enc.emplace_back(stats, i);
  }
  std::vector<int> bdbuf;
  size_t uncompressed = 0;
  for (size_t i = length; i > 0; i -= 3) {
    int bd = input[i - 3];
    uncompressed += bd * 2;
    bdbuf.push_back(bd);
    if (bd > 0) {
      enc[bd].EncodeSymbol(bse[bd], input[i - 1]);
      enc[bd].EncodeSymbol(bse[bd], input[i - 2]);
    }
  }
  std::vector<size_t> num_bytes(maxbd + 1, 0);
  for (int i = 1; i <= maxbd; ++i) {
    num_bytes[i] = enc[i].DoneEncoding(bse[i]);
    bse[i].terminate_stream();
  }
  //  size_t num_bytes = enc.DoneEncoding(bse);

  enc_time = get_elapsed_time(fstart);

  // Reverse compressed stream
  for (int i = 1; i <= maxbd; ++i) {
    std::reverse(bse[i].stream.begin(), bse[i].stream.end());
  }
  // reverse bd info
  std::reverse(bdbuf.begin(), bdbuf.end());

  // Decoding
  fstart = std::chrono::high_resolution_clock::now();
  std::vector<BitstreamDecoder> bsd;  //(&bse.stream[0], num_bytes);
  std::vector<rANSDecode> dec;
  for (int i = 0; i <= maxbd; ++i) {
    bsd.emplace_back(&bse[i].stream[0], num_bytes[i]);
    dec.emplace_back(i);
  }
  for (int i = maxbd; i > 0; --i) {
    dec[i].init(bsd[i], stats, i);
  }

  //  rANSDecode dec(bsd, stats);
  for (int i = 0, j = 0; i < length; i += 3, ++j) {
    int bd          = bdbuf[j];
    decoded_data[i] = bd;
    int32_t s       = -1;
    if (bd > 0) {
      s                   = dec[bd].DecodeSymbol(bsd[bd], stats.cum2sym[bd]);
      decoded_data[i + 1] = static_cast<int16_t>(s);
      s                   = dec[bd].DecodeSymbol(bsd[bd], stats.cum2sym[bd]);
      decoded_data[i + 2] = static_cast<int16_t>(s);
    } else {
      decoded_data[i + 1] = static_cast<int16_t>(s);
      decoded_data[i + 2] = static_cast<int16_t>(s);
    }
  }
  dec_time = get_elapsed_time(fstart);

  // Verification
  for (size_t i = 0; i < length; ++i) {
    if (input_data[i] != decoded_data[i]) {
      printf("ERROR at [%d]\n", i);
      return EXIT_FAILURE;
    }
  }
  double output_bytes = std::reduce(num_bytes.begin(), num_bytes.end());
  printf("\nInput data are perfectly decoded!\n\n");
  printf("Codestream bytes = %zu\nbpp = %f\n", (size_t)output_bytes, output_bytes * 8 / (4096 * 4096));
  double ratio = (double)output_bytes * 8 / uncompressed;
  printf("ratio %f\n", ratio);
  printf("Encoding time %-6.4lf[ms], %f [MB/s]\n", enc_time,
         (double)(length * sizeof(int16_t)) / enc_time / 1000);
  printf("Decoding time %-6.4lf[ms], %f [MB/s]\n", dec_time,
         (double)(length * sizeof(int16_t)) / dec_time / 1000);
  return EXIT_SUCCESS;
}
