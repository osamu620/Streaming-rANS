#pragma once
#include <numeric>
#include <vector>
#include <cstdint>
#include "BitStream.hpp"

#define RANS_L 1ull << 31
constexpr uint32_t prob_bits  = 14;
constexpr uint32_t prob_scale = 1 << prob_bits;

static inline uint64_t Rans64MulHi(uint64_t a, uint64_t b) {
  return (uint64_t)(((unsigned __int128)a * b) >> 64);
}

// Encoder symbol description
// This (admittedly odd) selection of parameters was chosen to make
// RansEncPutSymbol as cheap as possible.
struct enc_symbol_info {
  uint64_t rcp_freq;   // Fixed-point reciprocal frequency
  uint32_t freq;       // Symbol frequency
  uint32_t bias;       // Bias
  uint32_t cmpl_freq;  // Complement of frequency: (1 << scale_bits) - freq
  uint32_t rcp_shift;  // Reciprocal shift
};

// Decoder symbols are straightforward.
struct dec_symbol_info {
  uint32_t start;  // Start of range.
  uint32_t freq;   // Symbol frequency.
};

struct symbol_stats {
  std::vector<uint32_t> freqs;
  std::vector<uint32_t> cum_freqs;
  symbol_stats() noexcept;
  void count_freqs(std::vector<uint16_t> &);
  void calc_cum_freqs();
  void normalize_freqs();
};

class rANSencoder {
  uint64_t state;
  enc_symbol_info esyms[256];

 public:
  std::vector<uint32_t> outbytes;
  rANSencoder(symbol_stats &);
  void encode_symbol(int);
  size_t flush();
};

class rANSdecoder {
  uint64_t state;
  dec_symbol_info dsyms[256];
  uint32_t *codestream;
  size_t pos;

 public:
  rANSdecoder(std::vector<uint32_t> &, symbol_stats &);
  uint32_t decode_symbol(uint8_t *);
};
