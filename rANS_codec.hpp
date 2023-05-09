#pragma once
#include <numeric>
#include <vector>
#include <cstdint>
#include "BitStream.hpp"

#define LOW_LEVEL 2
#define RANGE_FACTOR 8

class St_rANSencoder {
 private:
  size_t total_counts;
  std::vector<size_t> cum_freqs;
  std::vector<uint32_t> freqs;
  int x;  // state
  int range_factor;

 public:
  BitStream bse;
  St_rANSencoder(std::vector<uint32_t> &symbol_counts)
      : total_counts(0),
        cum_freqs((symbol_counts.size() + 1)),
        freqs(symbol_counts),
        range_factor(RANGE_FACTOR),
        bse() {
    total_counts = std::reduce(symbol_counts.begin(), symbol_counts.end());
    // calculate total counts, M
    // for (int i = 0; i < symbol_counts.size(); i++) {
    //   total_counts += symbol_counts[i];
    // }
    x = LOW_LEVEL * total_counts;

    cum_freqs[0] = 0;
    // calculate cumulative sum, C
    // for (int i = 1; i < symbol_counts.size() + 1; i++) {
    //   cumul_counts[i] = cumul_counts[i - 1] + symbol_counts[i - 1];
    // }
    std::partial_sum(symbol_counts.begin(), symbol_counts.end(), cum_freqs.begin() + 1,
                     std::plus<size_t>());
    // normalize_freqs();
  }
  void encode_symbol(uint16_t s) {
    // Output bits to the stream to bring the state in the range for the next encoding
    while (x >= range_factor * freqs[s]) {
      bse.put_bit(x % 2);
      // printf("%d", state % 2);
      x >>= 1;
    }
    rANSencode(s);  // The rANS encoding step
  }

  void rANSencode(int s) { x = (x / freqs[s]) * total_counts + cum_freqs[s] + (x % freqs[s]); }

  int get_state() { return x; }
};

class St_rANSdecoder {
 private:
  size_t total_counts;
  std::vector<size_t> cum_freqs;
  std::vector<uint32_t> freqs;
  int x;  // state
  int range_factor;

 public:
  BitStream bsd;
  St_rANSdecoder(int st, std::vector<uint32_t> &symbol_counts, uint8_t *p, size_t psize, size_t plen)
      : total_counts(0),
        cum_freqs((symbol_counts.size() + 1)),
        freqs(symbol_counts),
        x(st),
        range_factor(RANGE_FACTOR >> 1),
        bsd(p, psize, plen) {
    // calculate total counts, M
    total_counts = std::reduce(symbol_counts.begin(), symbol_counts.end());
    cum_freqs[0] = 0;
    // calculate cumulative sum, C
    std::partial_sum(symbol_counts.begin(), symbol_counts.end(), cum_freqs.begin() + 1,
                     std::plus<size_t>());
  }
  int decode_symbol() {
    // perform the non-streaming rANS decoding
    int slot = x % total_counts;
    int s    = 0;
    for (int i = 0; i < cum_freqs.size(); ++i) {
      if (slot < cum_freqs[i]) {
        s = i - 1;
        break;
      }
    }
    x = (x / total_counts) * freqs[s] + slot - cum_freqs[s];

    // remap the state into the acceptable range
    while (x < range_factor * total_counts) {
      x <<= 1;
      x += bsd.get_bit();
    }
    return s;
  }
};