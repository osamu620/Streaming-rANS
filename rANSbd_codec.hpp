#pragma once
#include <numeric>
#include <vector>
#include <cstdint>
#include "BitStream.hpp"

#if __GNUC__ || __has_attribute(always_inline)
  #define FORCE_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
  #define FORCE_INLINE __forceinline
#else
  #define FORCE_INLINE inline
#endif

enum class ANS : uint64_t {
  ProbBits  = 12,
  MaxDepth  = (1 << ProbBits),
  RANS_L    = 1ull << 31,
  ProbScale = 1 << ProbBits,
  Mask      = (1ull << ProbBits) - 1
};

#if defined(_MSC_VER)
auto multiply_hi64 = [](uint64_t a, uint64_t b) { return __umulh(a, b); };
#else
auto multiply_hi64 = [](uint64_t a, uint64_t b) { return (uint64_t)(((unsigned __int128)a * b) >> 64); };
#endif

// Encoder symbol description
// This (admittedly odd) selection of parameters was chosen to make
// RansEncPutSymbol as cheap as possible.
struct EncSymbolInfo {
  uint64_t rcpFreq;   // Fixed-point reciprocal frequency
  uint32_t freq;      // Symbol frequency
  uint32_t bias;      // Bias
  uint32_t cmplFreq;  // Complement of frequency: (1 << scale_bits) - freq
  uint32_t rcpShift;  // Reciprocal shift
};

// Decoder symbols are straightforward.
struct DecSymbolInfo {
  uint32_t start;  // Start of range.
  uint32_t freq;   // Symbol frequency.
};

struct SymbolStats {
  int maxbd;
  std::vector<std::vector<int>> freqcoeff;
  std::vector<std::vector<int>> cumulatfreq;
  std::vector<std::vector<uint8_t>> cum2sym;
  SymbolStats(int, int16_t *, size_t) noexcept;
  void countFreqs(int16_t *in, size_t len);
  void calcCumFreqs();
  void calcCumFreqs(std::vector<int> &);
  void normalizeFreqs();
};

class rANSEncode {
  uint64_t state;
  std::vector<EncSymbolInfo> esyms;

 public:
  rANSEncode(SymbolStats &, int);
  void EncodeSymbol(BitstreamEncoder &, int);
  size_t DoneEncoding(BitstreamEncoder &);
};

class rANSDecode {
  bool initflag;
  uint64_t state;
  std::vector<DecSymbolInfo> dsyms;
  int remaining;
  uint64_t tmp;

 public:
  rANSDecode(int);
  void init(BitstreamDecoder &, SymbolStats &, int);
  uint32_t DecodeSymbol(BitstreamDecoder &, std::vector<uint8_t> &);
  bool isinit() { return initflag; }
};
