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
  MaxDepth  = (1 << 14),
  RANS_L    = 1ull << 31,
  ProbBits  = 14,
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
  std::vector<uint32_t> freqs;
  std::vector<uint32_t> cumFreqs;
  SymbolStats() noexcept;
  void countFreqs(int16_t *in, size_t len);
  void calcCumFreqs();
  void normalizeFreqs();
};

class rANSEncode {
  uint64_t state;
  EncSymbolInfo esyms[static_cast<size_t>(ANS::MaxDepth)];

 public:
  rANSEncode(SymbolStats &);
  void EncodeSymbol(BitstreamEncoder &, int);
  size_t DoneEncoding(BitstreamEncoder &);
};

class rANSDecode {
  uint64_t state;
  DecSymbolInfo dsyms[static_cast<size_t>(ANS::MaxDepth)];
  int remaining;
  uint64_t tmp;

 public:
  rANSDecode(BitstreamDecoder &, SymbolStats &);
  uint32_t DecodeSymbol(BitstreamDecoder &, uint8_t *cum2sym);
};

SymbolStats::SymbolStats() noexcept
    : freqs(static_cast<size_t>(ANS::MaxDepth), 0), cumFreqs(static_cast<size_t>(ANS::MaxDepth) + 1, 0) {}
void SymbolStats::countFreqs(int16_t *in, size_t len) {
  for (auto i = 0; i < len; ++i) {
    freqs[in[i]]++;
  }
}

void SymbolStats::calcCumFreqs() {
  cumFreqs[0] = 0;
  std::partial_sum(freqs.begin(), freqs.end(), cumFreqs.begin() + 1, std::plus<>());
}

void SymbolStats::normalizeFreqs() {
  assert(static_cast<uint32_t>(ANS::ProbScale) >= static_cast<size_t>(ANS::MaxDepth));

  calcCumFreqs();
  uint32_t cur_total = cumFreqs[static_cast<size_t>(ANS::MaxDepth)];

  // resample distribution based on cumulative freqs
  for (auto i = 1; i <= static_cast<size_t>(ANS::MaxDepth); i++)
    cumFreqs[i] = ((uint64_t) static_cast<uint32_t>(ANS::ProbScale) * cumFreqs[i]) / cur_total;

  // if we nuked any non-0 frequency symbol to 0, we need to steal
  // the range to make the frequency nonzero from elsewhere.
  //
  // this is not at all optimal, i'm just doing the first thing that comes to
  // mind.
  for (auto i = 0; i < static_cast<size_t>(ANS::MaxDepth); i++) {
    if (freqs[i] && cumFreqs[i + 1] == cumFreqs[i]) {
      // symbol i was set to zero freq

      // find best symbol to steal frequency from (try to steal from low-freq
      // ones)
      uint32_t best_freq = ~0u;
      int best_steal     = -1;
      for (auto j = 0; j < static_cast<size_t>(ANS::MaxDepth); j++) {
        uint32_t freq = cumFreqs[j + 1] - cumFreqs[j];
        if (freq > 1 && freq < best_freq) {
          best_freq  = freq;
          best_steal = j;
        }
      }
      assert(best_steal != -1);

      // and steal from it!
      if (best_steal < i) {
        for (int j = best_steal + 1; j <= i; j++) cumFreqs[j]--;
      } else {
        assert(best_steal > i);
        for (int j = i + 1; j <= best_steal; j++) cumFreqs[j]++;
      }
    }
  }

  // calculate updated freqs and make sure we didn't screw anything up
  assert(cumFreqs[0] == 0
         && cumFreqs[static_cast<size_t>(ANS::MaxDepth)] == static_cast<uint32_t>(ANS::ProbScale));
  for (auto i = 0; i < static_cast<size_t>(ANS::MaxDepth); i++) {
    if (freqs[i] == 0)
      assert(cumFreqs[i + 1] == cumFreqs[i]);
    else
      assert(cumFreqs[i + 1] > cumFreqs[i]);

    // calc updated freq
    freqs[i] = cumFreqs[i + 1] - cumFreqs[i];
  }
}

rANSEncode::rANSEncode(SymbolStats &stats) {
  state = static_cast<uint64_t>(ANS::RANS_L);
  EncSymbolInfo *s;
  uint32_t start, freq;
  for (auto i = 0; i < static_cast<size_t>(ANS::MaxDepth); ++i) {
    s     = &esyms[i];
    start = stats.cumFreqs[i];
    freq  = stats.freqs[i];
    assert(static_cast<uint32_t>(ANS::ProbBits) <= 31);
    assert(stats.cumFreqs[i] <= (1u << static_cast<uint32_t>(ANS::ProbBits)));
    assert(stats.freqs[i] <= (1u << static_cast<uint32_t>(ANS::ProbBits)) - stats.cumFreqs[i]);
    // Say M := 1 << scale_bits.
    //
    // The original encoder does:
    //   x_new = (x/freq)*M + start + (x%freq)
    //
    // The fast encoder does (schematically):
    //   q     = mul_hi(x, rcpFreq) >> rcpShift   (division)
    //   r     = x - q*freq                         (remainder)
    //   x_new = q*M + bias + r                     (new x)
    // plugging in r into x_new yields:
    //   x_new = bias + x + q*(M - freq)
    //        =: bias + x + q*cmplFreq             (*)
    //
    // and we can just precompute cmplFreq. Now we just need to
    // set up our parameters such that the original encoder and
    // the fast encoder agree.

    s->freq     = freq;
    s->cmplFreq = ((1 << static_cast<uint32_t>(ANS::ProbBits)) - freq);
    if (freq < 2) {
      // freq=0 symbols are never valid to encode, so it doesn't matter what
      // we set our values to.
      //
      // freq=1 is tricky, since the reciprocal of 1 is 1; unfortunately,
      // our fixed-point reciprocal approximation can only multiply by values
      // smaller than 1.
      //
      // So we use the "next best thing": rcpFreq=~0, rcpShift=0.
      // This gives:
      //   q = mul_hi(x, rcpFreq) >> rcpShift
      //     = mul_hi(x, (1<<64) - 1)) >> 0
      //     = floor(x - x/(2^64))
      //     = x - 1 if 1 <= x < 2^64
      // and we know that x>0 (x=0 is never in a valid normalization interval).
      //
      // So we now need to choose the other parameters such that
      //   x_new = x*M + start
      // plug it in:
      //     x*M + start                   (desired result)
      //   = bias + x + q*cmplFreq        (*)
      //   = bias + x + (x - 1)*(M - 1)    (plug in q=x-1, cmplFreq)
      //   = bias + 1 + (x - 1)*M
      //   = x*M + (bias + 1 - M)
      //
      // so we have start = bias + 1 - M, or equivalently
      //   bias = start + M - 1.
      s->rcpFreq  = ~0ull;
      s->rcpShift = 0;
      s->bias     = start + (1 << static_cast<uint32_t>(ANS::ProbBits)) - 1;
    } else {
      // Alverson, "Integer Division using reciprocals"
      // shift=ceil(log2(freq))
      uint32_t shift = 0;
      uint64_t x0, x1, t0, t1;
      while (freq > (1u << shift)) shift++;

      // long divide ((uint128) (1 << (shift + 63)) + freq-1) / freq
      // by splitting it into two 64:64 bit divides (this works because
      // the dividend has a simple form.)
      x0 = freq - 1;
      x1 = 1ull << (shift + 31);

      t1 = x1 / freq;
      x0 += (x1 % freq) << 32;
      t0 = x0 / freq;

      s->rcpFreq  = t0 + (t1 << 32);
      s->rcpShift = shift - 1;

      // With these values, 'q' is the correct quotient, so we
      // have bias=start.
      s->bias = start;
    }
  }
}

FORCE_INLINE void rANSEncode::EncodeSymbol(BitstreamEncoder &bse, int s) {
  EncSymbolInfo *sym = &esyms[s];
  assert(sym->freq != 0);  // can't encode symbol with freq=0

  // renormalize
  uint64_t x     = state;
  uint64_t x_max = (((static_cast<uint64_t>(ANS::RANS_L)) >> static_cast<uint32_t>(ANS::ProbBits)) << 32)
                   * sym->freq;  // turns into a shift
  if (x >= x_max) {
    bse.write_word(x, 32);
    x >>= 32;
  }

  // x = C(s,x)
  uint64_t q = multiply_hi64(x, sym->rcpFreq) >> sym->rcpShift;
  state      = x + sym->bias + q * sym->cmplFreq;
}

// Flushes the rANS encoder.
size_t rANSEncode::DoneEncoding(BitstreamEncoder &bse) {
  uint64_t x = state;
  bse.write_word(x >> 32, 32);
  bse.write_word(x, 32);

  size_t num_bytes = bse.bytesize_oracle();
  return num_bytes;
}

rANSDecode::rANSDecode(BitstreamDecoder &bsd, SymbolStats &stats) : remaining(2), tmp(0) {
  tmp = bsd.read_unsigned_word(64);
  uint64_t upper, lower;
  if (!(tmp & 0xFFFFFFFF)) {
    tmp >>= 32;
    upper = lower = bsd.read_unsigned_word(64);
    upper >>= 32;
    lower &= 0xFFFFFFFF;
    tmp |= lower << 32;
    state = tmp;
    tmp   = upper;
    remaining--;
    bsd.nbits_valid_remaining += 32;
  } else {
    state = tmp;
    tmp   = 0;
    remaining -= 2;
  }

  DecSymbolInfo *s;
  uint32_t start, freq;
  for (auto i = 0; i < static_cast<size_t>(ANS::MaxDepth); ++i) {
    s     = &dsyms[i];
    start = stats.cumFreqs[i];
    freq  = stats.freqs[i];
    assert(start <= (1 << 31));
    assert(freq <= (1 << 31) - start);
    s->start = start;
    s->freq  = freq;
  }
}
// Returns the current cumulative frequency, map it to a symbol
FORCE_INLINE uint32_t rANSDecode::DecodeSymbol(BitstreamDecoder &bsd, uint8_t *cum2sym) {
  uint32_t s = cum2sym[state & ((1u << static_cast<uint32_t>(ANS::ProbBits)) - 1)];

  // s, x = D(x)
  uint64_t x = state;
  x = dsyms[s].freq * (x >> static_cast<uint32_t>(ANS::ProbBits)) + (x & static_cast<uint64_t>(ANS::Mask))
      - dsyms[s].start;

  // read bitstream, if any
  if (remaining == 0) {
    if (bsd.nbits_valid_remaining > 0) {
      tmp = bsd.read_unsigned_word(32);
      remaining++;
      if (bsd.nbits_valid_remaining > 0) {
        tmp <<= 32;
        tmp |= bsd.read_unsigned_word(32);
        remaining++;
      }
    }
  }

  // renormalize
  if (x < static_cast<uint64_t>(ANS::RANS_L)) {
    x = (x << 32) | (tmp & 0xFFFFFFFF);
    if (tmp >> 32) {
      tmp >>= 32;
    }
    remaining--;
    assert(x >= static_cast<uint64_t>(ANS::RANS_L));
  }

  state = x;
  return s;
}
