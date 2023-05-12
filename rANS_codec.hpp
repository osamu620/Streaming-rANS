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

#define RANS_L 1ull << 31
constexpr uint32_t prob_bits  = 14;
constexpr uint32_t prob_scale = 1 << prob_bits;
constexpr uint64_t mask       = (1ull << prob_bits) - 1;

#if defined(_MSC_VER)
auto multiply_hi64 = [](uint64_t a, uint64_t b) { return __umulh(a, b); };
#else
auto multiply_hi64 = [](uint64_t a, uint64_t b) { return (uint64_t)(((unsigned __int128)a * b) >> 64); };
#endif

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

class rANSEncode {
  uint64_t state;
  enc_symbol_info esyms[256];

 public:
  rANSEncode(symbol_stats &);
  void EncodeSymbol(BitstreamEncoder &, int);
  size_t DoneEncoding(BitstreamEncoder &);
};

class rANSDecode {
  uint64_t state;
  dec_symbol_info dsyms[256];
  int remaining;
  uint64_t tmp;

 public:
  rANSDecode(BitstreamDecoder &, symbol_stats &);
  uint32_t DecodeSymbol(BitstreamDecoder &, uint8_t *cum2sym);
};

symbol_stats::symbol_stats() noexcept : freqs(256, 0), cum_freqs(257, 0) {}
void symbol_stats::count_freqs(std::vector<uint16_t> &input) {
  for (auto &e : input) {
    freqs[e]++;
  }
}

void symbol_stats::calc_cum_freqs() {
  cum_freqs[0] = 0;
  std::partial_sum(freqs.begin(), freqs.end(), cum_freqs.begin() + 1, std::plus<>());
}

void symbol_stats::normalize_freqs() {
  assert(prob_scale >= 256);

  calc_cum_freqs();
  uint32_t cur_total = cum_freqs[256];

  // resample distribution based on cumulative freqs
  for (int i = 1; i <= 256; i++) cum_freqs[i] = ((uint64_t)prob_scale * cum_freqs[i]) / cur_total;

  // if we nuked any non-0 frequency symbol to 0, we need to steal
  // the range to make the frequency nonzero from elsewhere.
  //
  // this is not at all optimal, i'm just doing the first thing that comes to
  // mind.
  for (int i = 0; i < 256; i++) {
    if (freqs[i] && cum_freqs[i + 1] == cum_freqs[i]) {
      // symbol i was set to zero freq

      // find best symbol to steal frequency from (try to steal from low-freq
      // ones)
      uint32_t best_freq = ~0u;
      int best_steal     = -1;
      for (int j = 0; j < 256; j++) {
        uint32_t freq = cum_freqs[j + 1] - cum_freqs[j];
        if (freq > 1 && freq < best_freq) {
          best_freq  = freq;
          best_steal = j;
        }
      }
      assert(best_steal != -1);

      // and steal from it!
      if (best_steal < i) {
        for (int j = best_steal + 1; j <= i; j++) cum_freqs[j]--;
      } else {
        assert(best_steal > i);
        for (int j = i + 1; j <= best_steal; j++) cum_freqs[j]++;
      }
    }
  }

  // calculate updated freqs and make sure we didn't screw anything up
  assert(cum_freqs[0] == 0 && cum_freqs[256] == prob_scale);
  for (int i = 0; i < 256; i++) {
    if (freqs[i] == 0)
      assert(cum_freqs[i + 1] == cum_freqs[i]);
    else
      assert(cum_freqs[i + 1] > cum_freqs[i]);

    // calc updated freq
    freqs[i] = cum_freqs[i + 1] - cum_freqs[i];
  }
}

rANSEncode::rANSEncode(symbol_stats &stats) {
  state = RANS_L;
  enc_symbol_info *s;
  uint32_t start, freq;
  for (int i = 0; i < 256; ++i) {
    s     = &esyms[i];
    start = stats.cum_freqs[i];
    freq  = stats.freqs[i];
    assert(prob_bits <= 31);
    assert(stats.cum_freqs[i] <= (1u << prob_bits));
    assert(stats.freqs[i] <= (1u << prob_bits) - stats.cum_freqs[i]);
    // Say M := 1 << scale_bits.
    //
    // The original encoder does:
    //   x_new = (x/freq)*M + start + (x%freq)
    //
    // The fast encoder does (schematically):
    //   q     = mul_hi(x, rcp_freq) >> rcp_shift   (division)
    //   r     = x - q*freq                         (remainder)
    //   x_new = q*M + bias + r                     (new x)
    // plugging in r into x_new yields:
    //   x_new = bias + x + q*(M - freq)
    //        =: bias + x + q*cmpl_freq             (*)
    //
    // and we can just precompute cmpl_freq. Now we just need to
    // set up our parameters such that the original encoder and
    // the fast encoder agree.

    s->freq      = freq;
    s->cmpl_freq = ((1 << prob_bits) - freq);
    if (freq < 2) {
      // freq=0 symbols are never valid to encode, so it doesn't matter what
      // we set our values to.
      //
      // freq=1 is tricky, since the reciprocal of 1 is 1; unfortunately,
      // our fixed-point reciprocal approximation can only multiply by values
      // smaller than 1.
      //
      // So we use the "next best thing": rcp_freq=~0, rcp_shift=0.
      // This gives:
      //   q = mul_hi(x, rcp_freq) >> rcp_shift
      //     = mul_hi(x, (1<<64) - 1)) >> 0
      //     = floor(x - x/(2^64))
      //     = x - 1 if 1 <= x < 2^64
      // and we know that x>0 (x=0 is never in a valid normalization interval).
      //
      // So we now need to choose the other parameters such that
      //   x_new = x*M + start
      // plug it in:
      //     x*M + start                   (desired result)
      //   = bias + x + q*cmpl_freq        (*)
      //   = bias + x + (x - 1)*(M - 1)    (plug in q=x-1, cmpl_freq)
      //   = bias + 1 + (x - 1)*M
      //   = x*M + (bias + 1 - M)
      //
      // so we have start = bias + 1 - M, or equivalently
      //   bias = start + M - 1.
      s->rcp_freq  = ~0ull;
      s->rcp_shift = 0;
      s->bias      = start + (1 << prob_bits) - 1;
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

      s->rcp_freq  = t0 + (t1 << 32);
      s->rcp_shift = shift - 1;

      // With these values, 'q' is the correct quotient, so we
      // have bias=start.
      s->bias = start;
    }
  }
}

FORCE_INLINE void rANSEncode::EncodeSymbol(BitstreamEncoder &bse, int s) {
  enc_symbol_info *sym = &esyms[s];
  assert(sym->freq != 0);  // can't encode symbol with freq=0

  // renormalize
  uint64_t x     = state;
  uint64_t x_max = (((RANS_L) >> prob_bits) << 32) * sym->freq;  // turns into a shift
  if (x >= x_max) {
    bse.write_word(x, 32);
    x >>= 32;
  }

  // x = C(s,x)
  uint64_t q = multiply_hi64(x, sym->rcp_freq) >> sym->rcp_shift;
  state      = x + sym->bias + q * sym->cmpl_freq;
}

// Flushes the rANS encoder.
size_t rANSEncode::DoneEncoding(BitstreamEncoder &bse) {
  uint64_t x = state;
  bse.write_word(x >> 32, 32);
  bse.write_word(x, 32);

  size_t num_bytes = bse.bytesize_oracle();
  return num_bytes;
}

rANSDecode::rANSDecode(BitstreamDecoder &bsd, symbol_stats &stats) : remaining(2), tmp(0) {
  tmp   = bsd.read_unsigned_word(64);
  state = tmp;
  tmp >>= 64;
  remaining -= 2;

  dec_symbol_info *s;
  uint32_t start, freq;
  for (int i = 0; i < 256; ++i) {
    s     = &dsyms[i];
    start = stats.cum_freqs[i];
    freq  = stats.freqs[i];
    assert(start <= (1 << 31));
    assert(freq <= (1 << 31) - start);
    s->start = start;
    s->freq  = freq;
  }
}
// Returns the current cumulative frequency (map it to a symbol yourself!)
FORCE_INLINE uint32_t rANSDecode::DecodeSymbol(BitstreamDecoder &bsd, uint8_t *cum2sym) {
  uint32_t s = cum2sym[state & ((1u << prob_bits) - 1)];

  // s, x = D(x)
  uint64_t x = state;
  x          = dsyms[s].freq * (x >> prob_bits) + (x & mask) - dsyms[s].start;

  if (remaining == 0 && bsd.nbits_valid_remaining > 0) {
    tmp       = bsd.read_unsigned_word(64);
    remaining = 2;
  }
  // renormalize
  if (x < RANS_L) {
    x = (x << 32) | (tmp & 0xFFFFFFFF);
    tmp >>= 32;
    remaining--;
    assert(x >= RANS_L);
  }

  state = x;
  return s;
}
