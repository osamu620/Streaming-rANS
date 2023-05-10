#include "rANS_codec.hpp"
symbol_stats::symbol_stats() noexcept : freqs(256, 0), cum_freqs(257, 0) {}
void symbol_stats::count_freqs(std::vector<uint16_t> &input) {
  for (auto &e : input) {
    freqs[e]++;
  }
}

void symbol_stats::calc_cum_freqs() {
  cum_freqs[0] = 0;
  std::partial_sum(freqs.begin(), freqs.end(), cum_freqs.begin() + 1, std::plus<size_t>());
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

rANSencoder::rANSencoder(symbol_stats &stats) {
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

void rANSencoder::encode_symbol(int s) {
  enc_symbol_info *sym = &esyms[s];
  assert(sym->freq != 0);  // can't encode symbol with freq=0

  // renormalize
  uint64_t x     = state;
  uint64_t x_max = (((RANS_L) >> prob_bits) << 32) * sym->freq;  // turns into a shift
  if (x >= x_max) {
    outbytes.push_back((uint32_t)x);
    x >>= 32;
  }

  // x = C(s,x)
  uint64_t q = Rans64MulHi(x, sym->rcp_freq) >> sym->rcp_shift;
  state      = x + sym->bias + q * sym->cmpl_freq;
}

// Flushes the rANS encoder.
size_t rANSencoder::flush() {
  uint64_t x = state;

  outbytes.push_back((uint32_t)(x >> 0));
  outbytes.push_back((uint32_t)(x >> 32));
  size_t num_bytes = outbytes.size() * 4;
  return num_bytes;
}

rANSdecoder::rANSdecoder(std::vector<uint32_t> &c, symbol_stats &stats) : pos(c.size() - 1) {
  codestream = &c[0];
  state      = (uint64_t)(codestream[pos]) << 32;
  state |= (uint64_t)(codestream[pos - 1]) << 0;
  pos -= 2;

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
uint32_t rANSdecoder::decode_symbol(uint8_t *cum2sym) {
  uint32_t s    = cum2sym[state & ((1u << prob_bits) - 1)];
  uint64_t mask = (1ull << prob_bits) - 1;

  // s, x = D(x)
  uint64_t x = state;
  x          = dsyms[s].freq * (x >> prob_bits) + (x & mask) - dsyms[s].start;

  // renormalize
  if (x < RANS_L) {
    x = (x << 32) | codestream[pos];
    pos--;
    assert(x >= RANS_L);
  }

  state = x;
  return s;
}