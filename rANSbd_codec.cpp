#include "rANSbd_codec.hpp"

SymbolStats::SymbolStats(int m, int16_t *const input_data, const size_t length) noexcept : maxbd(m) {
  uint64_t levels;
  for (size_t i = 0; i <= static_cast<size_t>(maxbd); i++) {
    // Initalize bin
    levels = 1ull << i;  // pow(2, i);
    std::vector<int> bin(levels, 0);
    std::vector<int> cumbin(levels + 1, 0);  // +1 because cumulative bin
    freqcoeff.push_back(bin);
    cumulatfreq.push_back(cumbin);
  }
  //  for (auto i = 0; i < length; i += 3) {
  //    int bd = input_data[i];
  //    if (bd > 0) {
  //      freqcoeff[bd][input_data[i + 1]]++;
  //      freqcoeff[bd][input_data[i + 2]]++;
  //    } else {
  //      freqcoeff[bd][input_data[i + 1] + 1]++;
  //      freqcoeff[bd][input_data[i + 2] + 1]++;
  //    }
  //  }
  auto get_logstic = [](size_t in, double mu, double s) {
    double E = exp(-(in - mu) / s);
    return ceil((E / (s * (1 + E) * (1 + E))) * static_cast<double>(ANS::MaxDepth));
  };
  auto get_normal = [](size_t in, double mu, double sigma) {
    return ceil(static_cast<double>(ANS::MaxDepth) / (sqrt(2 * 3.141592653589793 * sigma * sigma))
                * exp(-((in - mu) * (in - mu)) / (2 * sigma * sigma)));
  };
  freqcoeff[0][0] = 100;
  freqcoeff[1][0] = 99;
  freqcoeff[1][1] = 101;
  double s[8]     = {0, 0, 0.7, 1, 2, 4, 8, 10};
  for (size_t i = 2; i <= maxbd; i++) {
    // Initalize bin
    levels = 1ull << i;  // pow(2, i);
    for (size_t j = 0; j < levels; ++j) {
      int bd          = i;
      double mu       = (1 << (bd - 1)) - 0.5;
      double sigma    = -0.1019 * (bd * bd * bd) + 1.7857 * (bd * bd) - 6.2553 * bd + 7.2381;
      freqcoeff[i][j] = get_logstic(j, mu, s[i]);
    }
  }
}
void SymbolStats::countFreqs(int16_t *in, size_t len) {
  //  for (auto i = 0; i < len; ++i) {
  //    freqs[in[i]]++;
  //  }
}

void SymbolStats::calcCumFreqs() {
  for (size_t i = 0; i <= static_cast<size_t>(maxbd); i++) {
    std::partial_sum(freqcoeff[i].begin(), freqcoeff[i].end(), cumulatfreq[i].begin() + 1, std::plus<>());
  }
}

void SymbolStats::normalizeFreqs() {
  assert(static_cast<uint32_t>(ANS::ProbScale) >= static_cast<size_t>(ANS::MaxDepth));
  for (int b = 0; b <= maxbd; ++b) {
    cum2sym.emplace_back(static_cast<size_t>(ANS::MaxDepth));
    auto cumFreqs      = cumulatfreq[b];
    auto freqs         = freqcoeff[b];
    const size_t len   = freqs.size();
    uint32_t cur_total = cumFreqs[len];

    // resample distribution based on cumulative freqs
    for (auto i = 1; i <= len; i++)
      cumFreqs[i] = ((uint64_t) static_cast<int32_t>(ANS::ProbScale) * cumFreqs[i]) / cur_total;

    // if we nuked any non-0 frequency symbol to 0, we need to steal
    // the range to make the frequency nonzero from elsewhere.
    //
    // this is not at all optimal, i'm just doing the first thing that comes to
    // mind.
    for (auto i = 0; i < len; i++) {
      if (freqs[i] && cumFreqs[i + 1] == cumFreqs[i]) {
        // symbol i was set to zero freq

        // find best symbol to steal frequency from (try to steal from low-freq
        // ones)
        uint32_t best_freq = ~0u;
        int best_steal     = -1;
        for (auto j = 0; j < len; j++) {
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
    assert(cumFreqs[0] == 0 && cumFreqs[len] == static_cast<uint32_t>(ANS::ProbScale));
    for (auto i = 0; i < len; i++) {
      if (freqs[i] == 0)
        assert(cumFreqs[i + 1] == cumFreqs[i]);
      else
        assert(cumFreqs[i + 1] > cumFreqs[i]);

      // calc updated freq
      freqs[i] = cumFreqs[i + 1] - cumFreqs[i];
    }
    cumulatfreq[b] = cumFreqs;
    freqcoeff[b]   = freqs;
    for (int sym = 0; sym < len; sym++)
      for (uint32_t i = cumFreqs[sym]; i < cumFreqs[sym + 1]; i++) cum2sym[b][i] = sym;
  }
}

rANSEncode::rANSEncode(SymbolStats &stats, int bd) {
  state = static_cast<uint64_t>(ANS::RANS_L);
  EncSymbolInfo *s;
  uint32_t start, freq;
  for (auto i = 0; i < stats.freqcoeff[bd].size(); ++i) {
    esyms.emplace_back();
    s     = &esyms[i];
    start = stats.cumulatfreq[bd][i];
    freq  = stats.freqcoeff[bd][i];
    assert(static_cast<uint32_t>(ANS::ProbBits) <= 31);
    assert(start <= (1u << static_cast<uint32_t>(ANS::ProbBits)));
    assert(freq <= (1u << static_cast<uint32_t>(ANS::ProbBits)) - start);
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

void rANSEncode::EncodeSymbol(BitstreamEncoder &bse, int s) {
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

rANSDecode::rANSDecode(int bd) : initflag(false), state(0), remaining(0), tmp(0) {
  for (int i = 0; i < (1u << bd); ++i) {
    dsyms.emplace_back();
  }
}
void rANSDecode::init(BitstreamDecoder &bsd, SymbolStats &stats, int bd) {
  remaining = 2;
  tmp       = 0;

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
  for (auto i = 0; i < (1u << bd); ++i) {
    s     = &dsyms[i];
    start = stats.cumulatfreq[bd][i];
    freq  = stats.freqcoeff[bd][i];
    assert(start <= (1 << 31));
    assert(freq <= (1 << 31) - start);
    s->start = start;
    s->freq  = freq;
  }
  initflag = true;
}
// Returns the current cumulative frequency, map it to a symbol
uint32_t rANSDecode::DecodeSymbol(BitstreamDecoder &bsd, std::vector<uint8_t> &cum2sym) {
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
