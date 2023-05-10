#pragma once
#include <iostream>
#include <bitset>
#include <cstdint>
#include <string>
#include <vector>

// Reflection on MarkerByte data
#define EXPAND(...) __VA_ARGS__
#define MARKERBYTE_MAPPING_MACRO \
  MRKMAP(SKIP, 0x00)             \
  MRKMAP(SOC, 0xB0)              \
  MRKMAP(HOC, 0xB1)              \
  MRKMAP(COD, 0xB2)              \
  MRKMAP(COC, 0xB3)              \
  MRKMAP(QCD, 0xB4)              \
  MRKMAP(QCC, 0xB5)              \
  MRKMAP(TPM, 0xB6)              \
  MRKMAP(CPM, 0xB7)              \
  MRKMAP(SOT, 0xB8)              \
  MRKMAP(STC, 0xB9)              \
  MRKMAP(SOB, 0xBA)              \
  MRKMAP(EOC, 0xBB)              \
  MRKMAP(NIL, 0xFF)

// Last byte of codestream markers
#define MRKMAP(enuminst, val) enuminst = val,
enum class MarkerByte : uint8_t { EXPAND(MARKERBYTE_MAPPING_MACRO) };
#undef MRKMAP

// Bytes in marker
using MarkerType = uint32_t;

constexpr bool INTRINSICS_MSVC =
#ifdef _MSC_VER
    true;
  #include <intrin.h>
#else
    false;
#endif
using btype = uint64_t;
// Plot a binary array, bit by bit
void plot_aligned_binary_array(const btype* ptr, int64_t count, bool printnewline = true);

// Shifts a 128-bit quantity to the left, represented as two 64-bit quantities
inline uint64_t shiftleft128(const uint64_t first, const uint64_t second, const unsigned char bitcount) {
  constexpr unsigned char COUNT = 64;
  if (bitcount == COUNT) return first;
#ifdef _MSC_VER
  return __shiftleft128(first, second, bitcount);
#else
  else
    return (second << bitcount) | (first >> (static_cast<int>(COUNT) - static_cast<int>(bitcount)));
#endif
}

// Reverses all bytes in an unsigned integer value
template <typename T>
inline T byteswap(T x) {
  static_assert(std::is_unsigned_v<T>);
  if constexpr (std::is_same_v<T, uint8_t>)
    return x;
  else if constexpr (INTRINSICS_MSVC) {
    if constexpr (std::is_same_v<T, uint16_t>)
      return _byteswap_ushort(x);
    else if constexpr (std::is_same_v<T, uint32_t>)
      return _byteswap_ulong(x);
    else if constexpr (std::is_same_v<T, uint64_t>)
      return _byteswap_uint64(x);
  } else {
    if constexpr (std::is_same_v<T, uint16_t>)
      return (x << 8) | (x >> 8);
    else if constexpr (std::is_same_v<T, uint32_t>) {
      x = (x & 0x0000FFFF) << 16 | (x & 0xFFFF0000) >> 16;
      return (x & 0x00FF00FF) << 8 | (x & 0xFF00FF00) >> 8;
    } else if constexpr (std::is_same_v<T, uint64_t>) {
      x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32;
      x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16;
      return (x & 0x00FF00FF00FF00FF) << 8 | (x & 0xFF00FF00FF00FF00) >> 8;
    }
  }
}

// Returns the bit position of the most significant bit which is 1
template <typename T>
inline unsigned long revbitscan(T x) {
  static_assert(std::is_unsigned_v<T>);
  unsigned long bitcount;
  if constexpr (INTRINSICS_MSVC && std::is_same_v<T, uint32_t>)
    _BitScanReverse(&bitcount, x);
  else if constexpr (INTRINSICS_MSVC && std::is_same_v<T, uint64_t>)
    _BitScanReverse64(&bitcount, x);
  else {
    bitcount = -1;
    for (auto y = x; y > 0; y >>= 1) ++bitcount;
  }
  return bitcount;
}

// General shared properties between all bitstream classes
struct BitstreamProps {
  using internaltype       = uint64_t;
  using signedinternaltype = std::make_signed<internaltype>::type;
  using baseint            = int32_t;
  static_assert(sizeof(internaltype) / sizeof(baseint) == 2);
  static constexpr int packedUnitBit = CHAR_BIT * sizeof(internaltype);
};

// Data types handling bit stream fractions
struct BitstreamObject : BitstreamProps {
 protected:
  internaltype tail = 0;  // incomplete last word of stream
  int remaining;          // remaining free bits

  BitstreamObject(int remaining_ = 0) : remaining(remaining_) {}
};

// Writes buffered binary stream to a vector in memory
class BitstreamEncoder : protected BitstreamObject {
  // commit tail, update states
  void commit_buffer();

  // shift word within tail
  void doubleshiftleft(internaltype word, int bitcount);

 public:
  // Don't use outside of get_state/move assignment operator
  std::vector<internaltype> stream;

  // Returns number of written bits
  inline size_t bitsize() const { return packedUnitBit * (stream.size() + 1) - remaining; }

  // Returns number of written bytes, if stream would get terminated now
  inline size_t bytesize_oracle() const {
    return sizeof(internaltype) * stream.size() + (packedUnitBit + CHAR_BIT - 1 - remaining) / CHAR_BIT;
  }

  inline bool isempty() const { return (tail == 0) && (stream.size() == 0); };

  // Initialize bit stream, optionally reserve memory if (approximate) final size is known
  //  optionally also add / concat streams (byte aligned byte streams)
  BitstreamEncoder(size_t reserve_bytes1 = 0, const std::byte* data1 = NULL, size_t reserve_bytes2 = 0,
                   const std::byte* data2 = NULL);

  //// TOOD: Add move constructor
  // BitstreamEncoder(BitstreamEncoder&& bse);

  // Fills pointers and values from current bse for implementation of copy/move instructor
  void get_state(internaltype& tail_, int& remaining) const;

  // Move assignment operator
  BitstreamEncoder& operator=(BitstreamEncoder&& bse) noexcept {
    bse.get_state(tail, remaining);
    stream = std::move(bse.stream);
    return *this;
  };

  // Pushing values
  // Write a single bit
  void write_bit(bool bit);

  // Write a word (only push 'bitcount' LSB of 'word' to bit stream)
  void write_word(internaltype word, int bitcount);

  // Push an Exponential-Golomb-coded integer (32-bit), except for UINT_MAX
  void write_exp_golomb(uint32_t x);

  // Push a signed Exponential-Golomb-coded integer (32-bit), except for INT_MIN
  void write_signed_exp_golomb(baseint x);

  // void prepend_bytes(const std::byte* bytes, const size_t nbytes); // TODO: Implement
  void prepend_word(internaltype word);

  // Termination
  // Pads to nearest byte with zeros, pushes word to stream, returns length in bytes
  size_t terminate_stream();

  // Access data as bytes
  inline auto access_data() const { return (const std::byte*)stream.data(); }

  // Appends data to char vector. Terminates stream. Specify whether escape sequence must be inserted.
  void append_to_stream(std::vector<char>& outbuffer, bool escapesequence = false);

  // Debugging
  // Print bitstream contents
  void plain_print() const;
};

// Process binary stream from buffered object
class BitstreamDecoder : protected BitstreamObject {
  // peek next element
  inline btype peek_next() {
    if (source) {
      // if(nbits_valid_remaining > packedUnitBit)
      return byteswap(*source);
    } else {
      return 0;
    }
  }

  // read next element
  inline void read_next() {
    tail = peek_next();
    if (remaining > -1)
      source++;
    else
      source = NULL;
  }

  // peek general word
  template <typename T>
  T peek_word(int bitcount) {
    assert(bitcount > 0 && bitcount <= packedUnitBit);
    internaltype tmp_tail = tail;
    int tmp_remaining     = remaining;
    if (tmp_remaining == 0) {
      tmp_tail      = peek_next();
      tmp_remaining = packedUnitBit;
    }

    T result = (T)tmp_tail >> (packedUnitBit - bitcount);
    if (bitcount > tmp_remaining) {
      bitcount -= tmp_remaining;
      tmp_remaining = packedUnitBit - bitcount;
      tmp_tail      = peek_next();
      result |= tmp_tail >> tmp_remaining;
    } else
      tmp_remaining -= bitcount;

    assert(bitcount < 65);
    if (bitcount > 64)
      throw std::runtime_error("BitstreamDecoder::peek_word: Too many bits " + std::to_string(bitcount)
                               + " > 64 remaining in internal statistics.");
    if (bitcount < 64) tmp_tail <<= bitcount;
    return result;
  }

  // read general word
  template <typename T>
  T read_word(int bitcount) {
    nbits_valid_remaining -= bitcount;
    assert(nbits_valid_remaining > -1);

    assert(bitcount > 0 && bitcount <= packedUnitBit);
    if (remaining == 0) {
      read_next();
      remaining = packedUnitBit;
    }

    T result = (T)tail >> (packedUnitBit - bitcount);
    if (bitcount > remaining) {
      bitcount -= remaining;
      remaining = packedUnitBit - bitcount;
      read_next();
      result |= tail >> remaining;
    } else
      remaining -= bitcount;

    assert(bitcount < 65);
    if (bitcount > 64)
      throw std::runtime_error("BitstreamDecoder::read_word: Too many bits " + std::to_string(bitcount)
                               + " > 64 remaining in internal statistics.");
    if (bitcount < 64) tail <<= bitcount;
    return result;
  }

 public:
  const internaltype* source;
  int64_t nbits_valid_remaining;  // TODO: Properly implement checking to avoid illegal memory access here
                                  // for e.g. 1byte Bitstream

  // Initialize decoder with buffer (and optionally discard some bits in the front)
  BitstreamDecoder(const internaltype* in, int64_t nbytes_valid, int discardbits = 0);

  // Read a single bit
  inline bool read_bit() {
    assert(--nbits_valid_remaining > -1);
    if (remaining-- == 0) {  // Note: remaining-- is called independently of the success of the check
      remaining = packedUnitBit - 1;  // Avoid remaining < 0
      read_next();
    }
    const bool result = (signedinternaltype)tail < 0;
    tail <<= 1;
    return result;
  }

  inline bool peek_bit() {
    int tmp_remaining     = remaining;
    internaltype tmp_tail = tail;
    if (tmp_remaining-- == 0) {  // Note: remaining-- is called independently of the success of the check
      tmp_tail = peek_next();
    }
    const bool result = (signedinternaltype)tmp_tail < 0;
    tmp_tail <<= 1;
    return result;
  }

  // peek an unsigned word made of 'bitcount' bits
  inline internaltype peek_unsigned_word(int bitcount) { return peek_word<internaltype>(bitcount); };

  // peek a signed word made of 'bitcount' bits
  inline signedinternaltype peek_signed_word(int bitcount) {
    return peek_word<signedinternaltype>(bitcount);
  };

  // Read an unsigned word made of 'bitcount' bits
  inline internaltype read_unsigned_word(int bitcount) { return read_word<internaltype>(bitcount); };

  // Read a signed word made of 'bitcount' bits
  inline signedinternaltype read_signed_word(int bitcount) {
    return read_word<signedinternaltype>(bitcount);
  };

  // Read an Exponential-Golomb-coded integer (32-bit)
  uint32_t read_unsigned_exp_golomb();

  // Read a signed Exponential-Golomb-coded integer (32-bit)
  baseint read_signed_exp_golomb();

  // Reset pointer to point to next unused *byte*, based on remaining bits
  // Useful for directly reaching next new bitstream starting at byte boundaries
  void flush_align_byte();
};

// class BitStream {
//   std::vector<uint8_t> buf;
//   uint8_t tmp;
//   uint8_t bits;
//   int32_t pos;
//
//  public:
//   BitStream() : tmp(0), bits(8), pos(0) {}
//   BitStream(uint8_t *p, size_t size, size_t len) : buf((size)), tmp(0), bits(0), pos(size - 1) {
//     for (int i = 0; i < size; ++i) {
//       buf[i] = p[i];
//     }
//     tmp = buf[size - 1];
//     tmp >>= len;
//     bits += len;
//   }
//   void put_bit(uint8_t b) {
//     if (bits == 0) {
//       bits = 8;
//       buf.push_back(tmp);
//       pos++;
//       tmp = 0;
//     }
//     bits--;
//     tmp = static_cast<uint8_t>(tmp + (b << bits));
//   }
//
//   void flush(std::vector<uint8_t> &out) {
//     if (bits > 0) {
//       buf.push_back(tmp);
//       pos++;
//     }
//     buf.push_back(bits);  // put the length info
//     out = buf;
//   }
//
//   int get_bit() {
//     if (bits == 8) {
//       bits = 0;
//       pos--;
//       tmp = buf[pos];
//     }
//     int b = tmp & 1;
//     tmp >>= 1;
//     bits++;
//     return b;
//   }
// };