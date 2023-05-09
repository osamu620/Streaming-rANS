#pragma once
#include <cstdint>
#include <vector>

class BitStream {
  std::vector<uint8_t> buf;
  uint8_t tmp;
  uint8_t bits;
  int32_t pos;

 public:
  BitStream() : tmp(0), bits(8), pos(0) {}
  BitStream(uint8_t *p, size_t size, size_t len) : buf((size)), tmp(0), bits(0), pos(size - 1) {
    for (int i = 0; i < size; ++i) {
      buf[i] = p[i];
    }
    tmp = buf[size - 1];
    tmp >>= len;
    bits += len;
  }
  void put_bit(uint8_t b) {
    if (bits == 0) {
      bits = 8;
      buf.push_back(tmp);
      pos++;
      tmp = 0;
    }
    bits--;
    tmp = static_cast<uint8_t>(tmp + (b << bits));
  }

  void flush(std::vector<uint8_t> &out) {
    if (bits > 0) {
      buf.push_back(tmp);
      pos++;
    }
    buf.push_back(bits);  // put the length info
    out = buf;
  }

  int get_bit() {
    if (bits == 8) {
      bits = 0;
      pos--;
      tmp = buf[pos];
    }
    int b = tmp & 1;
    tmp >>= 1;
    bits++;
    return b;
  }
};