/* Copyright 2021 Vrije Universiteit Brussel - imec
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided
 * that the following conditions are met:
 *
 * 1. The software is used for research and standardization purposes within SC29/WG1 - JPEG Pleno.
 *
 * 2. Redistributions of source code must retain the above copyright notice, this list of conditions and the
 * following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 4. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or
 * promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 * TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// Description: Bitstream manipulations.

#include "BitStream.hpp"

#include <cassert>

void plot_aligned_binary_array(const btype* ptr, int64_t count, bool printnewline) {
  std::cout << std::bitset<64>(*ptr);
  count -= 64;
  for (; count >= 64; count -= 64) {
    std::cout << std::bitset<64>(*(++ptr));
  }
  for (auto val = *ptr; count > 0; count--) {
    std::cout << ((static_cast<int64_t>(val) < 0) ? '1' : '0');
    val <<= 1;
  }
  if (printnewline) std::cout << '\n';
}

void BitstreamEncoder::commit_buffer() {
  stream.push_back(byteswap(tail));
  tail      = 0;
  remaining = packedUnitBit;
}

void BitstreamEncoder::commit_buffer_noswap() {
  stream.push_back(tail);
  tail      = 0;
  remaining = packedUnitBit;
}

void BitstreamEncoder::doubleshiftleft(internaltype word, int bitcount) {
  static_assert(std::is_same_v<uint64_t, internaltype>);
  tail = shiftleft128(word, tail, bitcount);
}

BitstreamEncoder::BitstreamEncoder(size_t reserve_bytes1, const std::byte* data1, size_t reserve_bytes2,
                                   const std::byte* data2)
    : BitstreamObject(packedUnitBit) {
  if (reserve_bytes1 > 0) {  // data1 == NULL -> plain pre-alloc
    if (data1 == NULL) {
      assert(reserve_bytes2 == 0);  // Sanity check
      stream.reserve(reserve_bytes1 / sizeof(internaltype));
    } else {  // data1 != NULL
      auto byte_copy = [&](const std::byte* data, const size_t nbyte) {
        for (size_t ibyte = 0; ibyte < nbyte;
             ibyte++) {  // TODO: Accelerate for special case of multiples of btype
          this->write_word((internaltype)data[ibyte], sizeof(std::byte) * 8);
        }
      };

      stream.reserve((reserve_bytes1 + reserve_bytes2) / sizeof(internaltype));
      byte_copy(data1, reserve_bytes1);

      if (data2 != NULL) {
        byte_copy(data2, reserve_bytes2);
      }
    }
  } else                          // do nothing
    assert(reserve_bytes2 == 0);  // Sanity check
}

// BitstreamEncoder::BitstreamEncoder(BitstreamEncoder&& bse) : BitstreamObject(packedUnitBit) {
//     // TODO: How to write without public interface
//     //std::vector<internaltype>* stream_ptr;
//     //this->stream = std::move(stream)
// };

void BitstreamEncoder::get_state(internaltype& tail_, int& remaining_) const {
  tail_      = tail;
  remaining_ = remaining;
}

void BitstreamEncoder::write_bit(bool bit) {
  tail = (tail << 1) | (internaltype)bit;
  if (--remaining == 0) commit_buffer();
}

void BitstreamEncoder::write_word(internaltype word, int bitcount) {
  assert(packedUnitBit >= bitcount);
  // push relevant bits up to MSB
  word <<= (internaltype)packedUnitBit - bitcount;

  // if word doesn't fit in tail
  if (bitcount >= remaining) {
    doubleshiftleft(word, remaining);

    bitcount -= remaining;
    word <<= remaining;
    commit_buffer();
  }

  // push remaining bits
  if (bitcount <= 0) return;
  remaining -= bitcount;
  doubleshiftleft(word, bitcount);
}

void BitstreamEncoder::write_word_noswap(BitstreamProps::internaltype word, int bitcount) {
  assert(packedUnitBit >= bitcount);
  // push relevant bits up to MSB
  word <<= (internaltype)packedUnitBit - bitcount;

  // if word doesn't fit in tail
  if (bitcount >= remaining) {
    doubleshiftleft(word, remaining);

    bitcount -= remaining;
    word <<= remaining;
    commit_buffer_noswap();
  }

  // push remaining bits
  if (bitcount <= 0) return;
  remaining -= bitcount;
  doubleshiftleft(word, bitcount);
}

void BitstreamEncoder::write_exp_golomb(uint32_t x) {
  assert(x != std::numeric_limits<uint32_t>::max());
  // push number with (xcount-1) preceding zeros
  write_word((internaltype)x + 1u, 2u * revbitscan((uint32_t)x + 1u) + 1u);
}

void BitstreamEncoder::write_signed_exp_golomb(baseint x) {
  if (x > 0)
    write_exp_golomb(2 * x - 1);
  else
    write_exp_golomb(-2 * x);
}

// void BitstreamEncoder::prepend_words(std::vector<internaltype> words) {
//     stream.reserve(stream.size() + words.size());
//     std::for_each(words.begin(), words.end(), std::back_inserter(stream)); // TODO: Use move operator
//     instead std::rotate(stream.rbegin(), stream.rbegin() + words.size(), stream.rend()); // Rotate last n
//     elements to beginning
// };

// void BitstreamEncoder::prepend_bytes(const std::byte* bytes, const size_t nbytes) {
//     // TODO: Implement
// };

void BitstreamEncoder::prepend_word(internaltype word) {
  stream.push_back(word);
  std::rotate(stream.rbegin(), stream.rbegin() + 1,
              stream.rend());  // Rotate last element to first position
};

void BitstreamEncoder::plain_print() const {
  plot_aligned_binary_array(stream.data(), (int)(packedUnitBit * stream.size()), false);
  const btype shifttail = tail << remaining;
  plot_aligned_binary_array(&shifttail, (size_t)std::max(packedUnitBit - remaining, 0));
}

size_t BitstreamEncoder::terminate_stream() {
  const size_t bytelength = bytesize_oracle();
  // If no bits remain in tail tail will stay 0 --> commit_buffer will be NOOP
  assert(remaining < 65);
  if (remaining > 64)
    throw std::runtime_error("BitstreamEncoder::terminate_stream: Too many bits "
                             + std::to_string(remaining) + " > 64 remaining in internal statistics.");
  if (remaining < 64) tail <<= remaining;
  commit_buffer();
  return bytelength;
}

size_t BitstreamEncoder::terminate_stream_noswap() {
  const size_t bytelength = bytesize_oracle();
  // If no bits remain in tail tail will stay 0 --> commit_buffer will be NOOP
  assert(remaining < 65);
  if (remaining > 64)
    throw std::runtime_error("BitstreamEncoder::terminate_stream: Too many bits "
                             + std::to_string(remaining) + " > 64 remaining in internal statistics.");
  if (remaining < 64) tail <<= remaining;
  commit_buffer_noswap();
  return bytelength;
}

void BitstreamEncoder::append_to_stream(std::vector<char>& outbuffer, bool escapesequence) {
  const auto bytelength = terminate_stream();
  auto arr              = (const uint8_t*)stream.data();
  auto arr_end          = arr + bytelength;
  if constexpr (sizeof(MarkerType) != 4) escapesequence = false;
  auto equals_to_FF = [](uint8_t x) { return x == 0xFF; };

  if (escapesequence) {
    outbuffer.reserve(outbuffer.size() + 2 * bytelength);
    // Insert skip byte 00 after every sequence of >= 3 consecutive FF
    for (auto it = arr;;) {
      auto ff_it = std::find_if(it, arr_end, equals_to_FF);  // find the next FF
      outbuffer.insert(end(outbuffer), it, ff_it);           // copy sequence until that point
      if (ff_it == arr_end) return;                          // if end is reached, terminate
      it = ff_it;                                            // update iterator

      auto nonff_it      = std::find_if_not(it, arr_end, equals_to_FF);  // find the next non-FF symbol
      const auto ffcount = nonff_it - it;                                // count sequence of FF found
      constexpr unsigned char ff = 0xFF;
      outbuffer.insert(end(outbuffer), ffcount, ff);  // push corresponding number of FF to stream
      if (nonff_it == arr_end)
        return;  // stop if end is reached (TODO: push last 0 anyways, or assume marker is next?)
      if (ffcount >= 3)
        outbuffer.push_back(
            (uint8_t)MarkerByte::SKIP);  // Only insert 00 skip byte if at least 3 consecutive FF are found
      it = nonff_it;                     // update iterator
    }
  } else {
    outbuffer.insert(end(outbuffer), arr, arr_end);
  }
}

BitstreamDecoder::BitstreamDecoder(const internaltype* in, int64_t nbytes_valid, int discardbits)
    : source(in), nbits_valid_remaining(nbytes_valid * CHAR_BIT) {
  assert(discardbits >= 0 && discardbits < packedUnitBit);
  if (discardbits > 0) {
    read_next();
    tail <<= discardbits;
    remaining = packedUnitBit - discardbits;
  }
}

void BitstreamDecoder::flush_align_byte() {
  assert(remaining >= 0 && remaining < packedUnitBit);
  const int lastbytecount = remaining / CHAR_BIT;
  nbits_valid_remaining -= remaining;
  assert(nbits_valid_remaining > -1);
  source    = reinterpret_cast<decltype(source)>(reinterpret_cast<const char*>(source) - lastbytecount);
  tail      = 0;
  remaining = 0;
}

uint32_t BitstreamDecoder::read_unsigned_exp_golomb() {
  // Count total number of zeros in run
  int totalcount = 0;
  if (tail == 0 || remaining == 0) {
    totalcount = remaining;
    read_next();
    remaining = packedUnitBit;
    assert(tail != 0);
  }
  // Count number of preceding zeros in tail
  const int zerocount = (packedUnitBit - 1) - revbitscan(tail);
  totalcount += zerocount;

  // remove extraneous zeros (guarantee that zerocount < remaining for valid bitstream)
  tail <<= zerocount;
  remaining -= zerocount;

  // read remainder as word
  assert(totalcount < packedUnitBit / 2);
  return (uint32_t)read_unsigned_word(totalcount + 1) - 1;
}

BitstreamProps::baseint BitstreamDecoder::read_signed_exp_golomb() {
  const auto result = read_unsigned_exp_golomb();
  return (result & 1) ? (result + 1) >> 1 : -(baseint)(result >> 1);
}
