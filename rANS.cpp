#include <cstdio>
#include <numeric>
#include <vector>

template <class T>
void show_data(std::vector<T> &data) {
  for (auto &e : data) {
    printf("%d ", e);
  }
  printf("\n");
}

class rANSencoder {
 private:
  size_t total_counts;
  std::vector<size_t> cumul_counts;
  std::vector<uint16_t> frequency;
  int state;

 public:
  rANSencoder(std::vector<uint16_t> &symbol_counts)
      : total_counts(0),
        cumul_counts((symbol_counts.size() + 1)),
        frequency(symbol_counts),
        state(0) {
    total_counts = std::reduce(symbol_counts.begin(), symbol_counts.end());
    // calculate total counts, M
    // for (int i = 0; i < symbol_counts.size(); i++) {
    //   total_counts += symbol_counts[i];
    // }
    cumul_counts[0] = 0;
    // calculate cumulative sum, C
    // for (int i = 1; i < symbol_counts.size() + 1; i++) {
    //   cumul_counts[i] = cumul_counts[i - 1] + symbol_counts[i - 1];
    // }
    std::partial_sum(symbol_counts.begin(), symbol_counts.end(),
                     cumul_counts.begin() + 1, std::plus<size_t>());
  }
  void encode_symbol(int s) {
    state = (state / frequency[s]) * total_counts + cumul_counts[s] +
            (state % frequency[s]);
  }
  int get_state() { return state; }
};

class rANSdecoder {
 private:
  size_t total_counts;
  std::vector<size_t> cumul_counts;
  std::vector<uint16_t> frequency;
  int state;

 public:
  rANSdecoder(int st, std::vector<uint16_t> &symbol_counts)
      : total_counts(0),
        cumul_counts((symbol_counts.size() + 1)),
        frequency(symbol_counts),
        state(st) {
    // calculate total counts, M
    total_counts = std::reduce(symbol_counts.begin(), symbol_counts.end());
    cumul_counts[0] = 0;
    // calculate cumulative sum, C
    std::partial_sum(symbol_counts.begin(), symbol_counts.end(),
                     cumul_counts.begin() + 1, std::plus<size_t>());
  }
  int decode_symbol() {
    int slot = state % total_counts;
    int s;
    for (int i = 0; i < cumul_counts.size(); i++) {
      if (slot < cumul_counts[i]) {
        s = i - 1;
        break;
      }
    }
    state = (state / total_counts) * frequency[s] + slot - cumul_counts[s];
    return s;
  }
};

int main() {
  std::vector<uint16_t> F = {2, 8, 2, 4};
  std::vector<uint16_t> input = {2, 1, 1, 1, 3, 1, 1, 3,
                                 0, 3, 3, 2, 1, 0, 1, 1};

  printf("Input data:\n");
  show_data(input);

  rANSencoder enc(F);
  for (int i = 0; i < input.size(); ++i) {
    enc.encode_symbol(input[i]);
  }
  printf("Final state:%d = %f bits\n\n", enc.get_state(),
         log2(enc.get_state()));

  rANSdecoder dec(enc.get_state(), F);
  std::vector<uint16_t> rec(input.size());
  for (int i = 0; i < input.size(); ++i) {
    rec[rec.size() - 1 - i] = dec.decode_symbol();
  }

  printf("Decoded data:\n");
  show_data(rec);
  for (int i = 0; i < rec.size(); ++i) {
    if (input[i] != rec[i]) {
      printf("ERROR\n");
      return EXIT_FAILURE;
    }
  }
  printf("Input data are perfectly decoded.");
  return EXIT_SUCCESS;
}
