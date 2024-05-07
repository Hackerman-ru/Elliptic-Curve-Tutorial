#include "random.h"

#include "utils/csprng/csprng.hpp"

#include <random>

namespace ECG {
    static constexpr size_t c_size = sizeof(uint) / sizeof(uint32_t);

    uint generate_random_uint() {
        duthomhas::csprng rng;
        std::seed_seq sseq {228};
        rng.seed(sseq);
        uint result = 0;
        uint32_t x = rng(uint32_t());
        std::vector<uint32_t> values = rng(std::vector<uint32_t>(c_size));

        for (size_t i = 0; i < c_size; ++i) {
            result <<= 32;
            result += values[i];
        }

        return result;
    }

    uint generate_random_uint_modulo(const uint& modulus) {
        return generate_random_uint() % modulus;
    }

    uint generate_random_non_zero_uint_modulo(const uint& modulus) {
        return generate_random_uint_modulo(modulus - 1) + 1;
    }
}   // namespace ECG
