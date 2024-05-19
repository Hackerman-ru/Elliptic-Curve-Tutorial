#include "random.h"

#include "utils/csprng/csprng.hpp"

#include <random>

namespace elliptic_curve_guide::algorithm::random {
    static constexpr size_t c_size = uint_info::uint_bytes_number / sizeof(uint32_t);

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

    field::FieldElement generate_random_field_element(const field::Field& field) {
        uint result = generate_random_uint_modulo(field.modulus());
        return field.element(result);
    }

    field::FieldElement generate_random_non_zero_field_element(const field::Field& field) {
        uint result = generate_random_non_zero_uint_modulo(field.modulus());
        return field.element(result);
    }
}   // namespace elliptic_curve_guide::algorithm::random
