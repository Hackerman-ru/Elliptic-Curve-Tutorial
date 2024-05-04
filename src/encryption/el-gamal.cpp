#include "el-gamal.h"

#include "utils/bitsize.h"
#include "utils/random.h"

#include <map>

namespace ECG {
    static uint generate_non_zero_uint_modulo(const uint& modulus) {
        return generate_random_uint_modulo(modulus - 1) + 1;
    }

    ElGamal::ElGamal(const Curve& curve, const Point& generator, const uint& generator_order) :
        m_curve(curve), m_generator(generator), m_generator_order(generator_order) {}

    ElGamal::Keys ElGamal::generate_keys() const {
        uint private_key = generate_non_zero_uint_modulo(m_generator_order);
        Point public_key = private_key * m_generator;
        return Keys {.private_key = private_key, .public_key = public_key};
    }

    ElGamal::EncryptedMessage<ElGamal::EncryptionType::Standard>
        ElGamal::encrypt(const uint& message, const Point& public_key) const {
        Point point_message = map_to_curve(message);
        return encrypt(point_message, public_key);
    }

    ElGamal::EncryptedMessage<ElGamal::EncryptionType::Standard>
        ElGamal::encrypt(const Point& message, const Point& public_key) const {
        const uint k = generate_non_zero_uint_modulo(m_generator_order);
        const Point generator_degree = k * m_generator;
        const Point message_with_salt = message + k * public_key;
        return {.generator_degree = generator_degree, .message_with_salt = message_with_salt};
    }

    ElGamal::EncryptedMessage<ElGamal::EncryptionType::Hashed>
        ElGamal::encrypt(const uint& message, const Point& public_key,
                         const std::function<uint(const Point&)>& hash_function) const {
        const uint k = generate_non_zero_uint_modulo(m_generator_order);
        const Point generator_degree = k * m_generator;
        const uint message_with_salt = message ^ hash_function(k * public_key);
        return {.generator_degree = generator_degree, .message_with_salt = message_with_salt};
    }

    ElGamal::Point
        ElGamal::decrypt(const EncryptedMessage<ElGamal::EncryptionType::Standard>& encrypted_message,
                         const uint& private_key) const {
        return encrypted_message.message_with_salt - private_key * encrypted_message.generator_degree;
    }

    uint ElGamal::decrypt(const EncryptedMessage<ElGamal::EncryptionType::Hashed>& encrypted_message,
                          const uint& private_key,
                          const std::function<uint(const Point&)>& hash_function) const {
        return encrypted_message.message_with_salt
             ^ hash_function(private_key * encrypted_message.generator_degree);
    }

    static std::map<uint, uint> p_zero_mask;
    static constexpr uint full_bits = uint(0) - 1;

    ElGamal::Point ElGamal::map_to_curve(const uint& message) const {
        const Field& F = m_curve.get_field();
        const uint& p = F.modulus();

        if (!p_zero_mask.contains(p)) {
            const size_t l = actual_bit_size(p) >> 1;
            uint zero_mask = (full_bits >> l) << l;
            p_zero_mask.insert({p, zero_mask});
        }

        const uint& zero_mask = p_zero_mask.at(p);

        // Should take less than 3 iterations for large p: https://eprint.iacr.org/2013/373.pdf, page 5
        for (;;) {
            uint x = generate_random_uint_modulo(p);
            x &= zero_mask;
            x |= message ^ (message & zero_mask);
            auto opt = m_curve.point_with_x_equal_to(F.element(std::move(x)));

            if (opt.has_value()) {
                return opt.value();
            }
        }
    }

    uint ElGamal::map_to_uint(const Point& message) const {
        const uint& p = m_curve.get_field().modulus();

        if (!p_zero_mask.contains(p)) {
            const size_t l = actual_bit_size(p) >> 1;
            uint zero_mask = (full_bits >> l) << l;
            p_zero_mask.insert({p, zero_mask});
        }

        const uint& zero_mask = p_zero_mask.at(p);
        uint x = message.get_x().value();
        x ^= (x & zero_mask);
        return x;
    }
}   // namespace ECG
