#include "ecdsa.h"

#include "cthash/sha3/sha3-512.hpp"

namespace ECG {
    namespace {
        struct Coefficients {
            FieldElement a;
            FieldElement b;
        };

        struct EllipticCurveParameters {
            FieldElement a;
            FieldElement b;
        };
    }   // namespace

    static constexpr uint hash_length = 512;

    static uint bit_size(uint value) {
        uint result = 0;

        while (value != 0) {
            ++result;
            value >>= 1;
        }

        return result;
    }

    static uint hash_value(const uint& value) {
        const auto h = cthash::simple<cthash::sha3_512>(value.convert_to<std::string>());
        const size_t shift = sizeof(h[0]) * 8;
        uint result = 0;

        for (const auto& byte : h) {
            result <<= shift;
            result += byte;
        }

        return result;
    }

    static EllipticCurveParameters generate_random_elliptic_curve(const uint& p) {
        uint t = bit_size(p);
        uint s = (t - 1) / hash_length;
        size_t v = (t - (s << 9)).convert_to<size_t>();
        uint S = generate_random_uint();
        uint h = hash_value(S);
        uint r0 = (h << (512 - v)) >> v;
        uint R0 = (r0 << 1) >> 1;

        for (uint i = 1; i < s; ++i) {
        }
    }

    static Coefficients generate_random_a_b(const Field& F) {
        FieldElement a = F.element(0);
        FieldElement b = F.element(0);
        FieldElement zero = F.element(0);

        do {
            a = F.element(generate_random_uint());
            b = F.element(generate_random_uint());
            zero = (a.pow(3) << 2) + F.element(27) * b.pow(2);
        } while (!zero.is_invertible());

        return {a, b};
    }

    static uint get_greatest_divisor(uint value) {
        uint result = 1;

        for (uint divisor = 2; divisor * divisor <= value; ++divisor) {
            while (value % divisor == 0) {
                value /= divisor;
                result = divisor;
            }
        }

        return result;
    }

    static bool satisfies(const uint& n, const uint& field_order, const uint& security_level) {
        if (n == field_order) {
            return false;
        }

        if (bit_size(n) < security_level) {
            return false;
        }

        uint temp = field_order;

        for (int k = 1; k < 20; ++k) {
            if ((temp - 1) % n != 0) {
                return false;
            }

            temp *= field_order;
        }

        return true;
    }

    std::optional<ECDSA> ECDSA::generate(const uint& field_order, const uint& security_level) {
        uint bitsize = bit_size(field_order);

        if (security_level < 160 || security_level > bitsize) {
            return std::nullopt;
        }

        if (security_level < (2 + (bitsize >> 1))) {
            return std::nullopt;
        }

        Field F(field_order);
        Coefficients curve_coefficients = generate_random_a_b(F);
        EllipticCurve E(curve_coefficients.a, curve_coefficients.b, F);
        uint N = E.points_number();
        uint n = get_greatest_divisor(N);

        while (!satisfies(n, field_order, security_level)) {
            curve_coefficients = generate_random_a_b(F);
            E = EllipticCurve(curve_coefficients.a, curve_coefficients.b, F);
            N = E.points_number();
            n = get_greatest_divisor(N);
        }

        uint h = N / n;
        EllipticCurvePoint<> G = h * E.random_point<>();

        while (!G.is_zero()) {
            G = h * E.random_point<>();
        }

        return ECDSA(F, E, G, n, h);
    }

    Keys ECDSA::generate_keys() const {
        FieldElement d = m_parameters.m_Field.element(generate_random_uint());

        while (!d.is_invertible()) {
            d = m_parameters.m_Field.element(generate_random_uint());
        }

        EllipticCurvePoint<> Q = d.value() * m_parameters.m_generator;
        return {.public_key = Q, .private_key = d};
    }

    Signature ECDSA::generate_signature(const FieldElement& private_key, const uint& message) const {
        const Field& F = m_parameters.m_Field;
        FieldElement k = F.element(generate_random_uint());

        while (!k.is_invertible()) {
            k = F.element(generate_random_uint());
        }

        EllipticCurvePoint P = k.value() * m_parameters.m_generator;
        FieldElement r = P.get_x();

        while (!r.is_invertible()) {
            k = F.element(generate_random_uint());

            while (!k.is_invertible()) {
                k = F.element(generate_random_uint());
            }

            P = k.value() * m_parameters.m_generator;
            r = P.get_x();
        }

        FieldElement s = FieldElement::inverse(k) * (F.element(message) + private_key * r);
        while (!s.is_invertible()) {
            k = F.element(generate_random_uint());

            while (!k.is_invertible()) {
                k = F.element(generate_random_uint());
            }

            P = k.value() * m_parameters.m_generator;
            r = P.get_x();
            s = FieldElement::inverse(k) * (F.element(message) + private_key * r);
        }

        return {.r = r, .s = s};
    }

    bool ECDSA::is_correct_signature(const EllipticCurvePoint<>& public_key, const uint& message,
                                     const Signature& signature) const {
        const uint& n = m_parameters.m_n;
        const FieldElement& r = signature.r;
        const FieldElement& s = signature.s;

        if (r.modulus() != n || s.modulus() != n) {
            return false;
        }

        if (!r.is_invertible() || !s.is_invertible()) {
            return false;
        }

        const Field& F = m_parameters.m_Field;
        FieldElement w = FieldElement::inverse(s);
        FieldElement e = F.element(message);
        FieldElement u1 = e * w;
        FieldElement u2 = r * w;
        EllipticCurvePoint<> X = u1.value() * m_parameters.m_generator + u2.value() * public_key;

        if (X.is_zero()) {
            return false;
        }

        FieldElement v = X.get_x();
        return v == r;
    }

}   // namespace ECG
