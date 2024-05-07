#include "ecdsa.h"

#include "utils/random.h"

namespace ECG {
    namespace {
        struct Coefficients {
            FieldElement a;
            FieldElement b;
        };
    }   // namespace

    static constexpr uint c_hash_length = 512;
    static constexpr int c_attempts_number = 10000;

    //static uint bit_size(uint value) {
    //    uint result = 0;

    //    while (value != 0) {
    //        ++result;
    //        value >>= 1;
    //    }

    //    return result;
    //}

    //static uint hash_value(const uint& value) {
    //    const auto h = cthash::simple<cthash::sha3_512>(value.convert_to<std::string>());
    //    const size_t shift = sizeof(h[0]) * 8;
    //    uint result = 0;

    //    for (const auto& byte : h) {
    //        result <<= shift;
    //        result += byte;
    //    }

    //    return result;
    //}

    //static EllipticCurveParameters generate_random_elliptic_curve(const uint& p) {
    //    uint t = bit_size(p);
    //    uint s = (t - 1) / c_hash_length;
    //    size_t v = (t - (s << 9)).convert_to<size_t>();
    //    uint S = generate_random_uint();
    //    uint h = hash_value(S);
    //    uint r0 = (h << (512 - v)) >> v;
    //    uint R0 = (r0 << 1) >> 1;

    //    for (uint i = 1; i < s; ++i) {
    //    }
    //}

    //static Coefficients generate_random_a_b(const Field& F) {
    //    FieldElement a = F.element(0);
    //    FieldElement b = F.element(0);
    //    FieldElement non_zero = F.element(0);
    //    const uint& p = F.modulus();

    //    do {
    //        a = F.element(generate_random_uint_modulo(p));
    //        b = F.element(generate_random_uint_modulo(p));
    //        non_zero = (FieldElement::pow(a, 3) << 2) + F.element(27) * FieldElement::pow(b, 2);
    //    } while (!non_zero.is_invertible());

    //    return {a, b};
    //}

    //static uint get_greatest_prime_divisor(uint value) {
    //    uint result = 1;

    //    for (uint divisor = 2; divisor * divisor <= value; ++divisor) {
    //        while (value % divisor == 0) {
    //            value /= divisor;
    //            result = divisor;
    //        }
    //    }

    //    return result;
    //}

    //static bool satisfies(const uint& n, const uint& field_order, const uint& security_level) {
    //    if (n == field_order) {
    //        return false;
    //    }

    //    if (bit_size(n) < security_level) {
    //        return false;
    //    }

    //    uint temp = field_order;

    //    for (int k = 1; k < 20; ++k) {
    //        if ((temp - 1) % n != 0) {
    //            return false;
    //        }

    //        temp *= field_order;
    //    }

    //    return true;
    //}

    //ECDSA ECDSA::generate(const uint& field_order, const uint& security_level) {
    //    Field F(field_order);
    //    Coefficients curve_coefficients = generate_random_a_b(F);
    //    EllipticCurve E(curve_coefficients.a, curve_coefficients.b, F);
    //    uint N = E.points_number();
    //    uint n = get_greatest_prime_divisor(N);

    //    while (!satisfies(n, field_order, security_level)) {
    //        curve_coefficients = generate_random_a_b(F);
    //        E = EllipticCurve(curve_coefficients.a, curve_coefficients.b, F);
    //        N = E.points_number();
    //        n = get_greatest_prime_divisor(N);
    //    }

    //    uint h = N / n;
    //    EllipticCurvePoint<> G = h * E.random_point<>();

    //    while (!G.is_zero()) {
    //        G = h * E.random_point<>();
    //    }

    //    return ECDSA(F, E, G, n, h);
    //}

    Keys ECDSA::generate_keys() const {
        uint d = generate_random_non_zero_uint_modulo(m_parameters.m_n);
        EllipticCurvePoint<> Q = d * m_parameters.m_generator;
        return {.public_key = Q, .private_key = d};
    }

    Signature ECDSA::generate_signature(const uint& private_key, const uint& message) const {
        const Field F = Field(m_parameters.m_n);

        for (int i = 0; i < c_attempts_number; ++i) {
            const uint k = generate_random_non_zero_uint_modulo(m_parameters.m_n);

            const EllipticCurvePoint P = k * m_parameters.m_generator;
            const uint& r = P.get_x().value();

            if (r == 0) {
                continue;
            }

            const FieldElement edr = F.element(message) + F.element(private_key) * F.element(r);
            const uint s = (FieldElement::inverse(F.element(k)) * edr).value();

            if (s == 0) {
                continue;
            }

            return {.r = r, .s = s};
        }

        return {};
    }

    bool ECDSA::is_correct_signature(const EllipticCurvePoint<>& public_key, const uint& message,
                                     const Signature& signature) const {
        const uint& n = m_parameters.m_n;
        const uint& r = signature.r;
        const uint& s = signature.s;

        if (r == 0 || s == 0) {
            return false;
        }

        if (r >= n || s >= n) {
            return false;
        }

        const Field F = Field(m_parameters.m_n);
        const FieldElement w = FieldElement::inverse(F.element(s));
        const FieldElement u1 = F.element(message) * w;
        const FieldElement u2 = F.element(r) * w;
        const EllipticCurvePoint<> X = u1.value() * m_parameters.m_generator + u2.value() * public_key;

        if (X.is_zero()) {
            return false;
        }

        const uint& v = X.get_x().value();
        return v == r;
    }
}   // namespace ECG
