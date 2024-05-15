#include "ecdsa.h"

#include "utils/random.h"
#include "utils/schoof.h"
#include "utils/uint-algorithms.h"

namespace elliptic_curve_guide::algorithm::encryption {
    static elliptic_curve::EllipticCurve generate_random_elliptic_curve(const field::Field& field) {}

    static uint get_largest_prime_divisor(uint value) {
        uint result = 1;

        for (uint divisor = 2; divisor * divisor <= value; ++divisor) {
            while (value % divisor == 0) {
                value /= divisor;
                result = divisor;
            }
        }

        return result;
    }

    ECDSA ECDSA::generate(const uint& field_order, const uint& security_level) {
        const field::Field F(field_order);

        for (;;) {
            const elliptic_curve::EllipticCurve E = generate_random_elliptic_curve(F);
            const uint N = algorithm::schoof::points_number(E);
            const uint n = get_largest_prime_divisor(N);

            if (n == field_order) {
                continue;
            }

            if (actual_bit_size(n) <= security_level) {
                continue;
            }

            uint p_k = field_order;
            bool divide = false;

            for (int k = 1; k < 20; ++k) {
                if ((p_k - 1) % n != 0) {
                    divide = true;
                    break;
                }

                p_k *= field_order;
            }

            if (divide) {
                continue;
            }

            const uint h = N / n;
            elliptic_curve::EllipticCurvePoint P = h * E.random_point();

            while (P.is_zero()) {
                P = h * E.random_point();
            }

            return ECDSA(F, E, P, n, h);
        }
    }

    ECDSA::ECDSA(const Field& field, const Curve& elliptic_curve, const Point& generator, const uint& n,
                 const uint& h) :
        m_field(field), m_elliptic_curve(elliptic_curve), m_generator(generator), m_n(n), m_h(h) {};

    ECDSA::Keys ECDSA::generate_keys() const {
        uint d = random::generate_random_non_zero_uint_modulo(m_n);
        Point Q = d * m_generator;
        return {.public_key = Q, .private_key = d};
    }

    ECDSA::Signature ECDSA::generate_signature(const uint& private_key, const uint& message) const {
        const Field F(m_n);

        for (;;) {
            const uint k = random::generate_random_non_zero_uint_modulo(m_n);

            const Point P = k * m_generator;
            const uint& r = P.get_x().value();

            if (r == 0) {
                continue;
            }

            const Element edr = F.element(message) + F.element(private_key) * F.element(r);
            const uint s = (Element::inverse(F.element(k)) * edr).value();

            if (s == 0) {
                continue;
            }

            return {.r = r, .s = s};
        }

        return {};
    }

    bool ECDSA::is_correct_signature(const Point& public_key, const uint& message,
                                     const Signature& signature) const {
        const uint& r = signature.r;
        const uint& s = signature.s;

        if (r == 0 || s == 0) {
            return false;
        }

        if (r >= m_n || s >= m_n) {
            return false;
        }

        const Field F = Field(m_n);
        const Element w = Element::inverse(F.element(s));
        const Element u1 = F.element(message) * w;
        const Element u2 = F.element(r) * w;
        const Point X = u1.value() * m_generator + u2.value() * public_key;

        if (X.is_zero()) {
            return false;
        }

        const uint& v = X.get_x().value();
        return v == r;
    }
}   // namespace elliptic_curve_guide::algorithm::encryption
