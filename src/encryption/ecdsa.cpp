#include "ecdsa.h"

#include "utils/random.h"

namespace elliptic_curve_guide::algorithm::encryption {
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
            const Element k = random::generate_random_non_zero_field_element(F);

            const Point P = k * m_generator;
            const uint& r = P.get_x().value();

            if (r == 0) {
                continue;
            }

            const Element edr = F.element(message) + F.element(private_key) * F.element(r);
            const uint& s = (Element::inverse(k) * edr).value();

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
        const Point X = u1 * m_generator + u2 * public_key;

        if (X.is_zero()) {
            return false;
        }

        const uint& v = X.get_x().value();
        return v == r;
    }
}   // namespace elliptic_curve_guide::algorithm::encryption
