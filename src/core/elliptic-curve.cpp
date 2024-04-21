#include "elliptic-curve.h"

#include "map"

namespace ECG {
    EllipticCurve::EllipticCurve(const FieldElement& a, const FieldElement& b, Field F) :
        m_a {std::make_shared<const FieldElement>(a)},
        m_b {std::make_shared<const FieldElement>(b)},
        m_F {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(FieldElement&& a, const FieldElement& b, Field F) :
        m_a {std::make_shared<const FieldElement>(std::move(a))},
        m_b {std::make_shared<const FieldElement>(b)},
        m_F {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(const FieldElement& a, FieldElement&& b, Field F) :
        m_a {std::make_shared<const FieldElement>(a)},
        m_b {std::make_shared<const FieldElement>(std::move(b))},
        m_F {std::make_shared<const Field>(std::move(F))} {}

    EllipticCurve::EllipticCurve(FieldElement&& a, FieldElement&& b, Field F) :
        m_a {std::make_shared<const FieldElement>(std::move(a))},
        m_b {std::make_shared<const FieldElement>(std::move(b))},
        m_F {std::make_shared<const Field>(std::move(F))} {}

    uint EllipticCurve::points_number() const {   // TODO
        return uint();
    }

    static std::map<uint, FieldElement> b_cache;

    static FieldElement find_b(const FieldElement& one, const uint& degree) {
        if (b_cache.contains(degree)) {
            return b_cache.at(degree);
        }

        FieldElement b = one + one;

        while (b.pow(degree) == one) {
            b += one;
        }

        b_cache.insert({degree, b});
        return b;
    }

    namespace {
        struct Decomposition {
            size_t power_of_two;
            uint residue;
        };
    }   // namespace

    static Decomposition decompose(const uint& value) {
        size_t power_of_two = 0;
        uint residue = value;

        while (residue != 0 && (residue & 0b1) == 0) {
            ++power_of_two;
            residue >>= 1;
        }

        return {power_of_two, residue};
    }

    static size_t order_of_two(uint value) {
        size_t power_of_two = 0;

        while (value != 0 && (value & 0b1) == 0) {
            ++power_of_two;
        }

        return power_of_two;
    }

    std::optional<FieldElement> EllipticCurve::find_y(const FieldElement& x) const {
        FieldElement value = x.pow(3) + *m_a * x + *m_b;

        if (!value.is_invertible()) {
            return std::nullopt;
        }

        const uint& p = m_F->modulus();
        const FieldElement one = m_F->element(1);

        if (!m_cache) {
            uint degree = (p - 1) >> 1;
            Decomposition decomposition = decompose(p - 1);
            FieldElement b = find_b(one, degree);
            FieldElement inverse_two = one / m_F->element(2);
            size_t e = decomposition.power_of_two;

            std::vector<FieldElement> second_powers = {b << 1};
            second_powers.reserve(e - 1);

            for (size_t i = 1; i < e - 1; ++i) {
                second_powers.emplace_back(second_powers[i - 1] << 1);
            }

            std::vector<FieldElement> second_u_powers = {b.pow(decomposition.residue)};
            second_u_powers.reserve(e);

            for (size_t i = 1; i < e; ++i) {
                second_u_powers.emplace_back(second_u_powers[i - 1] << 1);
            }

            Cache cache = {.degree = std::move(degree),
                           .power_of_two = std::move(decomposition.power_of_two),
                           .residue = std::move(decomposition.residue),
                           .b = std::move(b),
                           .inverse_two = std::move(inverse_two),
                           .b_second_powers = std::move(second_powers),
                           .b_second_u_powers = std::move(second_u_powers)};

            m_cache.emplace(std::move(cache));
        }

        const Cache& cache = m_cache.value();
        FieldElement current_z_u_2 = value.pow(cache.degree);

        if (current_z_u_2 != one) {
            return std::nullopt;
        }

        if ((p & 0b11) == 3) {
            return value.pow((cache.degree + 1) >> 1);
        }

        const FieldElement& inverse_two = cache.inverse_two;
        const size_t& e = cache.power_of_two;

        std::vector<size_t> two_orders;
        size_t current_r = e;

        while (current_r > 0 && current_z_u_2 == one) {
            current_z_u_2 *= inverse_two;
            --current_r;
        }

        if (current_z_u_2 != one) {
            ++current_r;
        }

        two_orders.emplace_back(current_r);
        FieldElement current_z = value;

        while (current_r != 0) {
            current_z *= cache.b_second_powers[e - current_r];
            current_z_u_2 *= cache.b_second_u_powers[e - 1];
            current_r -= 1;

            assert(current_z_u_2 == one && "EllipticCurve::find_y : wrong implementation");

            while (current_r > 0 && current_z_u_2 == one) {
                current_z_u_2 *= inverse_two;
                --current_r;
            }

            if (current_z_u_2 != one) {
                ++current_r;
            }

            two_orders.emplace_back(current_r);
        }

        FieldElement current_x = current_z.pow((cache.residue + 1) >> 1);
        const size_t n = two_orders.size();

        for (size_t i = 0; i + 1 < n; ++i) {
            current_x /= cache.b_second_powers[e - two_orders[n - i - 2] - 1];
        }

        assert((current_x << 1) == value && "EllipticCurve::find_y : wrong implementation");
        return current_x;
    }
}   // namespace ECG
