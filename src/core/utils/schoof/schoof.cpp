#include "schoof.h"

#include "division_poly.h"
#include "endomorphism.h"
#include "polynomial.h"
#include "ring.h"
#include "utils/modulo_inversion.h"
#include "utils/primes.h"

namespace elliptic_curve_guide::algorithm::schoof {
    static constexpr size_t c_prime_number_list_size =
        sizeof(primes::prime_number_list) / sizeof(primes::prime_number_list[0]);

    static constexpr size_t c_max_needed_prime_number = 400;

    static std::vector<polynomial::DivisionPoly> get_division_polynomials(const polynomial::Poly& curve_poly,
                                                                          const field::Field& field) {
        using namespace polynomial;
        using namespace field;
        using Element = FieldElement;

        const Element& a = curve_poly[1];
        const Element& b = curve_poly[0];
        const Element a_squared = Element::pow(a, 2);
        const Element b_squared = Element::pow(b, 2);
        auto curve_poly_ptr = std::make_shared<const Poly>(curve_poly);

        DivisionPoly psi_0 = {Poly(field, {0}), curve_poly_ptr, 0};

        DivisionPoly psi_1 = {Poly(field, {1}), curve_poly_ptr, 0};

        DivisionPoly psi_2 = {Poly(field, {2}), curve_poly_ptr, 1};

        Poly psi_3_poly = Poly(
            field,
            {-a_squared, field.element(12) * b, field.element(6) * a, field.element(0), field.element(3)});
        DivisionPoly psi_3 = {std::move(psi_3_poly), curve_poly_ptr, 0};

        Poly psi_4_poly = Poly(field,
                               {-(b_squared << 3) - a_squared * a,
                                -(a * b << 4),
                                -field.element(5) * a_squared,
                                field.element(20) * b,
                                field.element(5) * a,
                                field.element(0),
                                field.element(1)})
                        * field.element(4);
        DivisionPoly psi_4 = {std::move(psi_4_poly), curve_poly_ptr, 1};

        std::vector<DivisionPoly> result = {
            std::move(psi_0), std::move(psi_1), std::move(psi_2), std::move(psi_3), std::move(psi_4)};

        for (size_t i = 5; i < 8; ++i) {
            const size_t n = i >> 1;

            if (i % 2 == 1) {
                DivisionPoly lhs = result[n + 2] * DivisionPoly::pow(result[n], 3);
                lhs.reduce_y();
                DivisionPoly rhs = result[n - 1] * DivisionPoly::pow(result[n + 1], 3);
                rhs.reduce_y();
                DivisionPoly next = lhs - rhs;
                result.emplace_back(std::move(next));
            } else {
                DivisionPoly lhs = result[n + 2] * DivisionPoly::pow(result[n - 1], 2);
                DivisionPoly rhs = result[n - 2] * DivisionPoly::pow(result[n + 1], 2);
                DivisionPoly next = result[n] * (lhs - rhs);
                next.divide_by_y();
                next *= FieldElement::inverse(field.element(2));
                next.reduce_y();
                result.emplace_back(std::move(next));
            }
        }

        return result;
    }

    namespace {
        struct PrimesSet {
            size_t number_of_primes = 0;
            uint accumulated_product = 1;
        };
    }   // namespace

    static uint32_t trace_modulo(const elliptic_curve::EllipticCurve& curve, const uint32_t modulus) {
        using Element = ring::RingElement;
        using endomorphism::End;

        const field::Field& F = curve.get_field();
        const uint& p = F.modulus();
        polynomial::Poly curve_poly(F, {curve.get_b(), curve.get_a(), F.element(0), F.element(1)});

        if (modulus == 2) {
            bool has_2_torsion_points = has_root(curve_poly);
            return has_2_torsion_points ? 0 : 1;
        }

        static std::vector<polynomial::DivisionPoly> division_polynomials =
            get_division_polynomials(curve_poly, F);

        polynomial::Poly h = division_polynomials[modulus].get_x_poly();

        for (;;) {
            ring::Ring R(h);

            Element x = R.element(polynomial::Poly(F, {0, 1}));
            Element one = R.element(polynomial::Poly(F, {1}));

            Element pi_a = Element::pow(x, p);
            auto curve_poly_ptr = std::make_shared<const Element>(R.element(curve_poly));
            Element pi_b = Element::pow(*curve_poly_ptr, (p - 1) >> 1);

            End pi(R, pi_a, pi_b, curve_poly_ptr);
            End pi_squared = pi * pi;
            End id(R, x, one, curve_poly_ptr);
            End::AdditionResult var = id * p;

            if (std::holds_alternative<polynomial::Poly>(var)) {
                h = std::get<polynomial::Poly>(var);
                continue;
            }

            End q = std::get<End>(var);
            var = pi_squared + q;

            if (std::holds_alternative<polynomial::Poly>(var)) {
                h = std::get<polynomial::Poly>(var);
                continue;
            }

            End sum = std::get<End>(var);
            uint32_t c = 0;
            End temp = id;
            bool bad_denominator = false;

            while (temp != sum) {
                ++c;
                assert(c < modulus);
                var = temp + pi;

                if (std::holds_alternative<polynomial::Poly>(var)) {
                    h = std::get<polynomial::Poly>(var);
                    bad_denominator = true;
                    break;
                }

                temp = std::get<End>(var);
            }

            if (bad_denominator) {
                continue;
            }

            return c;
        }
    }

    // Main logic was inspired by https://math.mit.edu/classes/18.783/2015/LectureNotes9.pdf
    uint points_number(const elliptic_curve::EllipticCurve& curve) {
        const uint& p = curve.get_field().modulus();

        uint M = 1;
        uint t = 0;
        size_t pos = 0;
        const uint edge = p << 4;

        while (M * M <= edge) {
            const uint32_t& l = primes::prime_number_list[pos++];
            const uint t_l = trace_modulo(curve, l);

            if (M == 1) {
                t = 1;
            } else {
                uint a = inverse_modulo(M.convert_to<uint32_t>() % l, l);
                a *= M * t_l;
                uint b = inverse_modulo(static_cast<uint>(l), M);
                b *= l * t;
                t = a + b;
            }

            M *= l;

            if (t >= M) {
                t %= M;
            }
        }

        if (t > M >> 1) {
            t = M - t;
        }

        return p + 1 - t;
    }
}   // namespace elliptic_curve_guide::algorithm::schoof
