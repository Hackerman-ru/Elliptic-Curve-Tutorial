#include "schoof.h"

#include "division_poly.h"
#include "modulo_inversion.h"
#include "polynomial.h"
#include "primes.h"
#include "ring.h"

namespace elliptic_curve_guide::algorithm::schoof {
    static constexpr size_t prime_number_list_size =
        sizeof(primes::prime_number_list) / sizeof(primes::prime_number_list[0]);

    static std::vector<polynomial::DivisionPoly> get_division_polynomials(const polynomial::Poly& curve_poly,
                                                                          const field::Field& field) {
        using namespace polynomial;
        using namespace field;
        using Element = FieldElement;

        const Element& a = curve_poly[1];
        const Element& b = curve_poly[0];
        const Element a_squared = Element::pow(a, 2);
        auto curve_poly_ptr = std::make_shared<const Poly>(curve_poly);

        DivisionPoly psi_0 = {Poly(field, {0}), curve_poly_ptr, field.element(1), 0};

        DivisionPoly psi_1 = {Poly(field, {1}), curve_poly_ptr, field.element(1), 0};

        DivisionPoly psi_2 = {Poly(field, {1}), curve_poly_ptr, field.element(2)};

        Poly psi_3_poly = Poly(
            field, {-a_squared, field.element(12) * b, field.element(6) * a_squared, 0, field.element(3)});
        DivisionPoly psi_3 = {std::move(psi_3_poly), curve_poly_ptr, field.element(1), 0};

        Poly psi_4_poly = Poly(field,
                               {(-Element::pow(b, 2) << 3) - a_squared * a,
                                -(a * b << 2),
                                -field.element(5) * a_squared,
                                field.element(20) * b,
                                field.element(5) * a,
                                0,
                                1});
        DivisionPoly psi_4 = {std::move(psi_4_poly), curve_poly_ptr, field.element(4)};

        std::vector<DivisionPoly> result = {
            std::move(psi_0), std::move(psi_1), std::move(psi_2), std::move(psi_3), std::move(psi_4)};

        for (size_t i = 5; i < prime_number_list_size; ++i) {
            const size_t n = i >> 1;

            if (i % 2 == 1) {
                DivisionPoly next = result[n + 2] * DivisionPoly::pow(result[n], 3)
                                  - result[n - 1] * DivisionPoly::pow(result[n + 1], 3);
                next.reduce_y();
                result.emplace_back(std::move(next));
            } else {
                DivisionPoly next = result[n]
                                  * (result[n + 2] * DivisionPoly::pow(result[n - 1], 2)
                                     - result[n - 2] * DivisionPoly::pow(result[n + 1], 2));
                next.divide_by_2_y();
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
        const field::Field& F = curve.get_field();
        const uint& p = F.modulus();
        polynomial::Poly curve_poly(F, {curve.get_b(), curve.get_a(), F.element(0), F.element(1)});

        if (modulus == 2) {
            bool has_2_torsion_points = has_root(curve_poly);
            return has_2_torsion_points ? 0 : 1;
        }

        static std::vector<polynomial::DivisionPoly> division_polynomials =
            get_division_polynomials(curve_poly, F);
    }

    // Main logic was inspired by https://math.mit.edu/classes/18.783/2015/LectureNotes9.pdf
    uint points_number(const elliptic_curve::EllipticCurve& curve) {
        const uint& p = curve.get_field().modulus();

        uint M = 1;
        uint t = 0;
        size_t pos = 0;
        const uint edge = p << 4;

        while (M * M <= edge) {
            const uint32_t& l = primes::prime_number_list[pos];
            const uint t_l = trace_modulo(curve, l);

            uint a = inverse_modulo(M.convert_to<uint32_t>() % l, l);
            a *= M * t_l;
            uint b = inverse_modulo(static_cast<uint>(l), M);
            b *= l * t;

            t = a + b;
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
