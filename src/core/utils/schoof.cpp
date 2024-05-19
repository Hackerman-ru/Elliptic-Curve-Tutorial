#include "schoof.h"

#include "modulo_inversion.h"
#include "polynomial.h"
#include "primes.h"
#include "ring.h"

namespace elliptic_curve_guide::algorithm::schoof {
    namespace {
        struct PrimesSet {
            size_t number_of_primes = 0;
            uint accumulated_product = 1;
        };
    }   // namespace

    static constexpr size_t prime_number_list_size =
        sizeof(primes::prime_number_list) / sizeof(primes::prime_number_list[0]);

    static uint32_t trace_modulo(const elliptic_curve::EllipticCurve& curve, const uint32_t modulus) {
        using field::Field;
        using polynomial::Poly;
        using ring::Ring;
        using ring::RingElement;

        const field::Field& F = curve.get_field();
        const uint& p = F.modulus();

        if (modulus == 2) {
            Poly f(F, {curve.get_b(), curve.get_a(), F.element(0), F.element(1)});
            bool has_2_torsion_points = has_root(f);
            return has_2_torsion_points ? 0 : 1;
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
