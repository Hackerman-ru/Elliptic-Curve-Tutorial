#include "schoof.h"

#include "polynomial.h"
#include "primes.h"
#include "ring.h"

namespace elliptic_curve_guide::algorithm::schoof {
    namespace {
        struct S {
            size_t L = 0;
            uint N = 1;
        };
    }   // namespace

    static constexpr size_t prime_number_list_size =
        sizeof(primes::prime_number_list) / sizeof(primes::prime_number_list[0]);

    static constexpr S count_primes_number(const uint& p) {
        uint edge_value = p << 4;
        S result;

        do {
            const uint l = primes::prime_number_list[result.L];

            if (l == p) {
                continue;
            }

            result.N *= l;
            ++result.L;

            if (result.L >= prime_number_list_size) {
                break;
            }
        } while (result.N <= edge_value);

        return result;
    }

    uint points_number(const elliptic_curve::EllipticCurve& curve) {
        using field::Field;
        using polynomial::Poly;
        using ring::Ring;
        using ring::RingElement;

        const field::Field& F = curve.get_field();
        const uint& p = F.modulus();
        S primes_set = count_primes_number(p);
        const size_t& L = primes_set.L;
        const uint& N = primes_set.N;

        std::vector<uint32_t> modulo_list;
        modulo_list.reserve(L);

        for (size_t i = 0; i < L; ++i) {
            const uint32_t& l = primes::prime_number_list[i];

            if (l == 2) {
                Poly f(F, {curve.get_b(), curve.get_a(), F.element(0), F.element(1)});
                bool has_2_torsion_points = has_root(f);
                modulo_list.emplace_back(has_2_torsion_points ? 0 : 1);
                continue;
            }
        }

        return 0;   // TODO
    }
}   // namespace elliptic_curve_guide::algorithm::schoof
