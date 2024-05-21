// clang-format off
#include "pch.h"
// clang-format on

#include "elliptic-curve.h"
#include "field.h"
#include "utils/primes.h"
#include "utils/random.h"

using namespace elliptic_curve_guide;
using namespace field;
using namespace elliptic_curve;
using namespace algorithm::random;
using enum CoordinatesType;

static constexpr size_t c_primes_n = 150;
static constexpr uint c_good_p = 53617;   // e = 4
static constexpr uint c_big_p = "0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff";

static constexpr size_t c_correctness_test_find_y_n = 20;
static constexpr size_t c_correctness_test_random_n = 20;
static constexpr size_t c_correctness_test_zero_n = 20;
static constexpr size_t c_correctness_test_addition_n = 20;
static constexpr size_t c_correctness_test_kp_n = 20;
static constexpr uint c_correctness_max_k_n = 100;

static constexpr size_t c_stress_test_find_y_n = 50;
static constexpr size_t c_stress_test_kp_n = 25;
static constexpr uint c_stress_max_k_n = "0xffffff12fffabcffff34fffffffffffff123ff";

struct Coordinates {
    std::string x;
    std::string y;
};

bool operator==(const Coordinates& lhs, const Coordinates& rhs) {
    return lhs.x == rhs.x && lhs.y == rhs.y;
}

template<CoordinatesType type>
Coordinates get_coordinates(const EllipticCurvePoint<type>& point) {
    return Coordinates {.x = point.get_x().value().convert_to<std::string>(),
                        .y = point.get_y().value().convert_to<std::string>()};
}

template<CoordinatesType type>
testing::AssertionResult has_point(const std::optional<EllipticCurvePoint<type>>& opt) {
    if (opt.has_value()) {
        return testing::AssertionSuccess();
    } else {
        return testing::AssertionFailure() << opt.has_value() << "curve has point at this x";
    }
}

#define UINT_EQ(lhs, rhs)                                    \
    if (lhs != rhs) {                                        \
        std::string lhs_str = lhs.convert_to<std::string>(); \
        std::string rhs_str = rhs.convert_to<std::string>(); \
        ASSERT_EQ(lhs_str, rhs_str);                         \
    }

#define FIELD_EQ(lhs, rhs)                                           \
    if (lhs != rhs) {                                                \
        std::string lhs_str = lhs.value().convert_to<std::string>(); \
        std::string rhs_str = rhs.value().convert_to<std::string>(); \
        ASSERT_EQ(lhs_str, rhs_str);                                 \
    }

#define FIELD_NEQ(lhs, rhs)                                          \
    if (lhs == rhs) {                                                \
        std::string lhs_str = lhs.value().convert_to<std::string>(); \
        std::string rhs_str = rhs.value().convert_to<std::string>(); \
        ASSERT_NE(lhs_str, rhs_str);                                 \
    }

#define POINT_EQ(lhs, rhs)                            \
    if (lhs != rhs) {                                 \
        Coordinates lhs_coord = get_coordinates(lhs); \
        Coordinates rhs_coord = get_coordinates(rhs); \
        ASSERT_EQ(lhs_coord, rhs_coord);              \
    }

// Common tests
// Simple tests
TEST(SimpleTest, Creation) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    auto opt1 = E.point_with_x_equal_to(F.element(4));
    auto opt2 = E.point(F.element(4), F.element(18));
    auto P = E.random_point();
}

TEST(SimpleTest, FindY) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    auto opt = E.point_with_x_equal_to(F.element(4));
    ASSERT_TRUE(has_point(opt));
    const auto& point = opt.value();
    const uint correct_y = 18;
    UINT_EQ(point.get_y().value(), correct_y);
}

// Correctness tests
TEST(CorrectnessTest, FindY) {
    for (size_t i = 0; i < c_correctness_test_find_y_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        uint p = primes::prime_number_list[pos];
        Field F(p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        for (uint j = 0; j < p; ++j) {
            FieldElement x = F.element(j);
            auto opt = E.point_with_x_equal_to(x);

            if (!x.is_invertible()) {
                ASSERT_TRUE(opt.has_value());
                const auto& P = opt.value();
                ASSERT_TRUE(P.is_zero());
                continue;
            }

            FieldElement value = FieldElement::pow(x, 3) + a * x + b;

            if (opt.has_value()) {
                const auto& P = opt.value();
                FieldElement y = P.get_y();
                FieldElement y2 = FieldElement::pow(y, 2);
                FIELD_EQ(y2, value);
            } else {
                value.pow((p - 1) >> 1);
                FIELD_NEQ(value, F.element(1));
            }
        }
    }

    {
        Field F(c_good_p);
        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        for (uint j = 1; j < c_good_p; j <<= 1) {
            FieldElement x = F.element(j);
            auto opt = E.point_with_x_equal_to(x);

            if (!x.is_invertible()) {
                ASSERT_TRUE(opt.has_value());
                const auto& P = opt.value();
                ASSERT_TRUE(P.is_zero());
                continue;
            }

            FieldElement value = FieldElement::pow(x, 3) + a * x + b;

            if (opt.has_value()) {
                const auto& P = opt.value();
                FieldElement y = P.get_y();
                FieldElement y2 = FieldElement::pow(y, 2);
                FIELD_EQ(y2, value);
            } else {
                value.pow((c_good_p - 1) >> 1);
                FIELD_NEQ(value, F.element(1));
            }
        }
    }

    {
        Field F(c_big_p);
        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        for (uint j = 1; j < c_big_p; j <<= 32) {
            FieldElement x = F.element(j);
            auto opt = E.point_with_x_equal_to(x);

            if (!x.is_invertible()) {
                ASSERT_TRUE(opt.has_value());
                const auto& P = opt.value();
                ASSERT_TRUE(P.is_zero());
                continue;
            }

            FieldElement value = FieldElement::pow(x, 3) + a * x + b;

            if (opt.has_value()) {
                const auto& P = opt.value();
                FieldElement y = P.get_y();
                FieldElement y2 = FieldElement::pow(y, 2);
                FIELD_EQ(y2, value);
            } else {
                value.pow((c_big_p - 1) >> 1);
                FIELD_NEQ(value, F.element(1));
            }
        }
    }
}

// Stress tests
TEST(StressTest, FindY) {
    for (size_t i = 0; i < c_stress_test_find_y_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        uint p = primes::prime_number_list[pos];
        Field F(p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        for (uint j = 0; j < p; ++j) {
            FieldElement x = F.element(j);
            auto opt = E.point_with_x_equal_to(x);
        }
    }

    {
        Field F(c_good_p);
        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        for (uint j = 1; j < c_good_p; j <<= 1) {
            FieldElement x = F.element(j);
            auto opt = E.point_with_x_equal_to(x);
        }
    }

    {
        Field F(c_big_p);
        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        for (uint j = 1; j < c_big_p; j <<= 32) {
            FieldElement x = F.element(j);
            auto opt = E.point_with_x_equal_to(x);
        }
    }
}

// Normal Coordinates tests
// Correctness tests
TEST(CorrectnessTest, RandomNormal) {
    for (size_t i = 0; i < c_correctness_test_random_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Normal> P = E.random_point<Normal>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }

    {
        Field F(c_big_p);
        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Normal> P = E.random_point<Normal>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }
}

TEST(CorrectnessTest, ZeroNormal) {
    for (size_t i = 0; i < c_correctness_test_zero_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Normal> Z = E.null_point<Normal>();
        ASSERT_TRUE(Z.is_zero());
        FIELD_EQ(Z.get_x(), F.element(0));
        FIELD_EQ(Z.get_y(), F.element(1));
    }
}

TEST(CorrectnessTest, AdditionNormal) {
    for (size_t i = 0; i < c_correctness_test_addition_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<Normal> Z = E.null_point<Normal>();
        ASSERT_TRUE(Z.is_zero());

        EllipticCurvePoint<Normal> two_Z = Z + Z;
        ASSERT_TRUE(two_Z.is_zero());
        POINT_EQ(two_Z, Z);

        EllipticCurvePoint<Normal> P = E.random_point<Normal>();
        EllipticCurvePoint<Normal> Z_P = Z + P;
        POINT_EQ(Z_P, P);

        EllipticCurvePoint<Normal> Q = E.random_point<Normal>();
        EllipticCurvePoint<Normal> P_Q = P + Q;

        if (P.is_zero()) {
            POINT_EQ(P_Q, Q);
            continue;
        } else if (Q.is_zero()) {
            POINT_EQ(P_Q, P);
            continue;
        } else if (P == -Q) {
            ASSERT_TRUE(P_Q.is_zero());
            continue;
        } else if (P == Q) {
            const FieldElement& x = P.get_x();
            const FieldElement& y = P.get_y();

            if (!y.is_invertible()) {
                ASSERT_TRUE(P_Q.is_zero());
                continue;
            }

            FieldElement k = (F.element(3) * FieldElement::pow(x, 2) + a) / (y << 1);
            FieldElement correct_x = FieldElement::pow(k, 2) - (x << 1);
            FieldElement correct_y = k * (x - correct_x) - y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        } else {
            const FieldElement& lhs_x = P.get_x();
            const FieldElement& lhs_y = P.get_y();
            const FieldElement& rhs_x = Q.get_x();
            const FieldElement& rhs_y = Q.get_y();

            FieldElement k = (rhs_y - lhs_y) / (rhs_x - lhs_x);
            FieldElement correct_x = FieldElement::pow(k, 2) - lhs_x - rhs_x;
            FieldElement correct_y = k * (lhs_x - correct_x) - lhs_y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        }
    }
}

TEST(CorrectnessTest, kPNormal) {
    for (size_t i = 0; i < c_correctness_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<Normal> P = E.random_point<Normal>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<Normal> kP = k * P;

        EllipticCurvePoint<Normal> correct_kP = E.null_point<Normal>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<Normal> P = E.random_point<Normal>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<Normal> kP = k * P;

        EllipticCurvePoint<Normal> correct_kP = E.null_point<Normal>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }
}

// Stress tests
TEST(StressTest, kPNormal) {
    for (size_t i = 0; i < c_stress_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        uint p = primes::prime_number_list[pos];
        Field F(p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Normal> P = E.random_point<Normal>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Normal> P = E.random_point<Normal>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }
}

// Jacobi Coordinates tests
// Correctness tests
TEST(CorrectnessTest, RandomJacobi) {
    for (size_t i = 0; i < c_correctness_test_random_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Jacobi> P = E.random_point<Jacobi>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }

    {
        Field F(c_big_p);
        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Jacobi> P = E.random_point<Jacobi>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }
}

TEST(CorrectnessTest, ZeroJacobi) {
    for (size_t i = 0; i < c_correctness_test_zero_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Jacobi> Z = E.null_point<Jacobi>();
        ASSERT_TRUE(Z.is_zero());
        FIELD_EQ(Z.get_x(), F.element(0));
        FIELD_EQ(Z.get_y(), F.element(1));
    }
}

TEST(CorrectnessTest, AdditionJacobi) {
    for (size_t i = 0; i < c_correctness_test_addition_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<Jacobi> Z = E.null_point<Jacobi>();
        ASSERT_TRUE(Z.is_zero());

        EllipticCurvePoint<Jacobi> two_Z = Z + Z;
        ASSERT_TRUE(two_Z.is_zero());
        POINT_EQ(two_Z, Z);

        EllipticCurvePoint<Jacobi> P = E.random_point<Jacobi>();
        EllipticCurvePoint<Jacobi> Z_P = Z + P;
        POINT_EQ(Z_P, P);

        EllipticCurvePoint<Jacobi> Q = E.random_point<Jacobi>();
        EllipticCurvePoint<Jacobi> P_Q = P + Q;

        if (P.is_zero()) {
            POINT_EQ(P_Q, Q);
            continue;
        } else if (Q.is_zero()) {
            POINT_EQ(P_Q, P);
            continue;
        } else if (P == -Q) {
            ASSERT_TRUE(P_Q.is_zero());
            continue;
        } else if (P == Q) {
            const FieldElement& x = P.get_x();
            const FieldElement& y = P.get_y();

            if (!y.is_invertible()) {
                ASSERT_TRUE(P_Q.is_zero());
                continue;
            }

            FieldElement k = (F.element(3) * FieldElement::pow(x, 2) + a) / (y << 1);
            FieldElement correct_x = FieldElement::pow(k, 2) - (x << 1);
            FieldElement correct_y = k * (x - correct_x) - y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        } else {
            const FieldElement& lhs_x = P.get_x();
            const FieldElement& lhs_y = P.get_y();
            const FieldElement& rhs_x = Q.get_x();
            const FieldElement& rhs_y = Q.get_y();

            FieldElement k = (rhs_y - lhs_y) / (rhs_x - lhs_x);
            FieldElement correct_x = FieldElement::pow(k, 2) - lhs_x - rhs_x;
            FieldElement correct_y = k * (lhs_x - correct_x) - lhs_y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        }
    }
}

TEST(CorrectnessTest, kPJacobi) {
    for (size_t i = 0; i < c_correctness_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<Jacobi> P = E.random_point<Jacobi>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<Jacobi> kP = k * P;

        EllipticCurvePoint<Jacobi> correct_kP = E.null_point<Jacobi>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<Jacobi> P = E.random_point<Jacobi>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<Jacobi> kP = k * P;

        EllipticCurvePoint<Jacobi> correct_kP = E.null_point<Jacobi>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }
}

// Stress tests
TEST(StressTest, kPJacobi) {
    for (size_t i = 0; i < c_stress_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        uint p = primes::prime_number_list[pos];
        Field F(p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Jacobi> P = E.random_point<Jacobi>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Jacobi> P = E.random_point<Jacobi>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }
}

// Projective Coordinates tests
// Correctness tests
TEST(CorrectnessTest, RandomProjective) {
    for (size_t i = 0; i < c_correctness_test_random_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Projective> P = E.random_point<Projective>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }

    {
        Field F(c_big_p);
        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Projective> P = E.random_point<Projective>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }
}

TEST(CorrectnessTest, ZeroProjective) {
    for (size_t i = 0; i < c_correctness_test_zero_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Projective> Z = E.null_point<Projective>();
        ASSERT_TRUE(Z.is_zero());
        FIELD_EQ(Z.get_x(), F.element(0));
        FIELD_EQ(Z.get_y(), F.element(1));
    }
}

TEST(CorrectnessTest, AdditionProjective) {
    for (size_t i = 0; i < c_correctness_test_addition_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<Projective> Z = E.null_point<Projective>();
        ASSERT_TRUE(Z.is_zero());

        EllipticCurvePoint<Projective> two_Z = Z + Z;
        ASSERT_TRUE(two_Z.is_zero());
        POINT_EQ(two_Z, Z);

        EllipticCurvePoint<Projective> P = E.random_point<Projective>();
        EllipticCurvePoint<Projective> Z_P = Z + P;
        POINT_EQ(Z_P, P);

        EllipticCurvePoint<Projective> Q = E.random_point<Projective>();
        EllipticCurvePoint<Projective> P_Q = P + Q;

        if (P.is_zero()) {
            POINT_EQ(P_Q, Q);
            continue;
        } else if (Q.is_zero()) {
            POINT_EQ(P_Q, P);
            continue;
        } else if (P == -Q) {
            ASSERT_TRUE(P_Q.is_zero());
            continue;
        } else if (P == Q) {
            const FieldElement& x = P.get_x();
            const FieldElement& y = P.get_y();

            if (!y.is_invertible()) {
                ASSERT_TRUE(P_Q.is_zero());
                continue;
            }

            FieldElement k = (F.element(3) * FieldElement::pow(x, 2) + a) / (y << 1);
            FieldElement correct_x = FieldElement::pow(k, 2) - (x << 1);
            FieldElement correct_y = k * (x - correct_x) - y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        } else {
            const FieldElement& lhs_x = P.get_x();
            const FieldElement& lhs_y = P.get_y();
            const FieldElement& rhs_x = Q.get_x();
            const FieldElement& rhs_y = Q.get_y();

            FieldElement k = (rhs_y - lhs_y) / (rhs_x - lhs_x);
            FieldElement correct_x = FieldElement::pow(k, 2) - lhs_x - rhs_x;
            FieldElement correct_y = k * (lhs_x - correct_x) - lhs_y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        }
    }
}

TEST(CorrectnessTest, kPProjective) {
    for (size_t i = 0; i < c_correctness_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<Projective> P = E.random_point<Projective>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<Projective> kP = k * P;

        EllipticCurvePoint<Projective> correct_kP = E.null_point<Projective>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<Projective> P = E.random_point<Projective>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<Projective> kP = k * P;

        EllipticCurvePoint<Projective> correct_kP = E.null_point<Projective>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }
}

// Stress tests
TEST(StressTest, kPProjective) {
    for (size_t i = 0; i < c_stress_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        uint p = primes::prime_number_list[pos];
        Field F(p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Projective> P = E.random_point<Projective>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<Projective> P = E.random_point<Projective>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }
}

// ModifiedJacobi Coordinates tests
// Correctness tests
TEST(CorrectnessTest, RandomModifiedJacobi) {
    for (size_t i = 0; i < c_correctness_test_random_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<ModifiedJacobi> P = E.random_point<ModifiedJacobi>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }

    {
        Field F(c_big_p);
        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<ModifiedJacobi> P = E.random_point<ModifiedJacobi>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }
}

TEST(CorrectnessTest, ZeroModifiedJacobi) {
    for (size_t i = 0; i < c_correctness_test_zero_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<ModifiedJacobi> Z = E.null_point<ModifiedJacobi>();
        ASSERT_TRUE(Z.is_zero());
        FIELD_EQ(Z.get_x(), F.element(0));
        FIELD_EQ(Z.get_y(), F.element(1));
    }
}

TEST(CorrectnessTest, AdditionModifiedJacobi) {
    for (size_t i = 0; i < c_correctness_test_addition_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<ModifiedJacobi> Z = E.null_point<ModifiedJacobi>();
        ASSERT_TRUE(Z.is_zero());

        EllipticCurvePoint<ModifiedJacobi> two_Z = Z + Z;
        ASSERT_TRUE(two_Z.is_zero());
        POINT_EQ(two_Z, Z);

        EllipticCurvePoint<ModifiedJacobi> P = E.random_point<ModifiedJacobi>();
        EllipticCurvePoint<ModifiedJacobi> Z_P = Z + P;
        POINT_EQ(Z_P, P);

        EllipticCurvePoint<ModifiedJacobi> Q = E.random_point<ModifiedJacobi>();
        EllipticCurvePoint<ModifiedJacobi> P_Q = P + Q;

        if (P.is_zero()) {
            POINT_EQ(P_Q, Q);
            continue;
        } else if (Q.is_zero()) {
            POINT_EQ(P_Q, P);
            continue;
        } else if (P == -Q) {
            ASSERT_TRUE(P_Q.is_zero());
            continue;
        } else if (P == Q) {
            const FieldElement& x = P.get_x();
            const FieldElement& y = P.get_y();

            if (!y.is_invertible()) {
                ASSERT_TRUE(P_Q.is_zero());
                continue;
            }

            FieldElement k = (F.element(3) * FieldElement::pow(x, 2) + a) / (y << 1);
            FieldElement correct_x = FieldElement::pow(k, 2) - (x << 1);
            FieldElement correct_y = k * (x - correct_x) - y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        } else {
            const FieldElement& lhs_x = P.get_x();
            const FieldElement& lhs_y = P.get_y();
            const FieldElement& rhs_x = Q.get_x();
            const FieldElement& rhs_y = Q.get_y();

            FieldElement k = (rhs_y - lhs_y) / (rhs_x - lhs_x);
            FieldElement correct_x = FieldElement::pow(k, 2) - lhs_x - rhs_x;
            FieldElement correct_y = k * (lhs_x - correct_x) - lhs_y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        }
    }
}

TEST(CorrectnessTest, kPModifiedJacobi) {
    for (size_t i = 0; i < c_correctness_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<ModifiedJacobi> P = E.random_point<ModifiedJacobi>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<ModifiedJacobi> kP = k * P;

        EllipticCurvePoint<ModifiedJacobi> correct_kP = E.null_point<ModifiedJacobi>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<ModifiedJacobi> P = E.random_point<ModifiedJacobi>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<ModifiedJacobi> kP = k * P;

        EllipticCurvePoint<ModifiedJacobi> correct_kP = E.null_point<ModifiedJacobi>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }
}

// Stress tests
TEST(StressTest, kPModifiedJacobi) {
    for (size_t i = 0; i < c_stress_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        uint p = primes::prime_number_list[pos];
        Field F(p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<ModifiedJacobi> P = E.random_point<ModifiedJacobi>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<ModifiedJacobi> P = E.random_point<ModifiedJacobi>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }
}

// JacobiChudnovski Coordinates tests
// Correctness tests
TEST(CorrectnessTest, RandomJacobiChudnovski) {
    for (size_t i = 0; i < c_correctness_test_random_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<JacobiChudnovski> P = E.random_point<JacobiChudnovski>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }

    {
        Field F(c_big_p);
        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<JacobiChudnovski> P = E.random_point<JacobiChudnovski>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }
}

TEST(CorrectnessTest, ZeroJacobiChudnovski) {
    for (size_t i = 0; i < c_correctness_test_zero_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<JacobiChudnovski> Z = E.null_point<JacobiChudnovski>();
        ASSERT_TRUE(Z.is_zero());
        FIELD_EQ(Z.get_x(), F.element(0));
        FIELD_EQ(Z.get_y(), F.element(1));
    }
}

TEST(CorrectnessTest, AdditionJacobiChudnovski) {
    for (size_t i = 0; i < c_correctness_test_addition_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<JacobiChudnovski> Z = E.null_point<JacobiChudnovski>();
        ASSERT_TRUE(Z.is_zero());

        EllipticCurvePoint<JacobiChudnovski> two_Z = Z + Z;
        ASSERT_TRUE(two_Z.is_zero());
        POINT_EQ(two_Z, Z);

        EllipticCurvePoint<JacobiChudnovski> P = E.random_point<JacobiChudnovski>();
        EllipticCurvePoint<JacobiChudnovski> Z_P = Z + P;
        POINT_EQ(Z_P, P);

        EllipticCurvePoint<JacobiChudnovski> Q = E.random_point<JacobiChudnovski>();
        EllipticCurvePoint<JacobiChudnovski> P_Q = P + Q;

        if (P.is_zero()) {
            POINT_EQ(P_Q, Q);
            continue;
        } else if (Q.is_zero()) {
            POINT_EQ(P_Q, P);
            continue;
        } else if (P == -Q) {
            ASSERT_TRUE(P_Q.is_zero());
            continue;
        } else if (P == Q) {
            const FieldElement& x = P.get_x();
            const FieldElement& y = P.get_y();

            if (!y.is_invertible()) {
                ASSERT_TRUE(P_Q.is_zero());
                continue;
            }

            FieldElement k = (F.element(3) * FieldElement::pow(x, 2) + a) / (y << 1);
            FieldElement correct_x = FieldElement::pow(k, 2) - (x << 1);
            FieldElement correct_y = k * (x - correct_x) - y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        } else {
            const FieldElement& lhs_x = P.get_x();
            const FieldElement& lhs_y = P.get_y();
            const FieldElement& rhs_x = Q.get_x();
            const FieldElement& rhs_y = Q.get_y();

            FieldElement k = (rhs_y - lhs_y) / (rhs_x - lhs_x);
            FieldElement correct_x = FieldElement::pow(k, 2) - lhs_x - rhs_x;
            FieldElement correct_y = k * (lhs_x - correct_x) - lhs_y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        }
    }
}

TEST(CorrectnessTest, kPJacobiChudnovski) {
    for (size_t i = 0; i < c_correctness_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<JacobiChudnovski> P = E.random_point<JacobiChudnovski>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<JacobiChudnovski> kP = k * P;

        EllipticCurvePoint<JacobiChudnovski> correct_kP = E.null_point<JacobiChudnovski>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<JacobiChudnovski> P = E.random_point<JacobiChudnovski>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<JacobiChudnovski> kP = k * P;

        EllipticCurvePoint<JacobiChudnovski> correct_kP = E.null_point<JacobiChudnovski>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }
}

// Stress tests
TEST(StressTest, kPJacobiChudnovski) {
    for (size_t i = 0; i < c_stress_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        uint p = primes::prime_number_list[pos];
        Field F(p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<JacobiChudnovski> P = E.random_point<JacobiChudnovski>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<JacobiChudnovski> P = E.random_point<JacobiChudnovski>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }
}

// SimplifiedJacobiChudnovski Coordinates tests
// Correctness tests
TEST(CorrectnessTest, RandomSimplifiedJacobiChudnovski) {
    for (size_t i = 0; i < c_correctness_test_random_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<SimplifiedJacobiChudnovski> P = E.random_point<SimplifiedJacobiChudnovski>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }

    {
        Field F(c_big_p);
        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);

        while (!a.is_invertible() && !b.is_invertible()) {
            a = generate_random_field_element(F);
            b = generate_random_field_element(F);
        }

        EllipticCurve E(a, b, F);
        EllipticCurvePoint<SimplifiedJacobiChudnovski> P = E.random_point<SimplifiedJacobiChudnovski>();

        if (P.is_zero()) {
            FIELD_EQ(P.get_x(), F.element(0));
            FIELD_EQ(P.get_y(), F.element(1));
        } else {
            FieldElement rhs = FieldElement::pow(P.get_y(), 2);
            FieldElement lhs = FieldElement::pow(P.get_x(), 3) + P.get_x() * a + b;
            FIELD_EQ(rhs, lhs);
        }
    }
}

TEST(CorrectnessTest, ZeroSimplifiedJacobiChudnovski) {
    for (size_t i = 0; i < c_correctness_test_zero_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<SimplifiedJacobiChudnovski> Z = E.null_point<SimplifiedJacobiChudnovski>();
        ASSERT_TRUE(Z.is_zero());
        FIELD_EQ(Z.get_x(), F.element(0));
        FIELD_EQ(Z.get_y(), F.element(1));
    }
}

TEST(CorrectnessTest, AdditionSimplifiedJacobiChudnovski) {
    for (size_t i = 0; i < c_correctness_test_addition_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<Normal> Z = E.null_point<Normal>();
        ASSERT_TRUE(Z.is_zero());

        EllipticCurvePoint<Normal> two_Z = Z + Z;
        ASSERT_TRUE(two_Z.is_zero());
        POINT_EQ(two_Z, Z);

        EllipticCurvePoint<Normal> P = E.random_point<Normal>();
        EllipticCurvePoint<Normal> Z_P = Z + P;
        POINT_EQ(Z_P, P);

        EllipticCurvePoint<Normal> Q = E.random_point<Normal>();
        EllipticCurvePoint<Normal> P_Q = P + Q;

        if (P.is_zero()) {
            POINT_EQ(P_Q, Q);
            continue;
        } else if (Q.is_zero()) {
            POINT_EQ(P_Q, P);
            continue;
        } else if (P == -Q) {
            ASSERT_TRUE(P_Q.is_zero());
            continue;
        } else if (P == Q) {
            const FieldElement& x = P.get_x();
            const FieldElement& y = P.get_y();

            if (!y.is_invertible()) {
                ASSERT_TRUE(P_Q.is_zero());
                continue;
            }

            FieldElement k = (F.element(3) * FieldElement::pow(x, 2) + a) / (y << 1);
            FieldElement correct_x = FieldElement::pow(k, 2) - (x << 1);
            FieldElement correct_y = k * (x - correct_x) - y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        } else {
            const FieldElement& lhs_x = P.get_x();
            const FieldElement& lhs_y = P.get_y();
            const FieldElement& rhs_x = Q.get_x();
            const FieldElement& rhs_y = Q.get_y();

            FieldElement k = (rhs_y - lhs_y) / (rhs_x - lhs_x);
            FieldElement correct_x = FieldElement::pow(k, 2) - lhs_x - rhs_x;
            FieldElement correct_y = k * (lhs_x - correct_x) - lhs_y;

            FIELD_EQ(P_Q.get_x(), correct_x);
            FIELD_EQ(P_Q.get_y(), correct_y);
        }
    }
}

TEST(CorrectnessTest, kPSimplifiedJacobiChudnovski) {
    for (size_t i = 0; i < c_correctness_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        Field F(primes::prime_number_list[pos]);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<SimplifiedJacobiChudnovski> P = E.random_point<SimplifiedJacobiChudnovski>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<SimplifiedJacobiChudnovski> kP = k * P;

        EllipticCurvePoint<SimplifiedJacobiChudnovski> correct_kP =
            E.null_point<SimplifiedJacobiChudnovski>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);

        EllipticCurvePoint<SimplifiedJacobiChudnovski> P = E.random_point<SimplifiedJacobiChudnovski>();
        uint k = generate_random_uint_modulo(c_correctness_max_k_n);
        EllipticCurvePoint<SimplifiedJacobiChudnovski> kP = k * P;

        EllipticCurvePoint<SimplifiedJacobiChudnovski> correct_kP =
            E.null_point<SimplifiedJacobiChudnovski>();

        for (uint j = 0; j < k; ++j) {
            correct_kP += P;
        }

        POINT_EQ(kP, correct_kP);
    }
}

// Stress tests
TEST(StressTest, kPSimplifiedJacobiChudnovski) {
    for (size_t i = 0; i < c_stress_test_kp_n; ++i) {
        size_t pos = generate_random_uint_modulo(c_primes_n).convert_to<size_t>();
        uint p = primes::prime_number_list[pos];
        Field F(p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<SimplifiedJacobiChudnovski> P = E.random_point<SimplifiedJacobiChudnovski>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }

    {
        Field F(c_big_p);

        FieldElement a = generate_random_field_element(F);
        FieldElement b = generate_random_field_element(F);
        EllipticCurve E(a, b, F);
        EllipticCurvePoint<SimplifiedJacobiChudnovski> P = E.random_point<SimplifiedJacobiChudnovski>();
        uint k = generate_random_uint_modulo(c_stress_max_k_n);
        auto kP = P * k;
    }
}
