// clang-format off
#include "pch.h"
// clang-format on

#include "elliptic-curve.h"

using namespace elliptic_curve_guide;
using namespace field;
using namespace elliptic_curve;
using enum CoordinatesType;

static constexpr uint c_find_y_n = 53617;   // e = 4
static const uint c_kp_timing_n("0xFFFFFFFFF00FFFFFFFFFFF1FFFFFFFFFFFFFF9FFFFAFFFFFFFFFFFFFFFF");
static constexpr size_t c_kp_comparison_timing_n = 228133;
static constexpr size_t c_doubling_timing_n = 228133;

TEST(CorrectnessTest, FindY) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    auto opt = E.point_with_x_equal_to(F.element(4));
    ASSERT_TRUE(opt.has_value());
    ASSERT_EQ(opt.value().get_y().value(), 18);
}

TEST(StressTest, FindMultipleY) {
    Field F(c_find_y_n);
    static const FieldElement one = F.element(1);
    FieldElement a = F.element(1984);
    FieldElement b = F.element(1337);
    EllipticCurve E(a, b, F);
    uint point_number = 0;

    for (uint i = 0; i < c_find_y_n; ++i) {
        FieldElement x = F.element(i);
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
            ASSERT_EQ(y2, value);
            ++point_number;
        } else {
            value.pow((c_find_y_n - 1) >> 1);
            ASSERT_NE(value, one);
        }
    }

    ASSERT_GE(point_number, 1);

    a = F.element(81741);
    b = F.element(42);
    E = EllipticCurve(a, b, F);
    point_number = 0;

    for (uint i = 0; i < c_find_y_n; ++i) {
        FieldElement x = F.element(i);
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
            ASSERT_EQ(y2, value);
            ++point_number;
        } else {
            value.pow((c_find_y_n - 1) >> 1);
            ASSERT_NE(value, one);
        }
    }

    ASSERT_GE(point_number, 1);
}

// Normal
TEST(NormalCorrectnessTest, Creating) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    FieldElement x = F.element(4);
    EllipticCurvePoint<Normal> P = E.point_with_x_equal_to<Normal>(x).value();
    ASSERT_EQ(P.get_x(), x);
    P = E.random_point();
}

TEST(NormalCorrectnessTest, Zero) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    auto Z = E.null_point<Normal>();
    ASSERT_TRUE(Z.is_zero());
    auto two_Z = Z + Z;
    ASSERT_TRUE(two_Z.is_zero());
}

TEST(NormalCorrectnessTest, Addition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Normal> P = E.point_with_x_equal_to<Normal>(F.element(4)).value();
    auto Z = E.null_point<Normal>();
    ASSERT_EQ(P + Z, P);
    auto twoP = P + P;
    auto newP = twoP - P;
    ASSERT_EQ(newP, P);
}

TEST(NormalCorrectnessTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Normal> P = E.point_with_x_equal_to<Normal>(F.element(4)).value();
    size_t k = 1985;
    auto kP = P * k;
    EllipticCurvePoint<Normal> right_kP = P;
    for (size_t i = 1; i < k; ++i) {
        right_kP += P;
    }
    ASSERT_EQ(kP, right_kP);
}

TEST(NormalTimingTest, kPbyMultiplying) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Normal> P = E.point_with_x_equal_to<Normal>(F.element(4)).value();
    auto kP = P * c_kp_comparison_timing_n;
}

TEST(NormalTimingTest, kPbyAddition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Normal> P = E.point_with_x_equal_to<Normal>(F.element(4)).value();
    EllipticCurvePoint<Normal> kP = P;
    for (size_t i = 1; i < c_kp_comparison_timing_n; ++i) {
        kP += P;
    }
}

TEST(NormalTimingTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Normal> P = E.point_with_x_equal_to<Normal>(F.element(4)).value();
    auto kP = P * c_kp_timing_n;
}

TEST(NormalTimingTest, Doubling) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Normal> P = E.point_with_x_equal_to<Normal>(F.element(4)).value();

    for (size_t i = 0; i < c_doubling_timing_n; ++i) {
        P += P;
    }
}

// Projective
TEST(ProjectiveCorrectnessTest, Creating) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    FieldElement x = F.element(4);
    EllipticCurvePoint<Projective> P = E.point_with_x_equal_to<Projective>(x).value();
    ASSERT_EQ(P.get_x(), x);
}

TEST(ProjectiveCorrectnessTest, Zero) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    auto Z = E.null_point<Projective>();
    ASSERT_TRUE(Z.is_zero());
    auto two_Z = Z + Z;
    ASSERT_TRUE(two_Z.is_zero());
}

TEST(ProjectiveCorrectnessTest, Addition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Projective> P = E.point_with_x_equal_to<Projective>(F.element(4)).value();
    auto Z = E.null_point<Projective>();
    ASSERT_EQ(P + Z, P);
    auto twoP = P + P;
    auto newP = twoP - P;
    ASSERT_EQ(newP, P);
}

TEST(ProjectiveCorrectnessTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Projective> P = E.point_with_x_equal_to<Projective>(F.element(4)).value();
    size_t k = 1985;
    auto kP = P * k;
    EllipticCurvePoint<Projective> right_kP = P;
    for (size_t i = 1; i < k; ++i) {
        right_kP += P;
    }
    ASSERT_EQ(kP, right_kP);
}

TEST(ProjectiveTimingTest, kPbyMultiplying) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Projective> P = E.point_with_x_equal_to<Projective>(F.element(4)).value();
    auto kP = P * c_kp_comparison_timing_n;
}

TEST(ProjectiveTimingTest, kPbyAddition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Projective> P = E.point_with_x_equal_to<Projective>(F.element(4)).value();
    EllipticCurvePoint<Projective> kP = P;
    for (size_t i = 1; i < c_kp_comparison_timing_n; ++i) {
        kP += P;
    }
}

TEST(ProjectiveTimingTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Projective> P = E.point_with_x_equal_to<Projective>(F.element(4)).value();
    auto kP = P * c_kp_timing_n;
}

TEST(ProjectiveTimingTest, Doubling) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Projective> P = E.point_with_x_equal_to<Projective>(F.element(4)).value();

    for (size_t i = 0; i < c_doubling_timing_n; ++i) {
        P += P;
    }
}

// Jacobi
TEST(JacobiCorrectnessTest, Creating) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    FieldElement x = F.element(4);
    EllipticCurvePoint<Jacobi> P = E.point_with_x_equal_to<Jacobi>(x).value();
    ASSERT_EQ(P.get_x(), x);
}

TEST(JacobiCorrectnessTest, Zero) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    auto Z = E.null_point<Jacobi>();
    ASSERT_TRUE(Z.is_zero());
    auto two_Z = Z + Z;
    ASSERT_TRUE(two_Z.is_zero());
}

TEST(JacobiCorrectnessTest, Addition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Jacobi> P = E.point_with_x_equal_to<Jacobi>(F.element(4)).value();
    auto Z = E.null_point<Jacobi>();
    ASSERT_EQ(P + Z, P);
    auto twoP = P + P;
    auto newP = twoP - P;
    ASSERT_EQ(newP, P);
}

TEST(JacobiCorrectnessTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Jacobi> P = E.point_with_x_equal_to<Jacobi>(F.element(4)).value();
    size_t k = 1985;
    auto kP = P * k;
    EllipticCurvePoint<Jacobi> right_kP = P;
    for (size_t i = 1; i < k; ++i) {
        right_kP += P;
    }
    ASSERT_EQ(kP, right_kP);
}

TEST(JacobiTimingTest, kPbyMultiplying) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Jacobi> P = E.point_with_x_equal_to<Jacobi>(F.element(4)).value();
    auto kP = P * c_kp_comparison_timing_n;
}

TEST(JacobiTimingTest, kPbyAddition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Jacobi> P = E.point_with_x_equal_to<Jacobi>(F.element(4)).value();
    EllipticCurvePoint<Jacobi> kP = P;
    for (size_t i = 1; i < c_kp_comparison_timing_n; ++i) {
        kP += P;
    }
}

TEST(JacobiTimingTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Jacobi> P = E.point_with_x_equal_to<Jacobi>(F.element(4)).value();
    auto kP = P * c_kp_timing_n;
}

TEST(JacobiTimingTest, Doubling) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Jacobi> P = E.point_with_x_equal_to<Jacobi>(F.element(4)).value();

    for (size_t i = 0; i < c_doubling_timing_n; ++i) {
        P += P;
    }
}

// ModifiedJacobi
TEST(ModifiedJacobiCorrectnessTest, Creating) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    FieldElement x = F.element(4);
    EllipticCurvePoint<ModifiedJacobi> P = E.point_with_x_equal_to<ModifiedJacobi>(x).value();
    ASSERT_EQ(P.get_x(), x);
}

TEST(ModifiedJacobiCorrectnessTest, Zero) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    auto Z = E.null_point<ModifiedJacobi>();
    ASSERT_TRUE(Z.is_zero());
    auto two_Z = Z + Z;
    ASSERT_TRUE(two_Z.is_zero());
}

TEST(ModifiedJacobiCorrectnessTest, Addition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<ModifiedJacobi> P = E.point_with_x_equal_to<ModifiedJacobi>(F.element(4)).value();
    auto Z = E.null_point<ModifiedJacobi>();
    ASSERT_EQ(P + Z, P);
    auto twoP = P + P;
    auto newP = twoP - P;
    ASSERT_EQ(newP, P);
}

TEST(ModifiedJacobiCorrectnessTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<ModifiedJacobi> P = E.point_with_x_equal_to<ModifiedJacobi>(F.element(4)).value();
    size_t k = 1985;
    auto kP = P * k;
    EllipticCurvePoint<ModifiedJacobi> right_kP = P;
    for (size_t i = 1; i < k; ++i) {
        right_kP += P;
    }
    ASSERT_EQ(kP, right_kP);
}

TEST(ModifiedJacobiTimingTest, kPbyMultiplying) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<ModifiedJacobi> P = E.point_with_x_equal_to<ModifiedJacobi>(F.element(4)).value();
    auto kP = P * c_kp_comparison_timing_n;
}

TEST(ModifiedJacobiTimingTest, kPbyAddition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<ModifiedJacobi> P = E.point_with_x_equal_to<ModifiedJacobi>(F.element(4)).value();
    EllipticCurvePoint<ModifiedJacobi> kP = P;
    for (size_t i = 1; i < c_kp_comparison_timing_n; ++i) {
        kP += P;
    }
}

TEST(ModifiedJacobiTimingTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<ModifiedJacobi> P = E.point_with_x_equal_to<ModifiedJacobi>(F.element(4)).value();
    auto kP = P * c_kp_timing_n;
}

TEST(ModifiedJacobiTimingTest, Doubling) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<ModifiedJacobi> P = E.point_with_x_equal_to<ModifiedJacobi>(F.element(4)).value();

    for (size_t i = 0; i < c_doubling_timing_n; ++i) {
        P += P;
    }
}

// JacobiChudnovski
TEST(JacobiChudnovskiCorrectnessTest, Creating) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    FieldElement x = F.element(4);
    EllipticCurvePoint<JacobiChudnovski> P = E.point_with_x_equal_to<JacobiChudnovski>(x).value();
    ASSERT_EQ(P.get_x(), x);
}

TEST(JacobiChudnovskiCorrectnessTest, Zero) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    auto Z = E.null_point<JacobiChudnovski>();
    ASSERT_TRUE(Z.is_zero());
    auto two_Z = Z + Z;
    ASSERT_TRUE(two_Z.is_zero());
}

TEST(JacobiChudnovskiCorrectnessTest, Addition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<JacobiChudnovski> P = E.point_with_x_equal_to<JacobiChudnovski>(F.element(4)).value();
    auto Z = E.null_point<JacobiChudnovski>();
    ASSERT_EQ(P + Z, P);
    auto twoP = P + P;
    auto newP = twoP - P;
    ASSERT_EQ(newP, P);
}

TEST(JacobiChudnovskiCorrectnessTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<JacobiChudnovski> P = E.point_with_x_equal_to<JacobiChudnovski>(F.element(4)).value();
    size_t k = 1985;
    auto kP = P * k;
    EllipticCurvePoint<JacobiChudnovski> right_kP = P;
    for (size_t i = 1; i < k; ++i) {
        right_kP += P;
    }
    ASSERT_EQ(kP, right_kP);
}

TEST(JacobiChudnovskiTimingTest, kPbyMultiplying) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<JacobiChudnovski> P = E.point_with_x_equal_to<JacobiChudnovski>(F.element(4)).value();
    auto kP = P * c_kp_comparison_timing_n;
}

TEST(JacobiChudnovskiTimingTest, kPbyAddition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<JacobiChudnovski> P = E.point_with_x_equal_to<JacobiChudnovski>(F.element(4)).value();
    EllipticCurvePoint<JacobiChudnovski> kP = P;
    for (size_t i = 1; i < c_kp_comparison_timing_n; ++i) {
        kP += P;
    }
}

TEST(JacobiChudnovskiTimingTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<JacobiChudnovski> P = E.point_with_x_equal_to<JacobiChudnovski>(F.element(4)).value();
    auto kP = P * c_kp_timing_n;
}

TEST(JacobiChudnovskiTimingTest, Doubling) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<JacobiChudnovski> P = E.point_with_x_equal_to<JacobiChudnovski>(F.element(4)).value();

    for (size_t i = 0; i < c_doubling_timing_n; ++i) {
        P += P;
    }
}

// SimplifiedJacobiChudnovski
TEST(SimplifiedJacobiChudnovskiCorrectnessTest, Creating) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    FieldElement x = F.element(4);
    EllipticCurvePoint<SimplifiedJacobiChudnovski> P =
        E.point_with_x_equal_to<SimplifiedJacobiChudnovski>(x).value();
    ASSERT_EQ(P.get_x(), x);
}

TEST(SimplifiedJacobiChudnovskiCorrectnessTest, Zero) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    auto Z = E.null_point<SimplifiedJacobiChudnovski>();
    ASSERT_TRUE(Z.is_zero());
    auto two_Z = Z + Z;
    ASSERT_TRUE(two_Z.is_zero());
}

TEST(SimplifiedJacobiChudnovskiCorrectnessTest, Addition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<SimplifiedJacobiChudnovski> P =
        E.point_with_x_equal_to<SimplifiedJacobiChudnovski>(F.element(4)).value();
    auto Z = E.null_point<SimplifiedJacobiChudnovski>();
    ASSERT_EQ(P + Z, P);
    auto twoP = P + P;
    auto newP = twoP - P;
    ASSERT_EQ(newP, P);
}

TEST(SimplifiedJacobiChudnovskiCorrectnessTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<SimplifiedJacobiChudnovski> P =
        E.point_with_x_equal_to<SimplifiedJacobiChudnovski>(F.element(4)).value();
    size_t k = 1985;
    auto kP = P * k;
    EllipticCurvePoint<SimplifiedJacobiChudnovski> right_kP = P;
    for (size_t i = 1; i < k; ++i) {
        right_kP += P;
    }
    ASSERT_EQ(kP, right_kP);
}

TEST(SimplifiedJacobiChudnovskiTimingTest, kPbyMultiplying) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<SimplifiedJacobiChudnovski> P =
        E.point_with_x_equal_to<SimplifiedJacobiChudnovski>(F.element(4)).value();
    auto kP = P * c_kp_comparison_timing_n;
}

TEST(SimplifiedJacobiChudnovskiTimingTest, kPbyAddition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<SimplifiedJacobiChudnovski> P =
        E.point_with_x_equal_to<SimplifiedJacobiChudnovski>(F.element(4)).value();
    EllipticCurvePoint<SimplifiedJacobiChudnovski> kP = P;
    for (size_t i = 1; i < c_kp_comparison_timing_n; ++i) {
        kP += P;
    }
}

TEST(SimplifiedJacobiChudnovskiTimingTest, kP) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<SimplifiedJacobiChudnovski> P =
        E.point_with_x_equal_to<SimplifiedJacobiChudnovski>(F.element(4)).value();
    auto kP = P * c_kp_timing_n;
}

TEST(SimplifiedJacobiChudnovskiTimingTest, Doubling) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<SimplifiedJacobiChudnovski> P =
        E.point_with_x_equal_to<SimplifiedJacobiChudnovski>(F.element(4)).value();

    for (size_t i = 0; i < c_doubling_timing_n; ++i) {
        P += P;
    }
}
