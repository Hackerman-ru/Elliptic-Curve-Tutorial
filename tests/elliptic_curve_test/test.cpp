#include "pch.h"

#include <elliptic-curve.h>

using namespace ECG;
using enum CoordinatesType;

static constexpr uint c_find_y_n = 53617;   // e = 4
static constexpr size_t c_kp_timing_n = 2281337;

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
        FieldElement value = FieldElement::pow(x, 3) + a * x + b;
        auto opt = E.point_with_x_equal_to(x);

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
        FieldElement value = FieldElement::pow(x, 3) + a * x + b;
        auto opt = E.point_with_x_equal_to(x);

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
    auto kP = P * c_kp_timing_n;
}

TEST(NormalTimingTest, kPbyAddition) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<Normal> P = E.point_with_x_equal_to<Normal>(F.element(4)).value();
    EllipticCurvePoint<Normal> kP = P;
    for (size_t i = 1; i < c_kp_timing_n; ++i) {
        kP += P;
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

//TEST(ProjectiveCorrectnessTest, kP) {
//    Field F(29);
//    FieldElement a = F.element(0);
//    FieldElement b = F.element(28);
//    EllipticCurve E(a, b, F);
//    EllipticCurvePoint<Projective> P = E.point_with_x_equal_to<Projective>(F.element(4)).value();
//    size_t k = 1985;
//    auto kP = P * k;
//    EllipticCurvePoint<Projective> right_kP = P;
//    for (size_t i = 1; i < k; ++i) {
//        right_kP += P;
//    }
//    ASSERT_EQ(kP, right_kP);
//}
//
//TEST(ProjectiveTimingTest, kPbyMultiplying) {
//    Field F(29);
//    FieldElement a = F.element(0);
//    FieldElement b = F.element(28);
//    EllipticCurve E(a, b, F);
//    EllipticCurvePoint<Projective> P = E.point_with_x_equal_to<Projective>(F.element(4)).value();
//    auto kP = P * c_kp_timing_n;
//}
//
//TEST(ProjectiveTimingTest, kPbyAddition) {
//    Field F(29);
//    FieldElement a = F.element(0);
//    FieldElement b = F.element(28);
//    EllipticCurve E(a, b, F);
//    EllipticCurvePoint<Projective> P = E.point_with_x_equal_to<Projective>(F.element(4)).value();
//    EllipticCurvePoint<Projective> kP = P;
//    for (size_t i = 1; i < c_kp_timing_n; ++i) {
//        kP += P;
//    }
//}
