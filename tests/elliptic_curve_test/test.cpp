#include "pch.h"

#include <elliptic-curve.h>

using namespace ECG;
static constexpr uint find_y_prime = 99371;

TEST(SimpleTest, Creating) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<> P = E.point_with_x_equal_to(F.element(4)).value();
}

TEST(SimpleTest, Zero) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(a, b, F);
    EllipticCurvePoint<> P = E.null_point();
    ASSERT_TRUE(P.is_zero());
}

TEST(SimpleTest, FindY) {
    Field F(29);
    FieldElement a = F.element(0);
    FieldElement b = F.element(28);
    EllipticCurve E(a, b, F);
    auto opt = E.point_with_x_equal_to(F.element(4));
    ASSERT_TRUE(opt.has_value());
    ASSERT_EQ(opt.value().get_y().value(), 18);
}

TEST(HardTest, FindMultipleY) {
    Field F(find_y_prime);
    static const FieldElement one = F.element(1);
    FieldElement a = F.element(1984);
    FieldElement b = F.element(1337);
    EllipticCurve E(a, b, F);
    uint point_number = 0;

    for (uint i = 0; i < find_y_prime; ++i) {
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
            value.pow((find_y_prime - 1) >> 1);
            ASSERT_NE(value, one);
        }
    }

    ASSERT_GE(point_number, 1);

    a = F.element(81741);
    b = F.element(42);
    E = EllipticCurve(a, b, F);
    point_number = 0;

    for (uint i = 0; i < find_y_prime; ++i) {
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
            value.pow((find_y_prime - 1) >> 1);
            ASSERT_NE(value, one);
        }
    }

    ASSERT_GE(point_number, 1);
}
