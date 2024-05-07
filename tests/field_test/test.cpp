// clang-format off
#include "pch.h"
// clang-format on
#include "field.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/fwd.hpp>

using namespace ECG;

static const uint inversion_hard_n = 1599827;

TEST(SimpleTesting, Creating) {
    Field f("7");
    FieldElement a = f.element(10);
    ASSERT_EQ(a.value(), 3);
    a = f.element(99);
    ASSERT_EQ(a.value(), 1);
}

TEST(SimpleTesting, Multiplication) {
    Field f("7");
    FieldElement a = f.element(3);
    FieldElement b = a;
    uint c = (a * b).value();
    ASSERT_EQ(c, 2);
}

TEST(SimpleTesting, Negotiation) {
    Field f("7");
    FieldElement a = -f.element(3);
    uint c = FieldElement::pow(a, 10).value();
    ASSERT_EQ(c, 4);
}

TEST(SimpleTesting, Inversion) {
    Field f("1000000007");
    FieldElement a = f.element("999999999");
    FieldElement b = f.element(2);
    FieldElement inverse = a * b;
    inverse.inverse();
    uint result = inverse.value();
    uint correct_result("437500003");
    ASSERT_EQ(result, correct_result);
}

TEST(SimpleTesting, Comparison) {
    Field f("1000000007");
    FieldElement a = f.element("999999999");
    FieldElement b = f.element(2);
    ASSERT_EQ(a < b, false);
}

TEST(SimpleTesting, Shift) {
    Field f("1000000007");
    FieldElement a = f.element("999999999");
    FieldElement b = a << uint("30");
    uint result = b.value();
    uint correct_result("410065471");
    ASSERT_EQ(result, correct_result);
}

TEST(HardTest, Inversion) {
    Field f(inversion_hard_n);
    const FieldElement one = f.element(1);
    for (uint i = 1; i < inversion_hard_n; ++i) {
        FieldElement a = f.element(i);
        FieldElement inv_a = FieldElement::inverse(a);
        FieldElement result = a * inv_a;
        ASSERT_EQ(result, one);
    }
}
