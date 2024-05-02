#include "pch.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/fwd.hpp>
#include <field.h>

using namespace ECG;

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
    Field g("1000000007");
    FieldElement a = g.element("999999999");
    FieldElement b = g.element(2);
    FieldElement inverse = a * b;
    inverse.inverse();
    uint result = inverse.value();
    uint correct_result("437500003");
    ASSERT_EQ(result, correct_result);
}

TEST(SimpleTesting, Comparison) {
    Field g("1000000007");
    FieldElement a = g.element("999999999");
    FieldElement b = g.element(2);
    ASSERT_EQ(a < b, false);
}

TEST(SimpleTesting, Shift) {
    Field g("1000000007");
    FieldElement a = g.element("999999999");
    FieldElement b = a << uint("30");
    uint result = b.value();
    uint correct_result("410065471");
    ASSERT_EQ(result, correct_result);
}
