#include "pch.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/fwd.hpp>
#include <field.h>

using namespace ECG;

TEST(SimpleTesting, Creating) {
    uint p("7");
    Field f(p);
    FieldElement a = f(10);
    ASSERT_EQ(a.value(), 3);
    a = f(99);
    ASSERT_EQ(a.value(), 1);
}

TEST(SimpleTesting, Arithmetic) {
    uint p("7");
    Field f(p);
    FieldElement a = f(3);
    FieldElement b = a;
    uint c = (a * b).value();
    ASSERT_EQ(c, 2);
}

TEST(SimpleTesting, ArithmeticBig) {
    uint p = uint("1000000007");
    Field g(p);
    uint a_("999999999");
    uint b_ = 2;
    FieldElement a = g(a_);
    FieldElement b = g(b_);
    FieldElement inverse = a * b;
    inverse.inverse();
    uint result = inverse.value();
    uint correct_result("437500003");
    ASSERT_EQ(result, correct_result);
}
