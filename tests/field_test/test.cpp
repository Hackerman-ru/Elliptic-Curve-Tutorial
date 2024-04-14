#include "pch.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/fwd.hpp>
#include <field.h>

using namespace ECG;

TEST(SimpleTesting, Creating) {
    uint p("7");
    Field f(p);
    FieldElement a = f(3);
    FieldElement b = f(4);
    uint c = (a * b).value();
    int result = c.convert_to<int>();
    ASSERT_EQ(result, 5);
}
