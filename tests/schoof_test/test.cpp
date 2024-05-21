// clang-format off
#include "pch.h"
// clang-format on

#include "elliptic-curve.h"
#include "field.h"
#include "uint.h"
#include "utils/schoof.h"

using namespace elliptic_curve_guide;
using namespace field;
using namespace elliptic_curve;

TEST(SimpleTest, Counting) {
    uint p = 7;
    Field F(p);
    FieldElement a = F.element(2);
    FieldElement b = F.element(1);
    EllipticCurve E(a, b, F);
    uint points_number = algorithm::schoof::points_number(E);
    std::string str = points_number.convert_to<std::string>();
}
