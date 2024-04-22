#include "pch.h"

#include <elliptic-curve.h>

using namespace ECG;

TEST(SimpleTest, Creating) {
    Field F(1000000007);
    FieldElement a = F.element(0);
    FieldElement b = F.element(17);
    EllipticCurve E(std::move(a), std::move(b), F);
    EllipticCurvePoint<> P = E.point_with_x_equal_to(F.element(4)).value();
}
