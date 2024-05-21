// clang-format off
#include "pch.h"
// clang-format on

#include "elliptic-curve.h"
#include "field.h"
#include "uint.h"
#include "utils/schoof/polynomial.h"
#include "utils/schoof/schoof.h"

using namespace elliptic_curve_guide;
using namespace field;
using namespace elliptic_curve;
using namespace polynomial;

TEST(SimpleTest, GCD) {
    uint p = 7;
    Field F(p);
    Poly a(F, {3, 3});
    Poly b(F, {3});
    auto gcd = algorithm::gcd(a, b);
    Poly ans(F, {3});
    ASSERT_EQ(gcd, ans);

    a = Poly(F, {-F.element(1), F.element(0), F.element(1)});
    b = Poly(F, {-F.element(1), F.element(1)});
    gcd = algorithm::gcd(a, b);
    ans = Poly(F, {-F.element(1), F.element(1)});
    ASSERT_EQ(gcd, ans);
}

TEST(SimpleTest, Power) {
    uint p = 7;
    Field F(p);
    Poly a(F, {3, 3});
    Poly b = Poly::pow(a, 4);
    Poly ans(F, {81, 324, 486, 324, 81});
    ASSERT_EQ(b, ans);
}

TEST(SimpleTest, Counting) {
    uint p = 7;
    Field F(p);
    FieldElement a = F.element(2);
    FieldElement b = F.element(1);
    EllipticCurve E(a, b, F);
    uint points_number = algorithm::schoof::points_number(E);
    std::string str = points_number.convert_to<std::string>();
}
