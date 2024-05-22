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

//TEST(SimpleTest, GCD) {
//    uint p = 7;
//    Field F(p);
//    Poly a(F, {3, 3});
//    Poly b(F, {3});
//    auto gcd = algorithm::gcd(a, b);
//    Poly ans(F, {3});
//    ASSERT_EQ(gcd, ans);
//
//    a = Poly(F, {-F.element(1), F.element(0), F.element(1)});
//    b = Poly(F, {-F.element(1), F.element(1)});
//    gcd = algorithm::gcd(a, b);
//    ans = Poly(F, {-F.element(1), F.element(1)});
//    ASSERT_EQ(gcd, ans);
//}
//
//TEST(SimpleTest, Power) {
//    uint p = 7;
//    Field F(p);
//    Poly a(F, {3, 3});
//    Poly b = Poly::pow(a, 4);
//    Poly ans(F, {81, 324, 486, 324, 81});
//    ASSERT_EQ(b, ans);
//}
//
//TEST(SimpleTest, Counting1) {
//    uint p = 7;
//    Field F(p);
//    FieldElement a = F.element(2);
//    FieldElement b = F.element(1);
//    EllipticCurve E(a, b, F);
//    uint points_number = algorithm::schoof::points_number(E);
//    ASSERT_EQ(points_number, 5);
//}
//
//TEST(SimpleTest, Counting2) {
//    uint p = 17;
//    Field F(p);
//    FieldElement a = F.element(5);
//    FieldElement b = F.element(2);
//    EllipticCurve E(a, b, F);
//    uint points_number = algorithm::schoof::points_number(E);
//    ASSERT_EQ(points_number, 15);
//}
//
//TEST(CorrectnessTest, Counting) {
//    uint p = "0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff";
//
//    Field F(p);
//    FieldElement a = F.element("0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc");
//    FieldElement b = F.element("0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");
//
//    FieldElement G_x = F.element("0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296");
//    FieldElement G_y = F.element("0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");
//
//    EllipticCurve E(a, b, F);
//    uint points_number = algorithm::schoof::points_number(E);
//    uint correct_number = "0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551";
//    ASSERT_EQ(points_number, correct_number);
//}
