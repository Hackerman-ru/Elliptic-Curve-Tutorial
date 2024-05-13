#include "pch.h"

#include <ecdsa.h>
#include <elliptic-curve.h>
#include <field.h>

using namespace elliptic_curve_guide;
using namespace field;
using namespace elliptic_curve;
using namespace algorithm::encryption;

TEST(SimpleTest, Creation) {
    const uint p("0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff");
    const Field F(p);
    const uint a("0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc");
    const uint b("0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");
    const EllipticCurve E(F.element(a), F.element(b), F);
    const FieldElement G_x = F.element("0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296");
    const FieldElement G_y = F.element("0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");
    auto opt = E.point(G_x, G_y);
    ASSERT_TRUE(opt.has_value());
    const EllipticCurvePoint<>& G = opt.value();
    const uint n("0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551");
    const uint h("0x1");
    ECDSA sh(F, E, G, n, h);
}

TEST(SimpleTest, SignGeneration) {
    const uint p("0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff");
    const Field F(p);
    const uint a("0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc");
    const uint b("0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");
    const EllipticCurve E(F.element(a), F.element(b), F);
    const FieldElement G_x = F.element("0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296");
    const FieldElement G_y = F.element("0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");
    auto opt = E.point(G_x, G_y);
    ASSERT_TRUE(opt.has_value());
    const EllipticCurvePoint<>& G = opt.value();
    const uint n("0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551");
    const uint h("0x1");
    ECDSA sh(F, E, G, n, h);
    const auto keys = sh.generate_keys();
    const uint message("1234567890");
    const auto sig = sh.generate_signature(keys.private_key, message);
}

TEST(SimpleTest, SignCorrectness) {
    const uint p("0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff");
    const Field F(p);
    const uint a("0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc");
    const uint b("0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");
    const EllipticCurve E(F.element(a), F.element(b), F);
    const FieldElement G_x = F.element("0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296");
    const FieldElement G_y = F.element("0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");
    auto opt = E.point(G_x, G_y);
    ASSERT_TRUE(opt.has_value());
    const EllipticCurvePoint<>& G = opt.value();
    const uint n("0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551");
    const uint h("0x1");
    ECDSA sh(F, E, G, n, h);
    const auto keys = sh.generate_keys();
    const uint message("1234567890");
    const auto sig = sh.generate_signature(keys.private_key, message);
    ASSERT_TRUE(sh.is_correct_signature(keys.public_key, message, sig));
}
