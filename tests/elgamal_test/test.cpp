// clang-format off
#include "pch.h"
// clang-format on

#include "el-gamal.h"
#include "elliptic-curve.h"
#include "field.h"
#include "utils/random.h"

using namespace elliptic_curve_guide;
using namespace field;
using namespace elliptic_curve;
using namespace algorithm::encryption;
using namespace algorithm::random;

using Curve = EllipticCurve;

#define UINT_EQ(lhs, rhs)                                    \
    if (lhs != rhs) {                                        \
        std::string lhs_str = lhs.convert_to<std::string>(); \
        std::string rhs_str = rhs.convert_to<std::string>(); \
        ASSERT_EQ(lhs_str, rhs_str);                         \
    }

constexpr uint p = "0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff";

const Field F(p);
const FieldElement a = F.element("0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc");
const FieldElement b = F.element("0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b");

const FieldElement G_x = F.element("0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296");
const FieldElement G_y = F.element("0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5");

const Curve E(a, b, F);
const ElGamal::Point G = E.point<ElGamal::point_type>(G_x, G_y).value();
const uint n = "0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551";

const ElGamal EG(E, G, n);
ElGamal::Keys keys = EG.generate_keys();

static constexpr uint c_message_mask = (uint(1) << 128) - 1;
static constexpr size_t c_correctness_test_encryption_n = 100;
static constexpr size_t c_stress_test_encryption_n = 1000;

TEST(SimpleTest, Encryption) {
    ElGamal::Keys keys = EG.generate_keys();
    uint message = "0xFFF12341ABCBFFBBBE";
    ElGamal::EncryptedMessage enc = EG.encrypt(message, keys.public_key);
    uint decrypted_message = EG.decrypt_to_uint(enc, keys.private_key);
    UINT_EQ(message, decrypted_message);
}

TEST(CorrectnessTest, Encryption) {
    ElGamal::Keys keys = EG.generate_keys();

    for (size_t i = 0; i < c_correctness_test_encryption_n; ++i) {
        uint message = generate_random_uint() & c_message_mask;
        ElGamal::EncryptedMessage enc = EG.encrypt(message, keys.public_key);
        uint decrypted_message = EG.decrypt_to_uint(enc, keys.private_key);
        UINT_EQ(message, decrypted_message);
    }
}

TEST(StressTest, Encryption) {
    ElGamal::Keys keys = EG.generate_keys();

    for (size_t i = 0; i < c_stress_test_encryption_n; ++i) {
        uint message = generate_random_uint() & c_message_mask;
        ElGamal::EncryptedMessage enc = EG.encrypt(message, keys.public_key);
        uint decrypted_message = EG.decrypt_to_uint(enc, keys.private_key);
    }
}
