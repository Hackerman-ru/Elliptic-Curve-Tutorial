#ifndef ECG_EL_GAMAL_H
#define ECG_EL_GAMAL_H

#include "elliptic-curve.h"
#include "functional"

namespace elliptic_curve_guide {
    namespace algorithm {
        namespace encryption {
            class ElGamal {
                using Curve = elliptic_curve::EllipticCurve;
                using Point = elliptic_curve::EllipticCurvePoint<elliptic_curve::CoordinatesType::Normal>;

            public:
                struct Keys {
                    uint private_key;
                    Point public_key;
                };

                enum class EncryptionType {
                    Standard,
                    Hashed,
                };

                template<EncryptionType type = EncryptionType::Standard>
                struct EncryptedMessage;

                template<>
                struct EncryptedMessage<EncryptionType::Standard> {
                    Point generator_degree;
                    Point message_with_salt;
                };

                template<>
                struct EncryptedMessage<EncryptionType::Hashed> {
                    Point generator_degree;
                    uint message_with_salt;
                };

                ElGamal(const Curve& curve, const Point& generator, const uint& generator_order);

                Keys generate_keys() const;

                EncryptedMessage<EncryptionType::Standard> encrypt(const uint& message,
                                                                   const Point& public_key) const;
                EncryptedMessage<EncryptionType::Standard> encrypt(const Point& message,
                                                                   const Point& public_key) const;
                EncryptedMessage<EncryptionType::Hashed>
                    encrypt(const uint& message,
                            const Point& public_key,
                            const std::function<uint(const Point&)>& hash_function) const;

                Point decrypt(const EncryptedMessage<EncryptionType::Standard>& encrypted_message,
                              const uint& private_key) const;
                uint decrypt(const EncryptedMessage<EncryptionType::Hashed>& encrypted_message,
                             const uint& private_key,
                             const std::function<uint(const Point&)>& hash_function) const;

            private:
                Point map_to_curve(const uint& message) const;
                uint map_to_uint(const Point& message) const;

                Curve m_curve;
                Point m_generator;
                uint m_generator_order;
            };
        }   // namespace encryption
    }       // namespace algorithm
}   // namespace elliptic_curve_guide

#endif
