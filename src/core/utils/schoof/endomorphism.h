#ifndef ECG_ENDOMORPHISM_H
#define ECG_ENDOMORPHISM_H

#include "ring.h"
#include "utils/wnaf.h"

#include <variant>

namespace elliptic_curve_guide {
    namespace endomorphism {
        class End {
            using Poly = polynomial::Poly;
            using Ring = ring::Ring;
            using Field = field::Field;
            using RingElement = ring::RingElement;
            using FieldElement = field::FieldElement;

            friend End algorithm::wnaf_addition<End>(End end, const uint& value);

        public:
            struct Info {
                Ring ring;
                FieldElement a;
                RingElement curve_function;
            };

            End(const Poly& a, const Poly& b, std::shared_ptr<const Info> info);
            End(const RingElement& a, const RingElement& b, std::shared_ptr<const Info> info);
            End(RingElement&& a, RingElement&& b, std::shared_ptr<const Info> info);

            using AdditionResult = std::variant<End, Poly>;

            static AdditionResult twice(const End& end);

            friend AdditionResult operator+(const End& lhs, const End& rhs);

            friend AdditionResult operator-(const End& lhs, const End& rhs);

            friend End operator*(const End& lhs, const End& rhs);
            friend End operator*(End&& lhs, const End& rhs);
            friend End operator*(const End& lhs, End&& rhs);
            friend End operator*(End&& lhs, End&& rhs);

            friend AdditionResult operator*(const End& end, const uint& value);
            friend AdditionResult operator*(const uint& value, const End& end);

            bool operator==(const End& other) const;

            End operator-() const;

            End& operator*=(const End& other);

        private:
            friend static AdditionResult multiply(End value, const uint& n);

            void nullify();

            RingElement m_a_x;
            RingElement m_b_x;
            std::shared_ptr<const Info> m_info;
        };
    }   // namespace endomorphism
}   // namespace elliptic_curve_guide
#endif
