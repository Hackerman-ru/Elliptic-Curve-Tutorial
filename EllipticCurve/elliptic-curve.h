#ifndef ECG_ELLIPTIC_CURVE_H
#define ECG_ELLIPTIC_CURVE_H

#include "../Field/Field.h"

namespace ECG {
    class ECE {
    public:
        using Field = PFE;
        ECE() = default;
        explicit ECE(const Field& x);
        template<is_convertible<Field> T, is_convertible<Field> W>
        ECE(const T& x, const W& y) : m_x(Field(x)), m_y(Field(y)) {};
        static void set_params(const Field& a, const Field& b);

        ECE operator+(const ECE& other) const;
        ECE operator-(const ECE& other) const;
        ECE operator-() const;

        ECE& operator+=(const ECE& other);
        ECE& operator-=(const ECE& other);

        bool operator==(const ECE& other) const;

        bool is_inf() const;

        static ECE find_point(const std::string& str);

    private:
        Field m_x = Field(0);
        Field m_y = Field(0);

        static Field m_a;
        static Field m_b;

        static bool find_y(const Field& x, Field* y);
        static Field::uint generate_random_uint();
        static Field find_n();
    };

    ECE::Field ECE::m_a(0);
    ECE::Field ECE::m_b(0);
}   // namespace ECG

#endif
