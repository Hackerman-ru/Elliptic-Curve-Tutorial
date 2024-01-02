#pragma once

namespace ECG {
    using _Ty = int;

    template<_Ty p>
    class PrimeField {
        public:
        PrimeField(const _Ty& value) : m_value(value) {};
        PrimeField operator+(const PrimeField& other) const;
        PrimeField operator-(const PrimeField& other) const;
        PrimeField operator-() const;
        PrimeField operator*(const PrimeField& other) const;
        PrimeField operator/(const PrimeField& other) const;
        PrimeField inverse() const;
        PrimeField& operator+=(const PrimeField& other);
        PrimeField& operator-=(const PrimeField& other);
        PrimeField& operator*=(const PrimeField& other);
        PrimeField& operator/=(const PrimeField& other);
        PrimeField& operator/=(const PrimeField& other);

        private:
        _Ty m_value;
    };

    template<_Ty p, _Ty n>
    class Field {};
}   // namespace ECG
