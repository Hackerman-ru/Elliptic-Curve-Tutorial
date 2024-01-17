#ifndef ECG_LONGINT_H
#define ECG_LONGINT_H

#ifndef CHAR_BIT
    #define CHAR_BIT 8
#endif

#include "../util.h"

#include <array>
#include <cassert>
#include <complex>
#include <functional>
#include <string>

namespace ECG {
    template<size_t bits>
    class uint_t {
    public:
        constexpr uint_t() = default;

        template<is_convertible<uint32_t> T>
        constexpr explicit uint_t(const T& value) : m_buckets({static_cast<uint32_t>(value)}) {};

        template<ConvertibleContainer<uint32_t> T>
        uint_t(const T& other);

        uint_t(const std::string& str, StringType str_type = StringType::DECIMAL);
        uint_t(const std::string& str, std::function<uint32_t(char)> map, size_t shift);

        auto operator<=>(const uint_t& other) const = default;
        bool operator==(const uint_t& other) const;

        template<is_explicitly_convertible<uint_t> T>
        bool operator==(const T& other) const {
            return *this == uint_t(other);
        }

        uint_t operator+(const uint_t& other) const;
        uint_t operator-(const uint_t& other) const;
        uint_t operator*(const uint_t& other) const;
        uint_t operator*(uint32_t other) const;
        uint_t operator/(const uint_t& other) const;
        uint_t operator%(const uint_t& other) const;
        uint_t operator/(uint32_t other) const;
        uint_t operator%(uint32_t other) const;
        uint_t operator>>(size_t shift) const;
        uint_t operator<<(size_t shift) const;
        uint_t operator^(const uint_t& other) const;
        uint_t operator|(const uint_t& other) const;
        uint_t operator&(const uint_t& other) const;

        uint_t operator-() const;

        uint_t& operator+=(const uint_t& other);
        uint_t& operator-=(const uint_t& other);
        uint_t& operator*=(const uint_t& other);
        uint_t& operator*=(uint32_t other);
        uint_t& operator/=(const uint_t& other);
        uint_t& operator%=(const uint_t& other);
        uint_t& operator/=(uint32_t other);
        uint_t& operator%=(uint32_t other);
        uint_t& operator>>=(size_t shift);
        uint_t& operator<<=(size_t shift);
        uint_t& operator^=(const uint_t& other);
        uint_t& operator|=(const uint_t& other);
        uint_t& operator&=(const uint_t& other);

        uint_t& operator++(int);
        uint_t& operator++();
        uint_t& operator--();

        template<typename T>
        requires is_convertible<uint32_t, T> || std::is_same_v<T, std::string>
        T convert_to() const {
            size_t shift = sizeof(T) * CHAR_BIT;
            size_t bucket_number = shift / c_BUCKET_SIZE;
            T result = 0;

            for (size_t i = 0; i < c_BUCKET_NUMBER && i < bucket_number; ++i) {
                result |= T(m_buckets[i]) << (i * c_BUCKET_SIZE);
            }

            shift %= c_BUCKET_SIZE;

            result |= ((T(m_buckets[bucket_number]) << (c_BUCKET_SIZE - shift)) >> (c_BUCKET_SIZE - shift))
                   << (bucket_number * c_BUCKET_SIZE);

            return result;
        }

        std::string into_string(StringType str_type = StringType::DECIMAL) const;
        std::string into_string(std::function<char(uint32_t)> map, size_t shift) const;

        const uint32_t& operator[](size_t pos) const;
        uint32_t& operator[](size_t pos);

        static constexpr size_t size();
        size_t clz() const;

    private:
        static uint_t divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr);
        static uint_t divide(const uint_t& lhs, const uint32_t& rhs, uint_t* remainder = nullptr);
        static uint_t d_divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr);

        static constexpr size_t c_BUCKET_SIZE = sizeof(uint32_t) * CHAR_BIT;
        static constexpr size_t c_BUCKET_NUMBER = bits / c_BUCKET_SIZE;

        std::array<uint32_t, c_BUCKET_NUMBER> m_buckets = {};
    };
}   // namespace ECG

#include "uint.tpp"

#endif
