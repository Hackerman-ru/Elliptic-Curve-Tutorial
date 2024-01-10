#ifndef ECG_LONGINT_H
#define ECG_LONGINT_H

#ifndef CHAR_BIT
    #define CHAR_BIT 8
#endif

#include <array>
#include <cassert>
#include <complex>
#include <functional>
#include <string>

template<typename From, typename To>
concept Convertible = std::is_nothrow_convertible_v<From, To>;

namespace ECG {
    template<typename T, typename W>
    concept ConvertableContainer = requires(T t, size_t i) {
        { t[i] } -> std::convertible_to<W>;
        { t.size() } -> std::same_as<size_t>;
    };

    template<size_t bits = 512>
    class uint_t {
        using bucket_type = uint32_t;
        static_assert(std::is_same_v<bucket_type, uint32_t>, "The bucket type must be uint32_t");

    public:
        enum class StringType {
            BINARY = 0b1,
            DECIMAL = 10,
            HEXADECIMAL = 0xF,
        };

        constexpr uint_t() = default;
        constexpr explicit uint_t(const bucket_type& value);

        template<ConvertableContainer<bucket_type> T>
        uint_t(const T& other) {
            const size_t min_size = std::min(size(), other.size());

            for (size_t i = 0; i < min_size; i++) {
                m_buckets[i] = other[i];
            }
        }

        uint_t(const std::string& str, StringType str_type);
        uint_t(const std::string& str, std::function<bucket_type(char)> map, size_t shift);

        auto operator<=>(const uint_t& other) const = default;

        uint_t operator+(const uint_t& other) const;
        uint_t operator-(const uint_t& other) const;
        uint_t operator*(const uint_t& other) const;
        uint_t operator*(bucket_type other) const;
        uint_t operator/(const uint_t& other) const;
        uint_t operator%(const uint_t& other) const;
        uint_t operator/(bucket_type other) const;
        uint_t operator%(bucket_type other) const;
        uint_t operator>>(size_t shift) const;
        uint_t operator<<(size_t shift) const;
        uint_t operator^(const uint_t& other) const;
        uint_t operator|(const uint_t& other) const;
        uint_t operator&(const uint_t& other) const;

        uint_t operator-() const;

        uint_t& operator+=(const uint_t& other);
        uint_t& operator-=(const uint_t& other);
        uint_t& operator*=(const uint_t& other);
        uint_t& operator*=(bucket_type other);
        uint_t& operator/=(const uint_t& other);
        uint_t& operator%=(const uint_t& other);
        uint_t& operator/=(bucket_type other);
        uint_t& operator%=(bucket_type other);
        uint_t& operator>>=(size_t shift);
        uint_t& operator<<=(size_t shift);
        uint_t& operator^=(const uint_t& other);
        uint_t& operator|=(const uint_t& other);
        uint_t& operator&=(const uint_t& other);

        uint_t& operator++(int);
        uint_t& operator++();
        uint_t& operator--();

        explicit operator bucket_type() const;

        template<typename T>
        requires(std::is_nothrow_convertible_v<bucket_type, T>)
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

        std::string into_string(StringType str_type) const;
        std::string into_string(std::function<char(bucket_type)> map, size_t shift) const;

        const bucket_type& operator[](size_t pos) const {
            return m_buckets[pos];
        }

        bucket_type& operator[](size_t pos) {
            return m_buckets[pos];
        }

        constexpr size_t size() const;
        size_t clz() const;

    private:
        static constexpr size_t c_BUCKET_SIZE = sizeof(bucket_type) * CHAR_BIT;
        static constexpr size_t c_BUCKET_NUMBER = bits / c_BUCKET_SIZE;

        std::array<bucket_type, c_BUCKET_NUMBER> m_buckets = {};

        static uint_t divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr);
        static uint_t divide(const uint_t& lhs, const bucket_type& rhs, uint_t* remainder = nullptr);
        static uint_t d_divide(const uint_t& lhs, const uint_t& rhs, uint_t* remainder = nullptr);
    };
}   // namespace ECG

#include "long-int.tpp"

#endif
