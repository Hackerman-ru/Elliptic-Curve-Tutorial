#ifndef ECG_UTIL_H
#define ECG_UTIL_H

#include <concepts>

namespace ECG {
    template<typename From, typename To>
    concept is_explicitly_convertible = !std::is_same_v<From, To> && requires(From f) { static_cast<To>(f); };

    template<typename From, typename To>
    concept is_convertible_to = requires(From f) { static_cast<To>(f); };

    template<typename To, typename From>
    concept is_convertible_from = is_convertible_to<From, To>;

    template<typename T, typename W>
    concept ConvertibleContainer = requires(T t, size_t i) {
        { t[i] } -> is_convertible_to<W>;
        { t.size() } -> std::same_as<size_t>;
    };

    enum class StringType {
        BINARY = 0b1,
        DECIMAL = 10,
        HEXADECIMAL = 0xF,
    };

    enum class CoordinatesType {
        Normal,
        Projective,
        Jacobi,
        ModifiedJacobi,
        JacobiChudnovski,
        SimplifiedJacobiChudnovski,
    };
}   // namespace ECG
#endif
