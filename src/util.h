#ifndef ECG_UTIL_H
#define ECG_UTIL_H

#include <concepts>

namespace ECG {
    template<typename From, typename To>
    concept is_explicitly_convertible = !std::is_same_v<From, To> && requires(From f) { static_cast<To>(f); };

    template<typename From, typename To>
    concept is_convertible_to = requires(From f) { static_cast<To>(f); };

    template<typename From, typename To>
    concept is_upcastable_to = requires(From f) {
        static_cast<To>(f);
        sizeof(From) <= sizeof(To);
    };

    template<typename From, typename To>
    concept is_downcastable_to = requires(From f) {
        static_cast<To>(f);
        sizeof(From) > sizeof(To);
    };

    template<typename To, typename From>
    concept is_convertible_from = is_convertible_to<From, To>;

    template<typename FromContainer, typename ToType>
    concept ConvertibleContainer = requires(FromContainer t, size_t i) {
        { t[i] } -> is_convertible_to<ToType>;
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
