#ifndef ECG_UTIL_H
#define ECG_UTIL_H

#include <concepts>

namespace ECG {
    template<typename From, typename To>
    concept is_convertible_to = requires(From f) {
        { static_cast<To>(f) } noexcept;
    };

    template<typename From, typename To>
    concept is_upcastable_to = sizeof(From) <= sizeof(To) && is_convertible_to<From, To>;

    template<typename From, typename To>
    concept is_downcastable_to = sizeof(From) > sizeof(To) && is_convertible_to<From, To>;

    template<typename Container, typename T>
    concept is_convertible_container = requires(Container t, size_t i) {
        { t[i] } -> is_convertible_to<T>;
        { t.size() } -> std::same_as<size_t>;
    };
}   // namespace ECG
#endif
