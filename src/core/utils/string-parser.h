#ifndef ECG_STRING_PARSER_H
#define ECG_STRING_PARSER_H

#include "concepts.h"

#include <cassert>

namespace elliptic_curve_guide {
    namespace algorithm {
        // Parsing invalid string of unsigned integer is UB
        template<typename T>
        requires concepts::is_integral<T>
        T parse_into_uint(const char* str) {
            assert(str != nullptr && "parse_into_uint got nullptr");

            T value = 0;
            uint16_t radix = 10;

            if (str[0] == '0' && str[1] == 'x') {
                radix = 16;
                str += 2;
            } else if (str[0] == '0' && str[1] == 'b') {
                radix = 2;
                str += 2;
            } else if (str[0] == '0') {
                radix = 8;
                ++str;
            }

            while (*str != '\0') {
                switch (radix) {
                case 16 :
                    value <<= 4;
                    break;
                case 8 :
                    value <<= 3;
                    break;
                case 2 :
                    value <<= 1;
                    break;
                default :
                    value *= static_cast<T>(radix);
                    break;
                }

                uint16_t symbol_value = radix + 1;

                if (*str >= '0' && *str <= '9') {
                    symbol_value = static_cast<uint16_t>(*str - '0');
                } else if (*str >= 'a' && *str <= 'f') {
                    symbol_value = static_cast<uint16_t>(*str - 'a') + 10;
                } else if (*str >= 'A' && *str <= 'F') {
                    symbol_value = static_cast<uint16_t>(*str - 'A') + 10;
                }

                if (symbol_value >= radix) {
                    assert(false && "parse_into_uint got incorrect string");
                }

                value += static_cast<T>(symbol_value);
                ++str;
            }

            return value;
        }
    }   // namespace algorithm
}   // namespace elliptic_curve_guide
#endif
