#ifndef ECG_STRING_PARSER_H
#define ECG_STRING_PARSER_H

#include "util.h"

#include <array>
#include <limits>
#include <string>

namespace ECG {
    enum class StringType {
        BINARY = 0b1,
        DECIMAL = 10,
        HEXADECIMAL = 0xF,
    };

    template<typename Block, size_t Size, typename T>
    constexpr std::array<Block, Size> split_into_blocks(const T& value);

    template<typename Block, size_t Size>
    constexpr std::array<Block, Size> split_into_blocks<Block, Size, char*>(const char*& str) {
        if (str == nullptr) {
            return {};
        }

        static constexpr bits = sizeof(Block) * Size;
        uint_t<bits> result;
        const char* pos = str;

        while (*str != '\0') {
            uint_t<bits> temp = (result <<= 1);
            result <<= 2;
            result += temp;
            result += uint_t<bits>(*str - '0');
            ++str;
        }

        return result.m_blocks;
    }

    template<typename Block, size_t Size>
    constexpr std::array<Block, Size> split_into_blocks<Block, Size, std::string>(const std::string& value);

    template<typename Block, size_t Size, typename T>
    requires std::numeric_limits<T>::is_integer && is_upcastable_to<T, Block>
    constexpr std::array<Block, Size> split_into_blocks(T value);

    template<typename Block, size_t Size, typename T>
    requires std::numeric_limits<T>::is_integer && is_downcastable_to<T, Block>
    static constexpr std::array<Block, Size> split_into_blocks(T value) {}

    static StringType find_type(const std::string& str) {
        if (str.size() < 2) {
            return StringType::DECIMAL;
        }

        switch (str[1]) {
        case 'x' :
            return StringType::HEXADECIMAL;
        case 'b' :
            return StringType::BINARY;
        default :
            return StringType::DECIMAL;
        }
    }
}   // namespace ECG
#endif
