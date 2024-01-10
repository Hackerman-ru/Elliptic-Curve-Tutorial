#include "LongInt/long-int.h"

#include <iostream>

int main() {
    ECG::uint_t<512> a("a123ff1fffffffff123", ECG::StringType::HEXADECIMAL);
    ECG::uint_t<512> b("b00af100fffffff", ECG::StringType::HEXADECIMAL);
    ECG::uint_t<512> c = a % b;
    std::cout << a.into_string(ECG ::StringType::HEXADECIMAL) << '\n';
    std::cout << b.into_string(ECG::StringType::HEXADECIMAL) << '\n';
    std::cout << c.into_string(ECG::StringType::HEXADECIMAL) << '\n';
}