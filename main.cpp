#include "Uint/uint.h"

#include <iostream>

int main() {
    ECG::uint_t<512> a = nullptr;
    ECG::uint_t<512> b = std::string("9");
    std::cout << (23 * b).into_string();
}
