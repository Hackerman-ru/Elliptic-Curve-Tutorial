#include "Uint/uint.h"

#include <iostream>

using namespace ECG;

int main() {
    uint_t<512> a = 5;
    a *= a;
    std::cout << a.into_string();
    return 0;
}
