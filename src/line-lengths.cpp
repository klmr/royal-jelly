#include <cstdlib>
#include <iostream>

auto main() -> int {
    std::ios_base::sync_with_stdio(false);

    auto chr = char{};
    auto line_length = int{};

    while (std::cin.get(chr)) {
        if (chr == '\n') {
            std::cout << line_length << '\n';
            line_length = 0;
        } else {
            ++line_length;
        }
    }

    if (std::cin.bad()) {
        std::cerr << "Error reading input.\n";
        return EXIT_FAILURE;
    }
}
