#include <random>
#include <iostream>

int main() {

    std::uniform_int_distribution<> dist(1, 100); // Distribution range

    // Generate and print 10 random numbers
    for (int i = 0; i < 10; ++i) {
    std::random_device rd;  // Seed for random number generator
    std::mt19937 g(rd());  // Standard Mersenne Twister engine
        std::cout << g << " ";
    }
    std::cout << std::endl;

    return 0;
}

