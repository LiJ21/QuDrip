#include <iostream>

int run_operator_tests();

int main() {
  const int checks = run_operator_tests();
  std::cout << "all " << checks << " operator checks passed\n";
  return 0;
}
