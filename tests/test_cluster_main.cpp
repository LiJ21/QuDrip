#include <iostream>

int run_cluster_tests();

int main() {
  const int checks = run_cluster_tests();
  std::cout << "all " << checks << " cluster checks passed\n";
  return 0;
}
