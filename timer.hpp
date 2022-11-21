#pragma once
#include <chrono>
#include <iostream>
class ScopeTimer {
public:
  ScopeTimer(std::string name)
      : m_name(std::move(name)),
        m_beg(std::chrono::high_resolution_clock::now()), m_tick(m_beg) {}
  void tick(const std::string &name) {
    auto end = std::chrono::high_resolution_clock::now();
    auto dur =
        std::chrono::duration_cast<std::chrono::microseconds>(end - m_tick);
    m_tick = std::chrono::high_resolution_clock::now();
    std::cout << ">>> " << m_name << " >>> " << name << ": " << dur.count()
              << " us\n";
  }
  ~ScopeTimer() {
    auto end = std::chrono::high_resolution_clock::now();
    auto dur =
        std::chrono::duration_cast<std::chrono::microseconds>(end - m_beg);
    std::cout << ">>> " << m_name << ": " << dur.count() << " us\n";
  }

private:
  std::string m_name;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_beg;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_tick;
};
