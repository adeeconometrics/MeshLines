/**
 * @file Vector.cpp
 * @author ddamiana
 * @brief representation of a mathematical vector
 * @version 0.1
 * @date 2021-12-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <initializer_list>

using std::vector;

template <typename T> class Vector{
    private:
        vector<T> m_vector;
    
    public:
        Vector(const vector<T>& _vector): m_vector(_vector) {}
        Vector(std::initializer_list<T> _list): m_vector ((vector<T>)_list) {}

        // Vector operator=(const Vector &rhs) const;
        // Vector operator=(Vector &&rhs) const;

        Vector operator*(const Vector &rhs) const{
          if (size() != rhs.size()) {
            throw std::domain_error(
                "size of vector is expected to have the same size");
          }

          vector<T> result(rhs.size());
          std::transform(m_vector.cbegin, m_vector.cend(),
                         rhs.m_vector.cbegin(), result.begin(),
                         [](const T &a, const T &b) -> T { return a * b; });

          return result;
        }

        Vector operator*(float scalar) const noexcept {
            vector<T> result (m_vector.size());
            std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                            [&scalar] (const T& a) -> T {return a * scalar;});
            return result;
        }

        Vector operator+(const Vector &rhs) const {
          if (size() != rhs.size()) {
            throw std::domain_error(
                "size of vector is expected to have the same size");
          }

          vector<T> result(rhs.size());
          std::transform(m_vector.cbegin, m_vector.cend(),
                         rhs.m_vector.cbegin(), result.begin(),
                         [](const T &a, const T &b) -> T { return a + b; });

          return result;
        }

        Vector operator+(float scalar) const noexcept{
          vector<T> result(m_vector.size());
          std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                         [&scalar](const T &a) -> T { return a + scalar; });
          return result;
        }

        Vector operator-(const Vector &rhs) const{
          if (size() != rhs.size()) {
            throw std::domain_error(
                "size of vector is expected to have the same size");
          }

          vector<T> result(rhs.size());
          std::transform(m_vector.cbegin, m_vector.cend(),
                         rhs.m_vector.cbegin(), result.begin(),
                         [](const T &a, const T &b) -> T { return a - b; });

          return result;
        }

        Vector operator-(float scalar) const noexcept{
          vector<T> result(m_vector.size());
          std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                         [&scalar](const T &a) -> T { return a - scalar; });
          return result;
        }

        Vector operator/(const Vector &rhs) const{
          if (size() != rhs.size()) {
            throw std::domain_error(
                "size of vector is expected to have the same size");
          }

          vector<T> result(rhs.size());
          std::transform(m_vector.cbegin, m_vector.cend(),
                         rhs.m_vector.cbegin(), result.begin(),
                         [](const T &a, const T &b) -> T { return a / b; });

          return result;
        }

        Vector operator/(float scalar) const noexcept{
          vector<T> result(m_vector.size());
          std::transform(m_vector.cbegin(), m_vector.cend(), result.begin(),
                         [&scalar](const T &a) -> T { return a / scalar; });
          return result;
        }

        bool operator==(const Vector &rhs) const noexcept{
          for (size_t i = 0; i < m_vector.size(); i++) {
            if (m_vector[i] != rhs.m_vector[i])
              return false;
          }
          return true;
        }

        bool operator!=(const Vector &rhs) const noexcept{
            return !(*this == rhs);
        }

        unsigned int size() const{
            return m_vector.size();
        }

        void print () {
            std::cout << '[';
            for(auto& i: m_vector){
                std::cout  << i << ' ';
            }
            std::cout << ']';
        }
};

int main(){
  vector<int> a{1, 2, 3};
  vector<int> b{1, 1, 1};
  vector<int> c{1, 1, 1};
  bool bc = b == c, ba = b == a;
  std::cout << bc << '\n';
  std::cout << ba;
}