#include"Matrix.hpp"

using namespace me;

int main()
{
    matrix<int> a(1, 2);
    matrix<int> b(1, 2,{1,2});
    std::cout << a << b;
    a = std::move(b);
    std::cout << a<<b;
}