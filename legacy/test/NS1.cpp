#include "NS1.h"
#include <iostream>

namespace ns { 
    auto test_function() -> void {
        std::cout << "test function is executed \n";
    }
}

int main(void) { 
    ns::test_function();
}