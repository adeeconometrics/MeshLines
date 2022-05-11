/**
 * @file Test.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-03-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#pragma once
#include <iostream>

namespace tools {
    auto assert(bool b, const char* c_str) -> void {
        if (b == false){
            std::cout << c_str << '\n'; 
            exit(1);
        }
        return;
    }
}