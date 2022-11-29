#!/bin/bash

[ -d ./build ] ||  mkdir build
g++ -o build/sparsify main.cpp -std=c++17 -lpthread -Wall -Wno-unused-result -Wno-sign-compare -Ofast -g -mavx2 -ltbb -flto
