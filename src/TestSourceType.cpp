/* 
 * File:   TestSourceType.cpp
 * Author: skiloop
 * 
 * Created on 2014年1月15日, 上午10:55
 */

#include "TestSourceType.h"

TestSourceType::TestSourceType(int n)
: mCount(n) {
}

TestSourceType::TestSourceType(const TestSourceType& orig)
: mCount(orig.mCount) {
}

TestSourceType::~TestSourceType() {
}

MyDataF TestSourceType::valueAtTime(MyDataF t) {
    if (mCount > 0) {
        mCount--;
        return 1.0;
    } else {
        return 0.0;
    }
}
