/* 
 * File:   TestSourceType.h
 * Author: skiloop
 *
 * Created on 2014年1月15日, 上午10:55
 */

#ifndef TESTSOURCETYPE_H
#define	TESTSOURCETYPE_H

#include "sourceType.h"

class TestSourceType:public sourceType {
public:
    TestSourceType(int n=1);
    TestSourceType(const TestSourceType& orig);
    virtual ~TestSourceType();
    MyDataF valueAtTime(MyDataF t);
private:
    int mCount;
};

#endif	/* TESTSOURCETYPE_H */

