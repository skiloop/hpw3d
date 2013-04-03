/* 
 * File:   data1d.h
 * Author: skiloop
 *
 * Created on October 12, 2012, 3:08 PM
 */

#ifndef DATA1D_H
#define	DATA1D_H
#include <cstdlib>
template<typename T>
class data1d {
public:

    data1d(unsigned num = 0, T val = 0) : p(NULL), n(num) {
        createArray(num);
        initArray(val);
    };

    data1d(const data1d& orig) : p(NULL), n(orig.n) {
        if (orig.p != NULL && orig.n != 0) {
            p = new T[orig.n];
            for (unsigned i = 0; i < orig.n; i++)p[i] = orig.p[i];
        } else if (n != 0) {
            createArray(n, 0);
        }
    };

    ~data1d() {
        if (p != NULL)delete []p;
        p = NULL;
        n = 0;
    };

    void createArray(unsigned num) {
        if (num > 0) {
            p = new T[num];
            n = num;
        }
    };

    void initArray(T initval = 0) {
        if (p == NULL)return;
        for (unsigned i = 0; i < n; i++)p[i] = initval;
    };

    void resetArray() {
        initArray();
    };

    void createArray(unsigned num, T val) {
        createArray(num);
        initArray(val);
    };
public:
    T* p;
    unsigned n;
};
#endif	/* DATA1D_H */

