/* 
 * File:   data1d.h
 * Author: skiloop
 *
 * Created on October 12, 2012, 3:08 PM
 */

#ifndef DATA1D_H
#define	DATA1D_H
template<typename T>
class data1d {
public:
    data1d();
    data1d(const data1d& orig);
    virtual ~data1d();
private:
    T* p;
    unsigned n;
};

#endif	/* DATA1D_H */

