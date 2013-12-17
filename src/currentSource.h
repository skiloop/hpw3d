/* 
 * File:   currentSource.h
 * Author: skiloop
 *
 * Created on 2013年12月16日, 下午12:05
 */

#ifndef CURRENTSOURCE_H
#define	CURRENTSOURCE_H
#include "source.h"

class currentSource : public source {
public:

    currentSource() {
    };
    currentSource(int direction, MyDataF R, Point lower, Point upper, MyDataF amp);
    currentSource(const currentSource& orig);
    virtual ~currentSource();
    void initUpdateCoefficients(data3d<MyDataF>&Ex, data3d<MyDataF>&Ey, data3d<MyDataF>&Ez,
            MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt);
    void updateSource(data3d<MyDataF>&Ex, data3d<MyDataF>&Ey, data3d<MyDataF>&Ez, MyDataF t);
private:

};

#endif	/* CURRENTSOURCE_H */

