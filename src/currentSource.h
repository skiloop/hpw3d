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
    
    void updateHardSource(data3d<MyDataF>&Ex, data3d<MyDataF>&Ey, data3d<MyDataF>&Ez, MyDataF t);

    void initUpdateCoefficients(data3d<MyDataF>&Cere, data3d<MyDataF>&Cerhw, data3d<MyDataF>&Cerhv, data3d<MyDataF>&Cervr,
            const data3d<MyDataF>&Beta, const data3d<MyDataF>&Ne, const MyDataF vm,
            const Point& nes, const Point& bes, unsigned m,
            MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt);

    virtual void initUpdateCoefficients(data3d<MyDataF>&Cere, data3d<MyDataF>&Cerhw, data3d<MyDataF>&Cerhv, data3d<MyDataF>&Cervr,
            const data3d<MyDataF>&Beta, const data3d<MyDataF>&Ne, const data3d<MyDataF>&Nu_c,
            const Point& nes, const Point& bes, unsigned m,
            MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt);
private:
    void desideXYZ(MyDataF &dr, MyDataF &dw, MyDataF &dv,const MyDataF dx, const MyDataF dy, const MyDataF dz);
};

#endif	/* CURRENTSOURCE_H */

