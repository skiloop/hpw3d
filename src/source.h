/* 
 * File:   source.h
 * Author: skiloop
 *
 * Created on 2013年12月16日, 下午12:10
 */

#ifndef SOURCE_H
#define	SOURCE_H

#include <vector>

#include "Point.h"
#include "data3d.h"
#include "sourceType.h"

class source {
public:

    source() {
    };
    source(int direction, MyDataF R, MyDataF amp);
    source(const source& orig);
    virtual ~source();
    /**
     * 
     * @param Cere
     * @param Cerhw
     * @param Cerhv
     * @param dx
     * @param dy
     * @param dz
     * @param dt
     */
    virtual void initCoefficients(data3d<MyDataF>&Cere, data3d<MyDataF>&Cerhw, data3d<MyDataF>&Cerhv,
            MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) = 0;

    /**
     * 
     * @param Cere
     * @param Cerhw
     * @param Cerhv
     * @param dx
     * @param dy
     * @param dz
     * @param dt
     */
    virtual void updateCoefficients(data3d<MyDataF>&Cere, data3d<MyDataF>&Cerhw, data3d<MyDataF>&Cerhv,
            MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) = 0;

    /**
     * 
     * @param Cere
     * @param Cerhw
     * @param Cerhv
     * @param Cervr
     * @param Beta
     * @param Ne
     * @param Nu_c
     * @param nes
     * @param bes
     * @param m
     * @param dx
     * @param dy
     * @param dz
     * @param dt
     */
    virtual void initCoefficients(data3d<MyDataF>&Cere, data3d<MyDataF>&Cerhw, data3d<MyDataF>&Cerhv, data3d<MyDataF>&Cervr,
            const data3d<MyDataF>&Beta, const data3d<MyDataF>&Ne, const data3d<MyDataF>&Nu_c,
            const Point& nes, const Point& bes, unsigned m,
            MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) = 0;

    /**
     * 
     * @param Cere
     * @param Cerhw
     * @param Cerhv
     * @param Cervr
     * @param Beta
     * @param Ne
     * @param Nu_c
     * @param nes
     * @param bes
     * @param m
     * @param dx
     * @param dy
     * @param dz
     * @param dt
     */
    virtual void updateCoefficients(data3d<MyDataF>&Cere, data3d<MyDataF>&Cerhw, data3d<MyDataF>&Cerhv, data3d<MyDataF>&Cervr,
            const data3d<MyDataF>&Beta, const data3d<MyDataF>&Ne, const data3d<MyDataF>&Nu_c,
            const Point& nes, const Point& bes, unsigned m,
            MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) = 0;

    /**
     * 
     * @param Cere
     * @param Cerhw
     * @param Cerhv
     * @param Cervr
     * @param Beta
     * @param Ne
     * @param vm
     * @param nes
     * @param bes
     * @param m
     * @param dx
     * @param dy
     * @param dz
     * @param dt
     */
    virtual void initCoefficients(data3d<MyDataF>&Cere, data3d<MyDataF>&Cerhw, data3d<MyDataF>&Cerhv, data3d<MyDataF>&Cervr,
            const data3d<MyDataF>&Beta, const data3d<MyDataF>&Ne, const MyDataF vm,
            const Point& nes, const Point& bes, unsigned m,
            MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) = 0;
    /**
     * 
     * @param Cere
     * @param Cerhw
     * @param Cerhv
     * @param Cervr
     * @param Beta
     * @param Ne
     * @param vm
     * @param nes
     * @param bes
     * @param m
     * @param dx
     * @param dy
     * @param dz
     * @param dt
     */
    virtual void updateCoefficients(data3d<MyDataF>&Cere, data3d<MyDataF>&Cerhw, data3d<MyDataF>&Cerhv, data3d<MyDataF>&Cervr,
            const data3d<MyDataF>&Beta, const data3d<MyDataF>&Ne, const MyDataF vm,
            const Point& nes, const Point& bes, unsigned m,
            MyDataF dx, MyDataF dy, MyDataF dz, MyDataF dt) = 0;



    /**
     * 
     * @param Ex
     * @param Ey
     * @param Ez
     * @param t
     */
    virtual void updateSource(data3d<MyDataF>&Ex, data3d<MyDataF>&Ey, data3d<MyDataF>&Ez, MyDataF t) = 0;

    /**
     * 
     * @param Ex
     * @param Ey
     * @param Ez
     * @param t
     */
    virtual void updateHardSource(data3d<MyDataF>&Ex, data3d<MyDataF>&Ey, data3d<MyDataF>&Ez, MyDataF t) = 0;

    const static int X = 1;
    const static int Y = 2;
    const static int Z = 3;

    /**
     * return direction of the source
     * @return mDirection
     */
    int getDirection() {
        return mDirection;
    };

    void setSourceType(sourceType*srcType) {
        mPSource = srcType;
    };

    void add(Point &point);

protected:
    int mDirection;
    MyDataF mResistancePerComponent;
    MyDataF mAmplitude;
    vector<Point> mSrcIndexes;
    vector<MyDataF>Cese;
    sourceType *mPSource;
};

#endif	/* SOURCE_H */

