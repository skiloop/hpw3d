/* 
 * File:   source.h
 * Author: skiloop
 *
 * Created on 2013年12月16日, 下午12:10
 */

#ifndef SOURCE_H
#define	SOURCE_H
#include "Point.h"
#include "data3d.h"
#include "sourceType.h"

class source {
public:

    source() {
    };
    source(int direction, MyDataF R, Point lower, Point upper, MyDataF amp);
    source(const source& orig);
    virtual ~source();
    /**
     * 
     * @param Cere
     * @param Cerhw
     * @param Cerhv
     * @param dr
     * @param dw
     * @param dv
     * @param dt
     */
    virtual void initUpdateCoefficients(data3d<MyDataF>&Cere, data3d<MyDataF>&Cerhw, data3d<MyDataF>&Cerhv,
            MyDataF dr, MyDataF dw, MyDataF dv, MyDataF dt) = 0;

    /**
     * 
     * @param Ex
     * @param Ey
     * @param Ez
     * @param t
     */
    virtual void updateSource(data3d<MyDataF>&Ex, data3d<MyDataF>&Ey, data3d<MyDataF>&Ez, MyDataF t) = 0;

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
protected:
    int mDirection;
    MyDataF mResistancePerComponent;
    MyDataF mAmplitude;
    Point mLowerIndex;
    Point mUpperIndex;
    data3d<MyDataF> Cese;
    sourceType *mPSource;
};

#endif	/* SOURCE_H */

