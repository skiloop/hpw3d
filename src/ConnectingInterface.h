/* 
 * File:   ConnectingInterface.h
 * Author: skiloop
 *
 * Created on 2013年12月17日, 下午3:39
 */

#ifndef CONNECTINGINTERFACE_H
#define	CONNECTINGINTERFACE_H

#include "data1d.h"
#include "data2d.h"
#include "data3d.h"
#include "Point.h"

class ConnectingInterface {
public:
    ConnectingInterface();
    ConnectingInterface(const ConnectingInterface& orig);
    virtual ~ConnectingInterface();
    void setLowerAndUpper(Point &lower, Point &upper);
    void setIncidentAngle(MyDataF theta, MyDataF phi, MyDataF psi, MyDataF lamda, MyDataF c, MyDataF dt, MyDataF ds);
    void initCoefficients(MyDataF ds, MyDataF dt);
    void updateESource(MyDataF source);
    void updateMSource();
    void updateEConnect(data3d<MyDataF>& Ex, const data3d<MyDataF> &Cexhy, const data3d<MyDataF> &Cexhz,
            data3d<MyDataF>& Ey, const data3d<MyDataF> &Ceyhx,const  data3d<MyDataF> &Ceyhz,
            data3d<MyDataF>& Ez, const data3d<MyDataF> &Cezhy, const data3d<MyDataF> &Cezhx);
    void updateMConnect(data3d<MyDataF> &Hx, data3d<MyDataF> &Chxey, const data3d<MyDataF> &Chxez,
            data3d<MyDataF>& Hy, const data3d<MyDataF> &Chyex, const data3d<MyDataF> &Chyez,
            data3d<MyDataF>& Hz, const data3d<MyDataF> &Chzey, const data3d<MyDataF> &Chzex);
    void invalidate();

    void saveEMInc(int step);
private:
    MyDataF mTheta;
    MyDataF mPhi;
    Point mLowerIndex;
    Point mUpperIndex;
    Point mOrigPoint;
    data1d<MyDataF> mEinc;
    data1d<MyDataF> mHinc;
    // amptidute of all components
    MyDataF mHx;
    MyDataF mHy;
    MyDataF mHz;
    MyDataF mEx;
    MyDataF mEy;
    MyDataF mEz;

    MyDataF mChe; // coefficient updating H field
    MyDataF mCeh; // coefficient updating E field
    MyDataF mCMur; // mur boundary condition coefficient
    MyDataF mPhaseVelocityRatio;

    //    // coefficients for  EM fields
    //    data3d<MyDataF> mCezhy;
    //    data3d<MyDataF> mCezhx;
    //    data3d<MyDataF> mCexhy;
    //    data3d<MyDataF> mCexhz;
    //    data3d<MyDataF> mCeyhz;
    //    data3d<MyDataF> mCeyhx;
    //    data3d<MyDataF> mChzey;
    //    data3d<MyDataF> mChzex;
    //    data3d<MyDataF> mChxey;
    //    data3d<MyDataF> mChxez;
    //    data3d<MyDataF> mChyez;
    //    data3d<MyDataF> mChyex;

    // delay in distance   
    data3d<MyDataF> mEzDelayFrontBack;
    data3d<MyDataF> mEzDelayLeftRight;

    data3d<MyDataF> mExDelayFrontBack;
    data3d<MyDataF> mExDelayTopBottom;

    data3d<MyDataF> mEyDelayLeftRight;
    data3d<MyDataF> mEyDelayTopBottom;

    data3d<MyDataF> mHzDelayFrontBack;
    data3d<MyDataF> mHzDelayLeftRight;

    data3d<MyDataF> mHxDelayFrontBack;
    data3d<MyDataF> mHxDelayTopBottom;

    data3d<MyDataF> mHyDelayLeftRight;
    data3d<MyDataF> mHyDelayTopBottom;

    MyDataF mDisC1;
    MyDataF mDisC2;
    MyDataF mDisC3;

    void initData();
    void setLowerIndex(Point &p);
    void setUpperIndex(Point &p);
    MyDataF distance(MyDataF i, MyDataF j, MyDataF k);
    void desideOrigPoint();
    MyDataF PhaseVelAngle(MyDataF phi,MyDataF theta, MyDataF lamba, MyDataF dt, MyDataF dx, MyDataF c);
};

#endif	/* CONNECTINGINTERFACE_H */

