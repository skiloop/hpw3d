/* 
 * File:   SineWaveSource.h
 * Author: skiloop
 *
 * Created on 2013年12月16日, 下午3:53
 */

#ifndef SINEWAVESOURCE_H
#define	SINEWAVESOURCE_H
#include "sourceType.h"

class SineWaveSource : public sourceType {
public:
    SineWaveSource(MyDataF omega);
    SineWaveSource(const SineWaveSource& orig);
    virtual ~SineWaveSource();
    MyDataF valueAtTime(MyDataF t);
private:
    MyDataF mOmega;
};

#endif	/* SINEWAVESOURCE_H */

