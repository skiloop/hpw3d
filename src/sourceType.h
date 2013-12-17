/* 
 * File:   source.h
 * Author: skiloop
 *
 * Created on 2013年5月10日, 下午9:58
 */

#ifndef SOURCE_TYPE_H
#define	SOURCE_TYPE_H

#include "common.h"

class sourceType {
public:
    sourceType();
    sourceType(const sourceType& orig);
    virtual ~sourceType();

    /**
     * Sine Pulse
     * @param t
     * @param omega
     * @param t_max
     * @param t_min
     * @return 
     */
    static MyDataF SinePulse(MyDataF t, MyDataF omega, MyDataF t_max, MyDataF t_min);

    /**
     * Sine Wave
     * @param t
     * @param omega
     * @return 
     */
    static MyDataF SineWave(MyDataF t, MyDataF omega);

    /**
     * 
     * @param t
     * @param tw
     * @return 
     */
    static MyDataF GaussianPulse(MyDataF t, MyDataF tw);
    
    /**
     * get value at time @c t
     * @param t
     * @return 
     */
    virtual MyDataF valueAtTime(MyDataF t);
private:

};

#endif	/* SOURCE_TYPE_H */

