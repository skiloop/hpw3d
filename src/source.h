/* 
 * File:   source.h
 * Author: skiloop
 *
 * Created on 2013年5月10日, 下午9:58
 */

#ifndef SOURCE_H
#define	SOURCE_H

class Source {
public:
    Source();
    Source(const Source& orig);
    virtual ~Source();

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
private:

};

#endif	/* SOURCE_H */

