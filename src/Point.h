/* 
 * File:   Point.h
 * Author: skiloop
 *
 * Created on 2013年11月22日, 下午4:43
 */

#ifndef POINT_H
#define	POINT_H

class Point {
public:
    Point();
    Point(const Point& orig);
    virtual ~Point();
    unsigned int x;
    unsigned int y;
    unsigned int z;
    
    bool  checkMax(const Point&maxPoint) const;
        bool  checkMax(unsigned xMax,unsigned yMax,unsigned zMax) const;

private:
    

};

#endif	/* POINT_H */

