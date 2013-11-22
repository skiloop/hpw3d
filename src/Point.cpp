/* 
 * File:   Point.cpp
 * Author: skiloop
 * 
 * Created on 2013年11月22日, 下午4:43
 */

#include "Point.h"

Point::Point()
: x(0), y(0), z(0) {
}

Point::Point(const Point& orig)
: x(orig.x)
, y(orig.y)
, z(orig.z) {

}

Point::~Point() {
}

bool Point::checkMax(const Point&maxPoint) const {
    return (x < maxPoint.x && y < maxPoint.y && z < maxPoint.z);
}

bool Point::checkMax(unsigned xMax, unsigned yMax, unsigned zMax) const {
    return (x < xMax && y < yMax && z < zMax);
}

