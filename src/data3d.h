#ifndef DATA3D_H
#define DATA3D_H


#ifdef MATLAB_SIMULATION

#include <engine.h>
#include <mex.h>
#ifdef printf
#undef printf
#endif

#endif // end MATLAB_SIMULATION
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;

//#include "microdef.h"
#define MAX_ARRAY_SIZE 300000
#include "common.h"
#include "Point.h"

/**
 * nx:point count in x direction
 * ny:point count in x direction
 * p:pointer to store data;
 */
template<class DataType>
class data3d {
public:
    unsigned int nx;
    unsigned int ny;
    unsigned int nz;
    DataType*** p;

public:
    static const std::string OUTPUT_FILE_NAME_TAIL;

#ifdef MATLAB_SIMULATION
    static unsigned int mMatlabFigureCount;
    static Engine *ep;
private:
    static bool mIsMatlabEngineStarted;
#endif
private:
    std::string mName; // name for this to save file

private:
#ifdef MATLAB_SIMULATION
    int mMatlabFigureIndex;
    mxArray *num;
    mxArray *mMatlabMXArray;
#endif
public:

    /**
     * if cx==0 then p = NULL
     * else p has cx pointers;
     *
     * if cy == 0 then all cx pointers of p is NULL
     * else p[i] has cy pointers with i from 0 to cx-1;
     *
     * when cannot create space for p and p[i],exit program;
     */
    data3d(unsigned int cx, unsigned int cy, unsigned cz)
    : nx(cx), ny(cy), nz(cz), p(NULL)
#ifdef MATLAB_SIMULATION
    , mMatlabFigureIndex(-1)
#endif
    {
        unsigned i, j;
        if (cx == 0 || cy == 0 || cz == 0) {
            return;
        }
        try {

            p = new DataType**[cx];
            for (i = 0; i < cx; i++) {
                p[i] = new DataType*[cy];
            }
            for (i = 0; i < cx; i++) {
                for (j = 0; j < cy; j++) {
                    p[i][j] = new DataType[cz];
                }
            }
        } catch (exception & e) {
            cerr << e.what() << endl;
            return;
        }

    };

    /**
     * default constructor set
     *  @c nx = 0
     *  @c ny = 0
     *  @c nz = 0
     *  @c p=NULL
     */
    data3d() : nx(0), ny(0), nz(0), p(NULL)
#ifdef MATLAB_SIMULATION
    , mMatlabFigureIndex(-1)
#endif
    {
    };

    /**
     * copy constructor
     * @param obj
     */
    data3d(const data3d< DataType > &obj);

    /**
     * deconstructor
     */
    ~data3d();

    /**
     * print p in struct data3d
     */
    void printArray();

    /**
     * free space created for data3d @c mst
     */
    void freeArray();

    /**
     * Set Data to val
     */
    int resetArray(DataType val = 0);

    /**
     * check p of data3d @c mst is valid
     * if p is not NULL and none of its subpointers,then
     * return true,otherwise false
     */
    bool checkArray();

    /**
     * Create Space for struct data3d and initialize its @c nx and @c ny
     */
    int create3DArray(unsigned nnx, unsigned nny, unsigned nnz);

    /**
     * Create Space for struct data3d and initialize its @c nx and @c ny
     */
    int create3DArray(unsigned nnx, unsigned nny, unsigned nnz, DataType initVal);

    /**
     * Copy all p in st to stpre
     * Dimensions of @c st and that of @c pstruct must macth,and both with valid
     * p
     */
    int backup3DArray(const data3d< DataType > &mstru);

    /**
     * @brief Save p of  data3d data skipping p rows and p columns
     *
     */
    void saveArrayData(const unsigned num, unsigned leap = 0);

    /**
     *
     * @param other
     */
    void operator=(data3d< DataType > const &other);

    /**
     *  return this.p[index.x][index.y][index.z]
     *
     * @param index
     * @return
     */
    DataType operator[](const Point index) const;

    /**
     * initial array to @c initVal
     * @param initVal
     */
    void initArray(DataType initVal = 0);

    /**
     * save array
     * @param leap
     * @param step
     */
    void saveData(unsigned leap, unsigned step);

    /**
     *
     * @param k
     * @param leap
     * @param step
     */
    void saveData(unsigned k, unsigned leap, unsigned step);

    /**
     *
     * @param i
     * @param leap
     * @param step
     */
    void saveXPlain(unsigned i, unsigned leap, unsigned step);

    /**
     *
     * @param j
     * @param leap
     * @param step
     */
    void saveYPlain(unsigned j, unsigned leap, unsigned step);

    /**
     *
     * @param k
     * @param leap
     * @param step
     */
    void saveZPlain(unsigned k, unsigned leap, unsigned step);

    /**
     *
     * @param k
     * @param leap
     * @param step
     * @param type
     */
    void savePlain(unsigned k, unsigned leap, unsigned step, int type);

    /**
     *  Save data at plain s=@c k where s=x,y or z which define by @c type
     * @param k
     * @param leap
     * @param step
     * @param type
     */
    void save(unsigned k, unsigned leap, unsigned step, int type);

    /**
     * save every @c leap cells data to file 
     * @param leap
     */
    void save(int leap = 1);

    /**
     * @brief Create a data3d with the same size;
     * @param stru the source data3d to be copied.
     * @return
     */
    int create3DArray(const data3d< DataType > &stru);

    /**
     * @brief Create a data3d with the same size as @c stru and initial all var to @c initVal;
     * @param stru the source data3d to be copied.
     */
    int create3DArray(const data3d< DataType > &stru, DataType initVal);

    /**
     * set @name to @sn
     * @param sn
     */
    void setName(const std::string &sn) {
        mName = sn;
    }

    /**
     * get name
     * @return @c name
     */
    string getName() {
        return mName;
    }

public:

    void clearMatlabEngineArray();
    void plotArrays();
    void preparePlotting();
public:
    static int initMatlabEngine();
    static int closeMatlabEngine();
    bool isNaN(unsigned i, unsigned j, unsigned k);
    bool isInf(unsigned i, unsigned j, unsigned k);
    bool isValid(unsigned i, unsigned j, unsigned k);
    /**
     * when value at (i,j,k) is larger than limit do something define by fun
     * @param i
     * @param j
     * @param k
     * @param limit
     * @param fun
     */
    void whenLargerThan(unsigned i, unsigned j, unsigned k, MyDataF limit, void(*fun)());
private:

    static void setMatlabEngineStarted(bool MatlabEngineStarted) {
#ifdef MATLAB_SIMULATION
        mIsMatlabEngineStarted = MatlabEngineStarted;
#endif
    }
public:

    static bool isMatlabEngineStarted() {
#ifdef MATLAB_SIMULATION
        return mIsMatlabEngineStarted;
#else
        return false;
#endif
    }

};

#include "data3d.hpp"

#endif // DATA3D_H
