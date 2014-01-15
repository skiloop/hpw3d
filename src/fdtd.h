#pragma once

#include "cpml.h"
#include "Point.h"
#include "source.h"
#include "ConnectingInterface.h"

class fdtd {
public:
    //#ifdef WITH_DENSITY
    fdtd(int useDensity, unsigned _totalTimeSteps = 500, unsigned _imax = 40, unsigned _jmax = 40, unsigned _kmax = 26,
            MyDataF _tw = 53.0e-12, MyDataF _dx = 1e-3, MyDataF _dy = 1e-3, MyDataF _dz = 1e-3,
            MyDataF _amp = 1000, unsigned _savemodulus = 100, unsigned _ksource = 12,
            unsigned _m = 3, unsigned _ma = 1, unsigned pmlw = 6, int useConnect = 0,
            unsigned _neGrid = 16, MyDataF maxNe = DEFAULT_DENSITY_MAX);

    void setPlasmaParam(MyDataF _rei, MyDataF _vm, MyDataF _p, int _ftype);
    //#else
    //    fdtd(unsigned _totalTimeSteps = 500, unsigned _imax = 40, unsigned _jmax = 40, unsigned _kmax = 26,
    //            MyDataF _tw = 53.0e-12, MyDataF _dx = 1e-3, MyDataF _dy = 1e-3, MyDataF _dz = 1e-3,
    //            MyDataF _amp = 1000, unsigned _savemodulus = 100, unsigned _ksource = 12,
    //            unsigned _m = 3, unsigned _ma = 1, unsigned pmlw = 6, int useConnect = 0);
    //#endif
    ~fdtd(void);

    static const int SOURCE_GAUSSIAN = GAUSSIAN_WAVE;
    static const int SOURCE_SINE = SINE_WAVE;
    static const int SOURCE_DERIVE_GAUSSIAN = DERIVATIVE_GAUSSIAN_WAVE;
    static const int SOURCE_ZERO = ZERO_TYPE;

    //Function prototype definitions
    void createFieldArray(); //Memory initialization

    /**
     * 
     * @param t_0
     * @param omega_
     * @param tUp
     * @param tDown
     * @param amptidute
     */
    void intSourceSinePulse(MyDataF t_0, MyDataF omega_, MyDataF tUp, MyDataF tDown, MyDataF amptidute);

    void setUp(); //Coefficients, parameters etc will get computed

    void compute(); //E & H Field update equation

    void startUp();

    /**
     * 
     * @param omega_
     */
    void defineSineSource(MyDataF omega_);

    /**
     * 
     * @param sourceType
     */
    void setSourceType(sourceType* pSrcType);

    void setSource(source*p) {
        pSource = p;
    };

    void setSrcType(int srcType);

    void desideDomainZone();

#ifdef MATLAB_SIMULATION
    int initMatlabSimulation();
    void doMatlabSimulation();
    void finishMatlabSimulation();
#endif

private:
    void buildObject(); //Creates the object geometry
    void yeeCube(unsigned, unsigned, unsigned, unsigned); //Sets material properties to a cell
    void writeField(unsigned); //Writes output
    void buildSphere(); //Builds a spherical object
    void buildDipole(); //Builds a dipole
    void printParam();

    int mIsUseDensity;
private:

    //  Specify Number of Time Steps and Grid Size Parameters
    unsigned mTotalTimeSteps; // total number of time steps

    // grid size corresponding to the number of Ez field components   
    Point mMaxIndex;

    //  Specify the Impulsive Source (Differentiated Gaussian) parameters
    MyDataF tw; //pulse width
    MyDataF mDt, mDx, mDy, mDz;

    MyDataF mAmplitude; // Amplitude
    // Specify the Time Step at which the data has to be saved for Visualization
    unsigned mSaveModulus;

    // Output recording 
    unsigned mKSource;
    // Specify the CPML Order and Other Parameters:
    unsigned mPMLOrder;
    unsigned mAlphaOrder;
    unsigned mPMLWidth;
    unsigned mAirBufferWidth;

    /////////////////////////////////////////////////////////
    // SOURCE PARAMETRS
    ////////////////////////////////////////////////////////
    MyDataF t_up; // up bound of Pulse for Sine Pulse
    MyDataF t_down; // down bound of Pulse for Sine Pulse

    //#ifdef WITH_DENSITY
    //how many fine grids per coarse grid 
    unsigned mNeGridSize;
    //initial plasma value
    MyDataF Ne0;
    //#endif
    // source type
    int mSrcType;
    // if use connecting interface
    bool mUseConnectingInterface;
    //permittivity, permeability and conductivity of different materials
    MyDataF *mEpsilon;
    MyDataF *mSigma;
    MyDataF *mMu;


    MyDataF t0; //delay    
    MyDataF mOmega; //wave angle frequency

    //  Specify the dipole Boundaries(A cuboidal rode- NOT as a cylinder)
    Point mNonPMLStartIndex;
    Point mNonPMLEndIndex;
    Point mDomainStartIndex;
    Point mDomainEndIndex;
    // source position
    Point mSourceIndex;
    Point mNeSrcPos;

    // Source 
    source *pSource;

    // source type
    sourceType *pSourceType;

    // common data
    MyDataF dtDivEps0DivDxyz;

    // H & E Field components
    data3d<MyDataF> Hx;
    data3d<MyDataF> Hy;
    data3d<MyDataF> Hz;
    data3d<MyDataF> Ex;
    data3d<MyDataF> Ey;
    data3d<MyDataF> Ez;

    data3d<unsigned> ID1; //medium definition array for Ex
    data3d<unsigned> ID2; //medium definition array for Ey
    data3d<unsigned> ID3; //medium definition array for Ez

    //E field update coefficients
    MyDataF *CA;
    MyDataF *CB;

    //H field update coefficients
    MyDataF DA;
    MyDataF DB;

    data3d<MyDataF> Cexe, Ceye, Ceze;
    data3d<MyDataF> Chxh, Chyh, Chzh;
    data3d<MyDataF> Cexhy, Ceyhz, Cezhx, Cexhz, Ceyhx, Cezhy;
    data3d<MyDataF> Chxey, Chyez, Chzex, Chxez, Chyex, Chzey;

    void initCoeficients();

    //#ifdef WITH_DENSITY
    // neutral gas density in cm^-3
    static const MyDataF mNeutralGasDensity;

    int mNiuType;
    //Fine Grid size
    MyDataF mDsFluid;
    //time step of plasma
    MyDataF mDtFluid;
    int mNeSkipStep;
    int mNeBoundWidth;

    // temporary variables that often used
    MyDataF mHalfDelta_t;
    MyDataF mHalf_e;
    MyDataF eMDtDiv2DivEps0;
    MyDataF dtDivEps0DivDx;
    MyDataF e2Dt2Div4DivEps0DivMe;
    MyDataF dtDivEps0DivDy;
    MyDataF dtDivEps0DivDz;
    MyDataF mCoeffVelocity;
    unsigned mHalfNeGridSize;

    //plasma variables
    MyDataF mNu_m; //collision frequency
    MyDataF mAirPressure; // air pressure
    MyDataF mDe;
    MyDataF mDa;
    MyDataF mRei;
    MyDataF mMu_i;
    MyDataF mMu_e;
    MyDataF mGamma;
    MyDataF mA;
    MyDataF mAlpha;

    // update coefficients
    //    MyDataF Chxey,Chxez,Chyez,Chyex,Chzex,Chzey;
    data3d<MyDataF> Cezvz, Ceyvy, Cexvx;
    data3d<MyDataF> Cvzvz, Cvyvy, Cvxvx;
    data3d<MyDataF> Cvxex, Cvyey, Cvzez;
    MyDataF mCvxex, mCvyey, mCvzez;

    //Plasma
    data3d<MyDataF> Ne, Ne_pre;
    //
    data3d<MyDataF> Eeff;
    // Beta
    data3d<MyDataF> Beta;

    // collision frequency for reused mode
    data3d<MyDataF> Nu_c;
    //velocity
    data3d<MyDataF> Vx;
    data3d<MyDataF> Vy;
    data3d<MyDataF> Vz;

    void createDensityArrays();
    //initials
    void initCoeffForDensity();
    void initDensity();
    void createCoeff();
    void updateCoeffWithDensity();
    void updateBeta();

    // Erms or Eeff operation
    void captureEFieldForEeff(void);
    void updateCollisionFrequency();
    void vec2Eeff();
    void updateEeff();
    void updateDensity(void);
    void updateVelocity(void);
    void wallCircleBound(data3d<MyDataF> &stru);
    //#endif
    void updateHx();
    void updateHy();
    void updateHz();
    void updateMagneitcFields();
    void updateEx();
    void updateEy();
    void updateEz();
    void updateElectricAndVeloityFields();
    void updateSource(unsigned n);
    cpml<MyDataF> mPML;
    ConnectingInterface mConnectingInterface;
private:

    /**
     * Ionization parameters at point (i,j,k) for density updating
     * @param i
     * @param j
     * @param k
     * @param va
     * @param vi
     * @param Deff
     */
    void calIonizationParam(int i, int j, int k, MyDataF &va, MyDataF &vi, MyDataF &Deff);
    
    void updateSourceCoeff();
    void updateSourceCoeff(const Point&nes,const Point&bes);
#ifdef DEBUG
    ofstream mOfNeCheck;
#endif
};


