#include "algebra.h"
#include "MersenneTwister.h"
#include<vector>
#include <iomanip>
#include<iostream>
#include <cstdint>
#include <fstream>
#include <Eigen/Dense>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace Eigen;
using namespace std;

class System;
class Geometry;

struct Nodes
{
    public:

    double w;
    double D;

    Vector2d C;
    Vector2d EG;

    // constructor
    Nodes(double x,double y, double vw_dist, double weight):D(vw_dist),w(weight)
    {
        C << x,y;
        EG << 0.0,0.0;
    }
};

struct lateralEdge
{
    public:

    double l_l;
    double m_l;
    double g_l;
    double alphL;

    // constructor
    lateralEdge(double length, double tension):m_l(1.0),g_l(1.0),l_l(length),alphL(tension)
    {
        return;
    }
};

class Parameters
{
    public:

    System* system;

    int fixedVertex;
    int vertex_Number;
    int movieFrameNumber;
    int movieFrameInterval;

    double gamma;
    double epsilon;
    double pressure;
    double betaCell;

    double cellRefArea;
    double default_alphaA;
    double default_alphaB;
    double lateralTension;
    double cellAspectRatio;
    double enegryTolerance;
    double cutOffEdgeLength;

    void initialize(const char *);
    void writeToFile(string);

    bool readFromFile(const char*, bool);

    // constructor
    Parameters(System*);
};

class Switches
{
    public:

    int refVertex;
    int myoCellNumber;
    int movieFrameNumber;
    int apicalAttachment;
    int movieFrameInterval;

    double delAlpha;
    double pressure;

    void initialize(const char *);

    bool readFromFile(const char*, bool);

    // constructor
    Switches(char*);
};

class System
{
    public:

    Parameters p;

    Geometry* geometry;

    char* vary_PARAMETERS_File;
    char* constant_PARAMETERS_File;

    ofstream movieFile;

    string outputDir;

    MTRand randomGen;

    int numNode;
    int CgMaxItPerDof;

    double L;
    double H;
    double F_x;
    double F_y;
    double Lmc2;
    double lambda;
    double v_Fx[2];
    double v_Fy[2];
    double energy(void);
    double semiMajorAxis;
    double semiMinorAxis;
    double laplasPressure;
    double default_alphaL;
    double midline_radius;
    double initialTension;
    double vitellineWidth;
    double epithelialWidth;
    double epithelialHeight;
    double LmMinMaxAbsSlope;
    double flatEpithelium_refX;
    double flatEpithelium_refY;
    double CgGradientTolPerDof;
    double CgnitialStepSizePerDof;

    void saveFrame(void);
    void run(string,bool);
    void closeFiles(void);
    double checkEnergy(void);
    void updateGeometry(void);
    void updateGradient(void);
    void updateParameter(void);
    void initializeGeometry(void);
    void initializeParameter(void);
    bool checkEnergyGradient(void);
    void createOutPutDirectory(void);
    tuple<double,double> calculate_numericalBendingRigidity(void);

    void initializeOutput(string);
    void saveConfiguration(string);
    void saveEmbryoConfiguration(string);

    vector<Nodes> Va;
    vector<Nodes> Vb;
    vector<Nodes> Vv;
    vector<int> myosin_Edges;
    vector<lateralEdge> lateral_Edg;

    // constructor
    System(char*,char*,string);

    ~System();

    private:
    System(const System&);
    System& operator=(const System&);
};

class crossSectionRegion
{
    friend class Geometry;

    public:
    Geometry* const geometry;

    double A;
    double l_a;
    double l_b;
    double m_a;
    double g_a;
    double m_b;
    double g_b;
    double alphA;
    double alphB;

    // constructor
    crossSectionRegion(Geometry* g):geometry(g),l_a(0.0),m_a(1.0),g_a(1.0),l_b(0.0),m_b(1.0),g_b(1.0),A(0.0),alphA(0.0),alphB(0.0)
    {
        return;
    }
    virtual ~crossSectionRegion()
    {
        return;
    };
};

class Geometry
{
    friend class System;

    public:

    System* const system;

    double A_Y;
    double A_V;
    double A_Y_ref;

    vector<crossSectionRegion> CELL;

    // constructor
    Geometry(System* s): system(s),A_Y(0.0),A_Y_ref(0.0),A_V(0.0)
    {
        return;
    };

    virtual ~Geometry()
    {
        return;
    };
};







