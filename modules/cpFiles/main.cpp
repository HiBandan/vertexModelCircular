#include "drosoSymBreak.h"

int main(int argc, char* argv[])
{
	//***************************//
    /// define-reference-system ///
    //***************************//
    System* embryo = new System(argv[1],argv[2],"./simulationData"); embryo->createOutPutDirectory(); embryo->saveConfiguration("refFrame");    //return(0);

    //***********************************************//
    /// energy-minimization-with-default-parameters ///
    //***********************************************//
    embryo->initializeParameter();  embryo->run("noUseFile",false); embryo->saveConfiguration("initialFrame");  //return(0);

    //********************//
    /// bending-rigidity ///
    //********************//
    auto [radius,energyPerCell] = embryo->calculate_numericalBendingRigidity();
    cout << radius << " (radius), " << energyPerCell << " (energy-per-cell)" << endl; //return (0);

    //***********************************************//
    /// energy-minimization-with-updated-parameters ///
    //***********************************************//
    embryo->updateParameter();  embryo->run("timeSeriesFrames",false);  embryo->saveConfiguration("finalFrame");    //return(0);

    //****************//
    /// reset-system ///
    //****************//
    delete embryo;  embryo = NULL;

    return (0);
}
