//
// Created by caleb on 23/01/23.
//

#ifndef JSB_STATE_DAWALLTEMPESTIMATION_H
#define JSB_STATE_DAWALLTEMPESTIMATION_H

#include "FGFDMExec.h"
#include "FGAtmosphere.h"

namespace JSBSim {



class DAWallTempEstimation {
public:
    DAWallTempEstimation(JSBSim::FGFDMExec *_fdmex);
    double GetWallTempEstimate();

private:
    FGFDMExec* fdmex;
    std::shared_ptr<FGAtmosphere> atmosphere;

    double HeatBalance(double estimate);
    double CalculateK();
    std::tuple<double, double> CalculateCpAndGamma(double estimate);
    double NewtonRaphson(double x, double xmin=10, double xmax=1000);

    double Earth_g0 = 9.80665;
    double Earth_radius = 6378137;
    double Earth_mu = 398600.44189 * 1e9;
    double air_R = 8.3145 / 0.0289645;

    double machSpeed;
    double altitude;


    double emissivity = 0.9;
    double noseDistance = 1.0;
    int flowType = 1;

};
}
#endif //JSB_STATE_DAWALLTEMPESTIMATION_H
