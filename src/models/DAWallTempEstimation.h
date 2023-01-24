//
// Created by caleb on 23/01/23.
//

#ifndef JSB_STATE_DAWALLTEMPESTIMATION_H
#define JSB_STATE_DAWALLTEMPESTIMATION_H

#include "FGFDMExec.h"
#include "FGAtmosphere.h"
#include "FGAircraft.h"

namespace JSBSim {



class DAWallTempEstimation {
public:
    explicit DAWallTempEstimation(JSBSim::FGFDMExec *_fdmex);
    double GetWallTempEstimate();

private:
    FGFDMExec* fdmex_;
    std::shared_ptr<FGAtmosphere> atmosphere_;
    std::shared_ptr<FGAircraft> aircraft_;

    double HeatBalance(double estimate);
    std::tuple<double, double> CalculateCpAndGamma(double estimate);
    double NewtonRaphson(double x, double xmin=10, double xmax=1000);

    double machSpeed_;
    double altitude_;
    double emissivity_ = 0.9;
    double noseDistance_;
    int flowType_ = 1;

};
}
#endif //JSB_STATE_DAWALLTEMPESTIMATION_H
