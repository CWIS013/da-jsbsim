//
// Created by caleb on 23/01/23.
//
#include "FGFDMExec.h"
#include "DAWallTempEstimation.h"

using namespace std;

namespace JSBSim {

    DAWallTempEstimation::DAWallTempEstimation(FGFDMExec* _fdmex){
      fdmex = _fdmex;
      atmosphere = fdmex->GetAtmosphere();
    }

    double DAWallTempEstimation::GetWallTempEstimate() {
      machSpeed = 1.8; //fdmex->GetPropertyValue("velocities/mach");
      return DAWallTempEstimation::NewtonRaphson(500);
    }


    double DAWallTempEstimation::HeatBalance(double estimate) {
      double Ta = atmosphere->GetTemperature()/1.8;
      double pa = atmosphere->GetPressure() * 47.880208;
      double nu = atmosphere->GetKinematicViscosity();
      double cp, gm;
      tie(cp, gm) = CalculateCpAndGamma(estimate);
      double Pr = atmosphere->GetAbsoluteViscosity() * cp / CalculateK();
      double sg = 5.67E-8;
      double w = 0.76;
      double aa = atmosphere->GetSoundSpeed()*0.3048;
      double Va = machSpeed * aa;

      double Re = Va * noseDistance / nu;
      double C;
      double N;
      double r;
      double a;
      double b;
      if (flowType == 0) {
        C = 0.664;
        N = 0.5;
        r = pow(Pr, 0.5);
        a = 0.032;
        b = 0.58;
      }
      else {
        C = 0.0592;
        N = 0.2;
        r = pow(Pr, (1 / 3));
        a = 0.035;
        b = 0.45;
      }


      double Twad = Ta * (1 + 0.5 * r * (gm - 1) * pow(machSpeed, 2)); //10 degree variance
      double q = 0.5 * gm * pa * pow(machSpeed, 2);
      double cfi = C / (pow(Re, N));
      double NN = 1 - N * (w + 1);
      double Tref = 1 + a * pow(machSpeed, 2) + b * (estimate / Ta - 1);
      double cf = cfi / (pow(Tref, NN));
      double St = 0.5 * cf / pow(Pr, (2 / 3));

      double conducted = St * r * q * Va * (Twad - estimate) / (Twad - Ta);
      double radiated = 0;
      if (machSpeed > 1.4) {
        radiated = emissivity * sg * pow((estimate - Ta), 4);
      }
      double balance = conducted - radiated;
      return balance;
    }

    double DAWallTempEstimation::CalculateK() {
      double T00 = 288.15;
      double T11 = 216.65;
      double T32 = 228.65;
      double T47 = 270.65;
      double T61 = 252.65;
      double T79 = 180.65;
      double Labda;
      double Trel;
      double T0;
      if (altitude <= 11E3) {
        Labda = -6.5E-3;
        Trel = 1 + Labda * altitude / T00;
        T0 = Trel * T00;
      } else if (altitude <= 20E3){
        T0 = T11;
      } else if (altitude <= 32E3){
        Labda = 1E-3;
        Trel = 1 + Labda * (altitude - 20E3) / T11;
        T0 = Trel * T11;
      } else if (altitude <= 47E3){
        Labda = 2.8E-3;
        Trel = 1 + Labda * (altitude - 32E3) / T32;
        T0 = Trel * T32;
      } else if (altitude <= 52E3){
        T0 = T47;
      } else if (altitude <= 61E3){
        Labda = -2E-3;
        Trel = 1 + Labda * (altitude - 52E3) / T47;
        T0 = Trel * T47;
      } else if ( altitude <= 79E3){
        Labda = -4E-3;
        Trel = 1 + Labda * (altitude - 61E3) / T61;
        T0 = Trel * T61;
      } else if (altitude <= 88.7E3){
        T0 = T79;
      } else if (altitude <= 120E3) {
        T0 = T79;
      }else {
        T0 = T79;
      }
      return 2.64638e-3 * pow(T0, 1.5) / (T0 + 245 * pow(10, (-12 / T0)));
    }

    std::tuple<double, double> DAWallTempEstimation::CalculateCpAndGamma(double estimate) {
      vector<int> T_arr = {250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500};
      vector<double> cp_arr = {1.003, 1.005, 1.008, 1.013, 1.02, 1.029, 1.04, 1.051, 1.063, 1.075, 1.087, 1.099, 1.121, 1.142, 1.155,
              1.173, 1.19, 1.204, 1.216};
      vector<double> cv_arr = {0.716, 0.718, 0.721, 0.726, 0.733, 0.742, 0.753, 0.764, 0.776, 0.788, 0.8, 0.812, 0.834, 0.855, 0.868,
              0.886, 0.903, 0.917, 0.929};
      double cp;
      double cv;
      int i = 0;
      for (int n =0; n< T_arr.size(); n++) {
        if (estimate < T_arr[n]) {
          i = n;
          break;
        }
      }
      if (i == 0) {
        cp = cp_arr[19];
        cv = cv_arr[19];
      }else if (i == 1){
        cp = cp_arr[1];
        cv = cv_arr[1];
      } else {
        double fraction = (estimate - T_arr[i - 1]) / (T_arr[i] - T_arr[i - 1]);
        cp = cp_arr[i - 1] + fraction * (cp_arr[i] - cp_arr[i - 1]);
        cv = cv_arr[i - 1] + fraction * (cv_arr[i] - cv_arr[i - 1]);
      }
      cp = 1000 * cp;
      cv = 1000 * cv;
      double gamma = cp / cv;
      return {cp, gamma};
    }

    double DAWallTempEstimation::NewtonRaphson(double x, double xmin, double xmax) {
      double tol=1E-2, dx=1E-1, maxiter=30, xs = x, x0, x1, f0, f1, df, fx;
      for (int i=0; i<maxiter; i++) {
        x0 = xs - dx;
        x1 = xs + dx;

        f0 = HeatBalance(x0);
        f1 = HeatBalance(x1);
        df = f1 - f0;
        df = 0.5 * df / dx;
        fx = HeatBalance(xs);
        if (abs(fx) < tol) {
          break;
        }
        if (abs(df) < 1E-16) {
          dx = 10 * dx;
        } else {
          xs = xs - fx / df;
        }
//        if (xs < xmin) {
//          xs = xmin;
//          break;
//        }
//        if (xs > xmax) {
//          xs = xmax;
//          break;
//        }
      }
      return xs;
    }
}