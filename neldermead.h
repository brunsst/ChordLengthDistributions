#ifndef NELDERMEAD_H
#define NELDERMEAD_H

#include <vector>
#include <string.h>
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <math.h>

#include "auxiliary.h"
#include "hdcommunication.h"

namespace curvefit
{
    class NelderMead
    {

    public:
        //Parameters stearing the simplex fit
        double alpha = 1.;
        double beta = 2.; // beta > alpha
        double gamma = .499; // 0 < gamma < .5
        double sigma = .5;

        uint64_t maxiter = 1000;
        double target_cost = 1e-12;
        bool bounded_cost = true;

        //Parameters for the fitting problem;
        std::vector<double> guess;
        std::vector<double> datapoints;
        std::vector<double> xaxis;
        std::vector<double> lower_bounds;
        std::vector<double> upper_bounds;

        //results of the fit
        std::vector<double> result_params;
        uint64_t result_iterations;
        double result_cost;

        /****************** Application specific functions ***************/
        bool fit_kgammafunction(std::vector<std::vector<double>> &densityestimate, double chordmean, double chordstd)
        {
            //extract the subvector to be fitted
            xaxis.assign(densityestimate[0].begin(), densityestimate[0].end());
            datapoints.assign(densityestimate[1].begin(), densityestimate[1].end());

            //define an initial guess
            guess.assign(3,0.); //mu, std
            guess[0] = chordmean;
            guess[1] = chordstd;
            guess[2] = 0.001;

            lower_bounds.assign(3, 0.);
            upper_bounds.assign(3,10000.);
            lower_bounds[0] = 1e-9;
            lower_bounds[1] = 0.1;
            lower_bounds[2] = 1e-6;

            bool successfull = run_bounded();

            //calculate goodness of fit
            double mean = std::accumulate(datapoints.begin(), datapoints.end(), 0.)/datapoints.size();
            double total_sumsquares = 0.;
            double residual_sumsquares = 0.;

            std::vector<double> y = kgamma_function(xaxis, result_params);

            for (uint64_t i = 0; i < datapoints.size(); i++)
            {
                double fitval = y[i];
                total_sumsquares += (datapoints[i] - mean) * (datapoints[i] - mean);
                residual_sumsquares += (datapoints[i] - fitval) * (datapoints[i] - fitval);
            }
            double r2 = 1.- residual_sumsquares/total_sumsquares;
            result_params.push_back(r2);

            return successfull;
        }
        bool fit_persistencelength(std::vector<std::vector<double>> &densityestimate, uint64_t modepos, double chordmean)
        {
            //extract the subvector to be fitted
            xaxis.assign(densityestimate[0].begin()+modepos, densityestimate[0].end());
            datapoints.assign(densityestimate[1].begin()+modepos, densityestimate[1].end());

            //define an initial guess
            guess.assign(3,0.); //y0, t
            guess[0] = datapoints[0];
            guess[1] = chordmean;
            guess[2] = densityestimate[0][modepos];

            bool successfull = run_unbounded();

            //calculate goodness of fit
            double mean = std::accumulate(datapoints.begin(), datapoints.end(), 0.)/datapoints.size();
            double total_sumsquares = 0.;
            double residual_sumsquares = 0.;

            for (uint64_t i = 0; i < datapoints.size(); i++)
            {
                double fitval = result_params[0] * exp(-result_params[1]*(xaxis[i]-result_params[2]));
                total_sumsquares += (datapoints[i] - mean) * (datapoints[i] - mean);
                residual_sumsquares += (datapoints[i] - fitval) * (datapoints[i] - fitval);
            }
            double r2 = 1.- residual_sumsquares/total_sumsquares;
            result_params.push_back(r2);

            return successfull;
        }
        std::vector<double> fitfunction(std::vector<double> &x, std::vector<double> &params)
        {
            std::vector<double> y(x.size(), 0.);
            for(uint64_t i = 0; i < x.size(); i++)
                y[i] = params[0]*exp(-params[1]*(x[i]-params[2]));
            return y;
        }
        std::vector<double> kgamma_function(std::vector<double> &x, std::vector<double> &params)
        {
            double mu = params[0];
            double k = (params[0]*params[0])/(params[1]*params[1]);
            double term1 = pow(k,k)/(tgamma(k)*pow(mu,k));
            double A = params[2];

            std::vector<double> y(x.size(), 0.);
            for(uint64_t i = 0; i < x.size(); i++)
                y[i] = A * term1 * pow(x[i], (k-1)) * exp(-k*x[i]/mu);
            return y;
        }
        double costfunction1(std::vector<double> &params)
        {
            double cost = 0.;
            double this_cost;
            std::vector<double> y = fitfunction(xaxis, params);
            for(uint64_t x = 0; x < xaxis.size(); x++)
            {
                this_cost = (y[x]-datapoints[x]);
                cost += this_cost*this_cost;
            }
            cost /= xaxis.size();
            return cost;
        }

    private:
        double PI = 3.14159265358979323846264338327950288419716939937510;
        std::vector<std::vector<double>> simplexes;

        /************ execute iterative fitting ******/
        bool run_unbounded();
        bool run_bounded();

        /********* simplex initialization **********/
        void _initialize_simplexes_type1();

        /************** auxiliary ******************/
        uint64_t _getsecondworsedsimplex(std::vector<double> &simplex_cost, uint64_t &worsed_idx);
    };
}

#endif // NELDERMEAD_H
