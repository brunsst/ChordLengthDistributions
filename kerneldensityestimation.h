#ifndef KERNELDENSITYESTIMATION_H
#define KERNELDENSITYESTIMATION_H

#include <vector>
#include <string.h>
#include <iostream>
#include <numeric>
#include <math.h>

using namespace std;

namespace analysis
{
    class KDE
    {
    public:
        double bandwidth = 1.;

        vector<vector<double>> run_kerneldensityestimation(vector<double> &yvalues, string kernel_function, string bandwidth_estimator, double stepsize)
        {
            std::sort(yvalues.begin(), yvalues.end());

            double N = yvalues.size();
            double minvalue = yvalues[0];
            double maxvalue = yvalues[N-1];
            double sum = std::accumulate(yvalues.begin(), yvalues.end(), 0.0);
            double mean = sum / N;
            double sq_sum = std::inner_product(yvalues.begin(), yvalues.end(), yvalues.begin(), 0.0);
            double stdev = std::sqrt(sq_sum / N - mean * mean);

            //set bandwidth
            if (bandwidth_estimator == "IQR")
            {
                double interquartile_range = yvalues[int(N*0.75)]-yvalues[int(N*0.25)];
                double A = min(stdev, interquartile_range);
                bandwidth = 0.9*A*pow(N,-0.2);
            }
            else if (bandwidth_estimator == "Gaussian")
                bandwidth = 1.06*stdev*pow(N,-0.2);

            vector<double> out_xaxis, out_yaxis;
            uint64_t first_valid_pos = 0;

            for (double xpos = floor(minvalue); xpos <= maxvalue; xpos += stepsize)
            {
                //estimate density at x-position
                double p = 0.;

                for (uint64_t ypos = first_valid_pos; ypos < N; ypos++)
                {
                    double u = (xpos-yvalues[ypos])/bandwidth;

                    //Apply kernel
                    if (kernel_function == "Epanechnikov")
                    {
                        if (u >= 1)
                            first_valid_pos = ypos;
                        else if (fabs(u) < 1)
                            p += 0.75*(1.-u*u);
                        else
                            break;
                    }
                }

                p = p/(N*bandwidth);

                out_xaxis.push_back(xpos);
                out_yaxis.push_back(p);
            }

            vector<vector<double>> out_kde;
            out_kde.push_back(out_xaxis);
            out_kde.push_back(out_yaxis);
            return out_kde;
        }
    };
}

#endif // KERNELDENSITYESTIMATION_H

