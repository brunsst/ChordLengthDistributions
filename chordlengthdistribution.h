#ifndef CHORDLENGTHDISTRIBUTION_H
#define CHORDLENGTHDISTRIBUTION_H

#include <vector>
#include <string.h>
#include <iostream>
#include "hdcommunication.h"

using namespace std;

namespace analysis
{
    class CLD
    {
    public:
        vector<vector<vector<int64_t>>> shift_lut; //dim0: chord direction, dim1: step, dim2: x, y, z, idx shifts

        double angular_precision = 1e-3;
        bool inplane_only = false; //set to true for 2D evaluation

        int color_void, color_material, color_interface;
        double spacing;

        CLD(double angular_spacing, int shape[3], uint8_t color_void_, uint8_t color_material_, uint8_t color_interface_)
        {
            spacing = angular_spacing;

            srand (time(NULL));

            color_void = int(color_void_);
            color_material = int(color_material_);
            color_interface = int(color_interface_);
            n_row = shape[0];
            n_slice = shape[0]*shape[1];

            //Calculate LUT for all possible chords:
            shift_lut.clear();
            double mininclination = 0.;
            double maxinclination = 180.;
            if (inplane_only){
                mininclination = 90.;
                maxinclination = 90.+eps;
            }

            for(double inclination =mininclination; inclination < maxinclination; inclination += angular_spacing)
            {
                for(double azimuth = 0.; azimuth < 180.; azimuth += angular_spacing)
                {
                    calculate_idxshifts(azimuth, inclination, shape);

                    if (inclination == 0)
                        break;
                }
            }

            cout << shift_lut.size() << " chord directions" << endl;
        }

        vector<double> evaluate_chords_alldirections(uint8_t* imagestack, int shape[3], int &x, int &y, int &z, int &phasevalue);
        double evaluate_singlerandomchord(uint8_t* imagestack, int shape[3], int &x, int &y, int &z, int &phasevalue);
        double evaluate_singledirection(uint8_t* imagestack, int shape[3], int x, int y, int z, uint64_t &direction, int phasevalue);
        void evaluate_chords_alldirections(uint8_t* imagestack, int shape[3], int &x, int &y, int &z, int &phasevalue, vector<vector<double>> &valid_chords);

        double get_harmonic_mean(vector<double> &chordlist);
        double get_geometric_mean(vector<double> &chordlist);

    private:
        double eps = 1.e-9;
        double pi = 3.141592653589793238462643383279502884197169399375105820974944;

        uint64_t n_row, n_slice;

        void calculate_idxshifts(double phi_deg, double theta_deg, int shape[3]);
        bool _forwardevaluation(int out_endpoint[3], vector<vector<int64_t>> &chord, uint8_t* imagestack, int shape[3],
                                const int &x, const int &y, const int &z, const int &phasevalue);
        bool _backwardevaluation(int out_endpoint[3], vector<vector<int64_t>> &chord, uint8_t* imagestack, int shape[3],
                                const int &x, const int &y, const int &z, const int &phasevalue);
    };
}

#endif // CHORDLENGTHDISTRIBUTION_H


