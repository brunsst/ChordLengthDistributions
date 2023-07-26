#include "chordlengthdistribution.h"
#include <math.h>
using namespace std;

namespace analysis
{
    /**************************** Analysis ***********************************/
    double CLD::evaluate_singlerandomchord(uint8_t* imagestack, int shape[3], int &x, int &y, int &z, int &phasevalue)
    {
        /*
         * returns the length of a chord in a random direction; -1 when chord is invalid
         */

        //choose a direction
        uint64_t direction = rand() % shift_lut.size();
        vector<vector<int64_t>> this_chord = shift_lut[direction];

        //evaluate the chord
        int endpoint1[3];
        int endpoint2[3];

        bool valid = _forwardevaluation(endpoint1, this_chord, imagestack, shape, x, y, z, phasevalue);
        if (valid == false) return -1.;

        valid = _backwardevaluation(endpoint2, this_chord, imagestack, shape, x, y, z, phasevalue);
        if (valid)
        {
            //measure the chord and return it
            double chordlength = (endpoint1[0]-endpoint2[0])*(endpoint1[0]-endpoint2[0]);
            chordlength += (endpoint1[1]-endpoint2[1])*(endpoint1[1]-endpoint2[1]);
            chordlength += (endpoint1[2]-endpoint2[2])*(endpoint1[2]-endpoint2[2]);
            return sqrt(chordlength);
        }
        else return -1.;
    }
    double CLD::evaluate_singledirection(uint8_t* imagestack, int shape[3], int x, int y, int z, uint64_t &direction, int phasevalue)
    {
        /*
         * returns the length of a chord in a given direction; -1 when chord is invalid
         */

        int x0 = x;
        int y0 = y;
        int z0 = z;

        //choose a direction
        vector<vector<int64_t>> this_chord = shift_lut[direction];

        //evaluate the chord
        int endpoint1[3];
        int endpoint2[3];

        bool valid = _forwardevaluation(endpoint1, this_chord, imagestack, shape, x0, y0, z0, phasevalue);
        if (valid == false) return -1.;

        valid = _backwardevaluation(endpoint2, this_chord, imagestack, shape, x0, y0, z0, phasevalue);
        if (valid)
        {
            //measure the chord and return it
            double chordlength = (endpoint1[0]-endpoint2[0])*(endpoint1[0]-endpoint2[0]);
            chordlength += (endpoint1[1]-endpoint2[1])*(endpoint1[1]-endpoint2[1]);
            chordlength += (endpoint1[2]-endpoint2[2])*(endpoint1[2]-endpoint2[2]);
            return sqrt(chordlength);
        }
        else return -1.;
    }

    void CLD::evaluate_chords_alldirections(uint8_t* imagestack, int shape[3], int &x, int &y, int &z, int &phasevalue, vector<vector<double>> &valid_chords)
    {
        /*
         * returns a vector with the length of all valid chords at a given location
         */

        for (uint64_t chord = 0; chord < shift_lut.size(); chord++)
        {
            vector<vector<int64_t>> this_chord = shift_lut[chord];
            int endpoint1[3];
            int endpoint2[3];

            bool valid = true;

            //forward evaluation
            valid = _forwardevaluation(endpoint1, this_chord, imagestack, shape, x, y, z, phasevalue);

            if (valid == false)
            {
                continue;
            }

            //backward evaluation
            valid = _backwardevaluation(endpoint2, this_chord, imagestack, shape, x, y, z, phasevalue);

            if (valid)
            {
                //measure the chord and keep it
                double chordlength = (endpoint1[0]-endpoint2[0])*(endpoint1[0]-endpoint2[0]);
                chordlength += (endpoint1[1]-endpoint2[1])*(endpoint1[1]-endpoint2[1]);
                chordlength += (endpoint1[2]-endpoint2[2])*(endpoint1[2]-endpoint2[2]);
                chordlength = sqrt(chordlength);
                valid_chords[chord].push_back(chordlength);
            }
        }
        return;
    }
    vector<double> CLD::evaluate_chords_alldirections(uint8_t* imagestack, int shape[3], int &x, int &y, int &z, int &phasevalue)
    {
        /*
         * returns a vector with the length of all valid chords at a given location
         */

        vector<double> valid_chords; //first value identifies the phase

        for (uint64_t chord = 0; chord < shift_lut.size(); chord++)
        {
            vector<vector<int64_t>> this_chord = shift_lut[chord];
            int endpoint1[3];
            int endpoint2[3];

            bool valid = true;

            //forward evaluation
            valid = _forwardevaluation(endpoint1, this_chord, imagestack, shape, x, y, z, phasevalue);

            if (valid == false)
                continue;

            //backward evaluation
            valid = _backwardevaluation(endpoint2, this_chord, imagestack, shape, x, y, z, phasevalue);

            if (valid)
            {
                //measure the chord and keep it
                double chordlength = (endpoint1[0]-endpoint2[0])*(endpoint1[0]-endpoint2[0]);
                chordlength += (endpoint1[1]-endpoint2[1])*(endpoint1[1]-endpoint2[1]);
                chordlength += (endpoint1[2]-endpoint2[2])*(endpoint1[2]-endpoint2[2]);
                chordlength = sqrt(chordlength);
                valid_chords.push_back(chordlength);
            }
        }
        return valid_chords;
    }

    /**************************** LUT generation ***********************************/
    void CLD::calculate_idxshifts(double phi_deg, double theta_deg, int shape[3])
    {
        /*
         * This function precalculates the trajectory of a given chord as a LUT
         */

        //to radian
        double phi_rad = pi/180.*phi_deg;
        double theta_rad = pi/180.*theta_deg;

        //origin for LUT
        double chordpos[3] = {0,0,0};
        vector<int64_t> shift = {0,0,0,0}; //x,y,z and idx shift
        double chordlength = 0.;

        //3D chord direction
        double direction[3] = {cos(phi_rad)*sin(theta_rad), sin(phi_rad)*sin(theta_rad), cos(theta_rad)};

        //std::cout << phi_deg << " / " << theta_deg << ": " << direction[0] << ",  " << direction[1] << ",  " << direction[2] << std::endl;

        double maxlength = sqrt(shape[0]*shape[0]+shape[1]*shape[1]+shape[2]*shape[2]);

        vector<vector<int64_t>> this_chord;
        this_chord.push_back(shift);

        while(chordlength <= maxlength)
        {
            chordlength += 1.;

            chordpos[0] = chordlength*direction[0];
            chordpos[1] = chordlength*direction[1];
            chordpos[2] = chordlength*direction[2];

            int new_shift[3] = {round(chordpos[0]), round(chordpos[1]), round(chordpos[2])};

            if ((new_shift[0] != shift[0]) || (new_shift[1] != shift[1]) ||(new_shift[2] != shift[2]))
            {
                shift[0] = new_shift[0];
                shift[1] = new_shift[1];
                shift[2] = new_shift[2];
                shift[3] = shift[0] + shift[1]*shape[0] + shift[2]*shape[0]*shape[1];

                this_chord.push_back(shift);
            }
        }

        shift_lut.push_back(this_chord);
    }

    /********************* Evaluation of specific means ****************************/
    double CLD::get_harmonic_mean(vector<double> &chordlist)
    {
        double harmonic_mean = 0.;
        double counter = 0;
        for (uint64_t i = 0; i < chordlist.size(); i++)
        {
            if (chordlist[i] != 0)
            {
                harmonic_mean += 1./chordlist[i];
                counter++;
            }
        }
        harmonic_mean = counter/harmonic_mean;
        return harmonic_mean;
    }
    double CLD::get_geometric_mean(vector<double> &chordlist)
    {
        double geometric_mean = 0.;
        double counter = 0;
        for (uint64_t i = 0; i < chordlist.size(); i++)
        {
            if (chordlist[i] != 0)
            {
                geometric_mean += log(chordlist[i]);
                counter++;
            }
        }
        geometric_mean = exp(geometric_mean/counter);

        return geometric_mean;
    }

    /**************************** Forward and backward evaluation of a single chord ***********************************/
    bool CLD::_forwardevaluation(int out_endpoint[3], vector<vector<int64_t>> &chord, uint8_t* imagestack, int shape[3],
                                const int &x, const int &y, const int &z, const int &phasevalue)
    {
        bool valid = true;
        int64_t origin_idx = x + y*n_row + z*n_slice;

        //forward evaluation
        for(uint64_t step = 1; step < chord.size(); step++)
        {
            vector<int64_t> this_step = chord[step];

            //check if forward step out of bounds
            if    ((x+this_step[0] < 0) || (x+this_step[0] >= shape[0])
            || (y+this_step[1] < 0) || (y+this_step[1] >= shape[1])
            || (z+this_step[2] < 0) || (z+this_step[2] >= shape[2]))
            {
                valid = false;
                break;
            }

            int this_color = int(imagestack[origin_idx+this_step[3]]);

            if (this_color == color_interface)
            {
                //decide if the chord passes or continues
                bool passed = rand() % 2;
                if (passed == 0)
                {
                    //hitting the other phase
                    out_endpoint[0] = x + chord[step-1][0];
                    out_endpoint[1] = y + chord[step-1][1];
                    out_endpoint[2] = z + chord[step-1][2];
                    break;
                }
            }
            else if (this_color != phasevalue)
            {
                //hitting the other phase
                out_endpoint[0] = x + chord[step-1][0]; //taking one step back on forward evaluation
                out_endpoint[1] = y + chord[step-1][1];
                out_endpoint[2] = z + chord[step-1][2];
                break;
            }
        }

        return valid;
    }
    bool CLD::_backwardevaluation(int out_endpoint[3], vector<vector<int64_t>> &chord, uint8_t* imagestack, int shape[3],
                                const int &x, const int &y, const int &z, const int &phasevalue)
    {
        bool valid = true;
        int64_t origin_idx = x + y*n_row + z*n_slice;

        //forward evaluation
        for(uint64_t step = 1; step < chord.size(); step++)
        {
            vector<int64_t> this_step = chord[step];

            //check if forward step out of bounds
            if    ((x-this_step[0] < 0) || (x-this_step[0] >= shape[0])
            || (y-this_step[1] < 0) || (y-this_step[1] >= shape[1])
            || (z-this_step[2] < 0) || (z-this_step[2] >= shape[2]))
            {
                valid = false;
                break;
            }

            int this_color = int(imagestack[origin_idx-this_step[3]]);

            if (this_color == color_interface)
            {
                //decide if the chord passes or continues
                bool passed = rand() % 2;
                if (passed == 0)
                {
                    //hitting the other phase
                    out_endpoint[0] = x - chord[step][0];
                    out_endpoint[1] = y - chord[step][1];
                    out_endpoint[2] = z - chord[step][2];
                    break;
                }
            }
            else if (this_color != phasevalue)
            {
                //hitting the other phase
                out_endpoint[0] = x - chord[step][0]; //not taking one step back on backward evaluation
                out_endpoint[1] = y - chord[step][1];
                out_endpoint[2] = z - chord[step][2];
                break;
            }
        }

        return valid;
    }
}
