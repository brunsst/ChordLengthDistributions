#include <iostream>
#include <omp.h>
#include <ctime>
#include "hdcommunication.h"
#include "chordlengthdistribution.h"
#include "auxiliary.h"
#include "kerneldensityestimation.h"
#include "neldermead.h"
#include <unistd.h>

using namespace std;

int main(int argc, char* argv[])
{
    /*
     * For diagonal chords a binsize of does not work. Use an even moving window to compensate
     */

    //"average": records average chord length for all directions sampled
    //"random": fastest mode, just probes one random direction
    //"minchord": evaluates all directions and keeps smallest chord
    //"directionality": dumped from spring8_directionality. Calculates a CLD for all available directions.
    string mode = "random"; //"random", "average", "minchord", "directionality"

    vector<string> tasks;

    uint64_t n_samples = 10000;

    double voxelsize = 1.0; //in nm//10.79;// //0.0405;

    uint8_t color_material = 255;
    uint8_t color_void = 0;
    uint8_t color_interface = 128;

    double angular_spacing = 15;//5.;

    //Kernel density estimate:
    string bandwidth_estimator = "IQR"; //"IQR", "Gaussian"
    string kernel = "Epanechnikov";
    double stepsize = 1.;

    //persistence_length
    int64_t persistence_offset = 10; //start the fit right of the mode

    bool material_only = false;
    bool fit_persistencelength = false;
    bool fit_kgamma = false; //not working yet
    bool export_allchords = false;

    ////////////////////////////////////////

    if ("extract comman line arguments)"){
        for (uint16_t i = 1; i < argc; i++)
        {
            if ((string(argv[i]) == "-i") || (string(argv[i]) == "-input"))
            {
                i++;
                tasks.push_back(string(argv[i]));
            }
            else if ((string(argv[i]) == "-vxl") || (string(argv[i]) == "-voxelsize"))
            {
                i++;
                voxelsize = atof(argv[i]);
            }
            else if ((string(argv[i]) == "-n") || (string(argv[i]) == "-samples"))
            {
                i++;
                n_samples = atoi(argv[i]);
            }
            else if (string(argv[i]) == "-mode")
            {
                i++;
                mode = string(argv[i]);
            }
            else if ((string(argv[i]) == "-angle"))
            {
                i++;
                angular_spacing = atof(argv[i]);
            }
            else if ((string(argv[i]) == "-cm") || (string(argv[i]) == "-color_material"))
            {
                i++;
                color_material = atof(argv[i]);
            }
            else if ((string(argv[i]) == "-cv") || (string(argv[i]) == "-color_void"))
            {
                i++;
                color_void = atof(argv[i]);
            }
            else if ((string(argv[i]) == "-ci") || (string(argv[i]) == "-color_interface"))
            {
                i++;
                color_interface = atof(argv[i]);
            }
        }
    }

    voxelsize /= 1000; //to um

    hdcom::HdCommunication hdcom;
    srand (time(NULL));

    for (uint16_t task = 0; task < tasks.size(); task++)
    {
        cout << "CLD for " << tasks[task] << endl;
        cout << "##################################" << endl;

        string rootpath = tasks[task].substr(0, tasks[task].rfind("/", tasks[task].length()-2)+1);
        int shape[3];
        std::vector<std::string> filelist = hdcom.GetFilelist(tasks[task], shape);
        uint8_t* imagestack = hdcom.Get3DTifSequence_8bitPointer(filelist, shape, true);

        analysis::CLD cld(angular_spacing, shape, color_void, color_material, color_interface);

        if (mode == "directionality") // <-- section simply copied from spring8_directionality to give Florian a allinone-solution
        {
           //3D-TwoPointCorrelationFunction --> test in 2D
            //////////////////////////////////////////////////////////////
            int nx = shape[0]; int ny = shape[1]; int nz = shape[2];
            analysis::KDE kde;

            int64_t n_directions = cld.shift_lut.size();

            vector<vector<double>> collected_chords;
            vector<double> dummy;
            for (int64_t i = 0; i < n_directions; i++)
                collected_chords.push_back(dummy);

            int64_t minvalid = 0;
            int64_t iter = 0;
            int n_locations = n_samples;


            while (minvalid < n_locations)
            {
                iter++;
                cout << "iter " << iter << "    \r";
                cout.flush();

                for (int64_t n = 0; n < n_locations/2; n++)
                {
                    if ((n+1) % 100 == 0)
                    {
                        cout << "iter " << iter << ": " << n << "     \r";
                        cout.flush();
                    }
                    //probe the image at desired amount of locations
                    int x1 = rand() % nx;
                    int y1 = rand() % ny;
                    int z1 = rand() % nz;
                    int64_t idx1 = z1*nx*ny+ y1*nx + x1;

                    if (imagestack[idx1] == color_void)
                    {
                        //probe every direction
                        #pragma omp parallel for
                        for (uint64_t i = 0; i < n_directions; i++)
                        {
                            if (collected_chords[i].size() == n_locations) continue;

                            double chordlength = cld.evaluate_singledirection(imagestack, shape, x1, y1, z1, i, color_void);

                            #pragma omp critical
                            {
                                if (chordlength > -1.)
                                    collected_chords[i].push_back(chordlength);
                            }
                        }
                    }
                    else
                        n--;
                }

                //do we have sufficient chords
                minvalid = n_locations;
                for (int64_t i = 0; i < n_directions; i++)
                {
                    if (collected_chords[i].size() < minvalid)
                    {
                        minvalid = collected_chords[i].size();
                    }
                    break;
                }

                cout << "iter " << iter << ": " << minvalid << " chords in worst direction"  << endl;
            }

            //Kernel density estimation and fit an exponential to every direction
            vector<vector<double>> output;
            vector<double> dummy2;
            output.push_back(dummy2);
            for (uint64_t i = 0; i < n_directions; i+=1)
            {
                vector<vector<double>> density_estimate = kde.run_kerneldensityestimation(collected_chords[i], kernel, bandwidth_estimator, stepsize);
                if (i == 0 || density_estimate[0].size() < output[0].size()) output[0] = density_estimate[0];
                output.push_back(density_estimate[1]);
            }
            hdcom.SaveColumnData(output, rootpath, "directionality");

            continue;
        }

        vector<double> all_materialchords;
        vector<double> all_voidchords;

        uint64_t n_slice = shape[0]*shape[1];
        uint64_t n_row = shape[0];
        uint64_t n_voidchords = 0;
        uint64_t n_materialchords = 0;

        uint64_t maxiter = 1e12;

        if (material_only)
        {
            //fill void chords with 0 such that only material is evaluated
            all_voidchords.resize(n_samples, 0.);
            n_voidchords = n_samples+1;
        }

        #pragma omp parallel
        {
            #pragma omp for
            for (uint64_t iteration = 0; iteration < maxiter; iteration++) //dummy iterator. Will break on sufficient chords
            {
                if ((n_voidchords)%10000 == 0)
                {
                    if (n_voidchords != 0){
                        cout << n_voidchords << " valid void chords            \r";
                        cout.flush();
                    }
                }
                if ((n_materialchords)%10000 == 0)
                {
                    if (n_voidchords != 0){
                        cout << n_materialchords << " valid material chords            \r";
                        cout.flush();
                    }
                }
                if ((n_voidchords >= n_samples) && (n_materialchords > n_samples))
                {
                    iteration = maxiter;
                    continue;
                }

                //get a random sample location
                int x = (rand() % (shape[0]-2))+1;
                int y = (rand() % (shape[1]-2))+1;
                int z = (rand() % (shape[2]-2))+1;

                int64_t origin_idx = x + y*n_row + z*n_slice;
                int phasevalue = int(imagestack[origin_idx]);

                if (phasevalue == color_interface)
                {
                    //decide which phase to evaluate
                    phasevalue = rand() % 2;
                    if (phasevalue == 0) phasevalue = color_void;
                    else phasevalue = color_material;
                }
                //don't collect too many chords in case one phase is dominant
                if (((phasevalue == color_void) && (n_voidchords >= n_samples)) || ((phasevalue == color_material) && (n_materialchords >= n_samples)))
                    continue;

                double this_chordlength = -1.;

                if (mode == "average"){
                    vector<double> chords = cld.evaluate_chords_alldirections(imagestack, shape, x, y, z, phasevalue);

                    if (chords.size() > 0)
                    {
                        //get the average chord
                        double sum = std::accumulate(chords.begin(), chords.end(), 0.0);
                        this_chordlength = sum / chords.size();
                    }
                }
                else if (mode == "random"){
                    this_chordlength = cld.evaluate_singlerandomchord(imagestack, shape, x, y, z, phasevalue);
                }
                else if (mode == "minchord"){
                    vector<double> chords = cld.evaluate_chords_alldirections(imagestack, shape, x, y, z, phasevalue);

                    //find the smallest chord
                    if (chords.size() > 0)
                        this_chordlength = *min_element(chords.begin(), chords.end());
                }

                if (this_chordlength == -1.) continue; //invalid chord

                #pragma omp critical
                {
                    if (phasevalue == color_material)
                    {
                        all_materialchords.push_back(this_chordlength);
                        n_materialchords++;
                    }
                    else
                    {
                        all_voidchords.push_back(this_chordlength);
                        n_voidchords++;

                    }
                }
            }
        }
        cout << "\n" << endl;
        //usleep(100);

        //calculate mean and median
        double mean_void = std::accumulate(all_voidchords.begin(), all_voidchords.end(), 0.0)/all_voidchords.size()*voxelsize;
        double mean_material = std::accumulate(all_materialchords.begin(), all_materialchords.end(), 0.0)/all_materialchords.size()*voxelsize;
        double median_void = aux::median(all_voidchords.begin(), all_voidchords.end())*voxelsize;
        double median_material = aux::median(all_materialchords.begin(), all_materialchords.end())*voxelsize;

        //calculate standard deviation
        vector<double> diff1(all_voidchords.size());
        vector<double> diff2(all_materialchords.size());
        transform(all_voidchords.begin(), all_voidchords.end(), diff1.begin(), std::bind2nd(minus<double>(), mean_void/voxelsize));
        transform(all_materialchords.begin(), all_materialchords.end(), diff2.begin(), std::bind2nd(minus<double>(), mean_material/voxelsize));
        double sq_sum1 = inner_product(diff1.begin(), diff1.end(), diff1.begin(), 0.0);
        double sq_sum2 = inner_product(diff2.begin(), diff2.end(), diff2.begin(), 0.0);
        double std_void = sqrt(sq_sum1 / all_voidchords.size())*voxelsize;
        double std_material = sqrt(sq_sum2 / all_materialchords.size())*voxelsize;

        //calculate different means
        double harmonic_mean_void = cld.get_harmonic_mean(all_voidchords)*voxelsize;
        double geometric_mean_void = cld.get_geometric_mean(all_voidchords)*voxelsize;
        double harmonic_mean_material = cld.get_harmonic_mean(all_materialchords)*voxelsize;
        double geometric_mean_material = cld.get_geometric_mean(all_materialchords)*voxelsize;

        //kernel density estimate
        analysis::KDE kde;
        vector<vector<double>> density_estimate_void = kde.run_kerneldensityestimation(all_voidchords, kernel, bandwidth_estimator, stepsize);
        vector<vector<double>> density_estimate_material = kde.run_kerneldensityestimation(all_materialchords, kernel, bandwidth_estimator, stepsize);

        //convert to physical length scale
        std::transform(density_estimate_void[0].begin(), density_estimate_void[0].end(), density_estimate_void[0].begin(), std::bind1st(std::multiplies<double>(), voxelsize));
        std::transform(density_estimate_material[0].begin(), density_estimate_material[0].end(), density_estimate_material[0].begin(), std::bind1st(std::multiplies<double>(), voxelsize));

        //get mode
        uint64_t modepos_void = distance(density_estimate_void[1].begin(), max_element(density_estimate_void[1].begin(), density_estimate_void[1].end()));
        double mode_void = density_estimate_void[0][modepos_void];
        uint64_t modepos_material = distance(density_estimate_material[1].begin(), max_element(density_estimate_material[1].begin(), density_estimate_material[1].end()));
        double mode_material = density_estimate_material[0][modepos_material];

        curvefit::NelderMead nelder;

        vector<double> persistence_params_void, persistence_params_material;
        if (fit_persistencelength){
            //get persistence length
            nelder.fit_persistencelength(density_estimate_void, modepos_void+persistence_offset, mean_void);
            persistence_params_void = nelder.result_params;
            nelder.fit_persistencelength(density_estimate_material, modepos_material+persistence_offset, mean_material);
            persistence_params_material = nelder.result_params;

            //plot persistence length
            vector<double> persistence_fit_void(density_estimate_void[0].size(), 0.);
            vector<double> persistence_fit_material(density_estimate_material[0].size(), 0.);
            for (uint64_t i = modepos_void; i < density_estimate_void[0].size(); i++)
                persistence_fit_void[i] = persistence_params_void[0]*exp(-persistence_params_void[1]*(density_estimate_void[0][i]-persistence_params_void[2]));
            for (uint64_t i = modepos_material; i < density_estimate_material[0].size(); i++)
                persistence_fit_material[i] = persistence_params_material[0]*exp(-persistence_params_material[1]*(density_estimate_material[0][i]-persistence_params_material[2]));
            density_estimate_void.push_back(persistence_fit_void);
            density_estimate_material.push_back(persistence_fit_material);
        }

        vector<double> kgamma_params_void, kgamma_params_material;
        if (fit_kgamma){
            //get kgamma fit
            nelder.fit_kgammafunction(density_estimate_void, mean_void, std_void);
            kgamma_params_void = nelder.result_params;
            nelder.fit_kgammafunction(density_estimate_material, mean_material, std_material);
            kgamma_params_material = nelder.result_params;

            //plot kgamma fit
            vector<double> kgamma_fit_void = nelder.kgamma_function(density_estimate_void[0], kgamma_params_void);
            vector<double> kgamma_fit_material = nelder.kgamma_function(density_estimate_material[0], kgamma_params_material);
            density_estimate_void.push_back(kgamma_fit_void);
            density_estimate_material.push_back(kgamma_fit_material);
        }

        //Output:
        ////////////////////////////////////////////////////////////////////
        cout << "pore space: " << all_voidchords.size() << " chords, mean: " << mean_void << ", mode: " << mode_void << ", median: " << median_void << ", std: " << std_void << endl;
        cout << "material  : " << all_materialchords.size() << " chords, mean: " << mean_material << ", mode: " << mode_material << ", median: " << median_material << ", std: " << std_material << endl;

        hdcom.SaveColumnData(density_estimate_void, tasks[task], "chord_densityestimate_void");
        hdcom.SaveColumnData(density_estimate_material, tasks[task], "chord_densityestimate_material");

        if (export_allchords){
            vector<vector<double>> allchords;
            allchords.push_back(all_voidchords);
            allchords.push_back(all_materialchords);
            hdcom.SaveColumnData(allchords, tasks[task], "allchords");
        }

        if ("append logfile")
        {
            time_t now = time(0);
            ofstream logfile;
            logfile.open(rootpath + "/logfile.txt", fstream::in | fstream::out | fstream::app);
            logfile << ctime(&now);
            logfile << "ran spring8_grainsizeanalysis for " << tasks[task] << ":\n";
            logfile << "-------------------------------------------------------------------------\n";
            logfile << "    n_samples: " << n_samples << "\n";
            logfile << "    mode: " << mode << "\n";
            logfile << "    voxelsize: " << voxelsize << "\n";
            logfile << "    CLD mean (void) [um]: " << mean_void << "\n";
            logfile << "    CLD mean (material) [um]: " << mean_material << "\n";
            logfile << "    CLD mode (void) [um]:  " << mode_void << "\n";
            logfile << "    CLD mode (material) [um]: " << mode_material << "\n";
            logfile << "    CLD median (void) [um]:  " << median_void << "\n";
            logfile << "    CLD median (material) [um]: " << median_material << "\n";
            logfile << "    CLD std (void) [um]:  " << std_void << "\n";
            logfile << "    CLD std (material) [um]: " << std_material << "\n";
            logfile << "    CLD harmonic mean (void) [um]: " << harmonic_mean_void << "\n";
            logfile << "    CLD harmonic mean (material) [um]: " << harmonic_mean_material << "\n";
            logfile << "    CLD geometric mean (void) [um]: " << geometric_mean_void << "\n";
            logfile << "    CLD geometric mean (material) [um]: " << geometric_mean_material << "\n";
            if (fit_persistencelength){
            logfile << "    CLD persistence length (void) [um]: " << persistence_params_void[1] << "    (r2: " << persistence_params_void[3] << ")\n";
            logfile << "    CLD persistence length (material) [um]: " << persistence_params_material[1] << "    (r2: " << persistence_params_material[3] << ")\n";
            }
            if (fit_kgamma){
            logfile << "    CLD kgamma_mu (void) [um]:  " << kgamma_params_void[0] << "\n";
            logfile << "    CLD kgamma_mu (material) [um]: " << kgamma_params_material[0] << "\n";
            logfile << "    CLD kgamma_k (void):  " << kgamma_params_void[1] << "    (r2: " << kgamma_params_void[3] << ")\n";
            logfile << "    CLD kgamma_k (material): " << kgamma_params_material[1] << "    (r2: " << kgamma_params_material[3] << ")\n";
            }
            logfile << "-------------------------------------------------------------------------\n\n";
            logfile.close();
        }

        free(imagestack);
        cout << "##################################" << endl;
    }

    return 0;
}
