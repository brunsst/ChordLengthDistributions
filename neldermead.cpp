#include "neldermead.h"

namespace curvefit
{
    /************ Fitting ****************/
    bool NelderMead::run_unbounded()
    {
        _initialize_simplexes_type1();

        uint64_t n_params = guess.size();
        uint64_t n_simplexes = n_params+1;
        uint64_t best_idx = -1;
        double cost_best = std::numeric_limits<double>::infinity();
        uint64_t iteration;
        bool successfull = false;

        std::vector<double> simplex_cost(n_simplexes, 0.);


        for (iteration = 0; iteration < maxiter; iteration++)
        {
            for (uint64_t j = 0; j < n_simplexes; j++)
                simplex_cost[j] = costfunction1(simplexes[j]);

            //determine position of worsed, second worsed and best simplex
            best_idx = std::min_element(simplex_cost.begin(), simplex_cost.end())-simplex_cost.begin();
            cost_best = simplex_cost[best_idx];
            if (cost_best <= target_cost)
            {
                successfull = true;
                break;
            }
            uint64_t worsed_idx = std::max_element(simplex_cost.begin(), simplex_cost.end())-simplex_cost.begin();
            uint64_t secondworsed_idx = _getsecondworsedsimplex(simplex_cost, worsed_idx);

            //determine centroid of the n best simplexes
            std::vector<double> centroid(n_params, 0.);
            for(uint64_t n = 0; n < n_params; n++)
            {
                for(uint64_t j = 0; j < n_simplexes; j++)
                    centroid[n] += simplexes[j][n];
                centroid[n] -= simplexes[worsed_idx][n];
                centroid[n] /= n_params;
                //std::cout << simplexes[best_idx][n] << "  ";
            }

            //reflected worsed point
            std::vector<double> reflected(n_params, 0.);
            for(uint64_t n = 0; n < n_params; n++)
                reflected[n] = centroid[n] + alpha*(centroid[n]-simplexes[worsed_idx][n]);

            double cost_2ndworsed = simplex_cost[secondworsed_idx];
            double cost_reflected = costfunction1(reflected);

            //std::cout << " <- "<< cost_best << std::endl;

            if ((cost_best <= cost_reflected) && (cost_reflected < cost_2ndworsed))
            {
                //reflect
                simplexes[worsed_idx] = reflected;
            }
            else if (cost_reflected < cost_best)
            {
                //expand
                std::vector<double> expansion(n_params, 0.);
                for(uint64_t n = 0; n < n_params; n++)
                    expansion[n] = reflected[n] + beta*(reflected[n]- centroid[n]);

                double cost_expansion = costfunction1(expansion);
                if (cost_expansion < cost_reflected)
                    simplexes[worsed_idx] = expansion;
                else
                    simplexes[worsed_idx] = reflected;
            }
            else //if (cost_reflected >= cost_2ndworsed)
            {
                double cost_worsed = simplex_cost[worsed_idx];

                //contraction
                std::vector<double> contraction(n_params, 0.);
                for(uint64_t n = 0; n < n_params; n++)
                    contraction[n] = centroid[n] + gamma*(simplexes[worsed_idx][n]-centroid[n]);


                double cost_contraction = costfunction1(contraction);

                if(cost_contraction < cost_worsed)
                    simplexes[worsed_idx] = contraction;
                else
                {
                    //reduce
                    for(uint64_t j = 0; j < n_simplexes; j++)
                    {
                        if (j == best_idx)
                            continue;
                        else
                        {
                            for(uint64_t n = 0; n < n_params; n++)
                                simplexes[j][n] = simplexes[best_idx][n] + sigma*(simplexes[j][n]-simplexes[best_idx][n]);
                        }
                    }
                }
            }
            //if((iteration+1) % 10 == 0)
            //    std::cin.get();
        }

        result_params.assign(n_params, 0.);
        for(uint64_t n = 0; n < n_params; n++)
            result_params[n] = simplexes[best_idx][n];
        result_iterations = iteration;
        result_cost = cost_best;
        return successfull;
    }
    bool NelderMead::run_bounded()
    {
        _initialize_simplexes_type1();

        uint64_t n_params = guess.size();
        uint64_t n_simplexes = n_params+1;
        uint64_t best_idx = -1;
        double cost_best = std::numeric_limits<double>::infinity();
        uint64_t iteration;
        bool successfull = false;

        std::vector<double> simplex_cost(n_simplexes, 0.);


        for (iteration = 0; iteration < maxiter; iteration++)
        {
            for (uint64_t j = 0; j < n_simplexes; j++)
                simplex_cost[j] = costfunction1(simplexes[j]);

            //determine position of worsed, second worsed and best simplex
            best_idx = std::min_element(simplex_cost.begin(), simplex_cost.end())-simplex_cost.begin();
            cost_best = simplex_cost[best_idx];
            if (cost_best <= target_cost)
            {
                successfull = true;
                break;
            }
            uint64_t worsed_idx = std::max_element(simplex_cost.begin(), simplex_cost.end())-simplex_cost.begin();
            uint64_t secondworsed_idx = _getsecondworsedsimplex(simplex_cost, worsed_idx);

            //determine centroid of the n best simplexes
            std::vector<double> centroid(n_params, 0.);
            for(uint64_t n = 0; n < n_params; n++)
            {
                for(uint64_t j = 0; j < n_simplexes; j++)
                    centroid[n] += simplexes[j][n];
                centroid[n] -= simplexes[worsed_idx][n];
                centroid[n] /= n_params;
                //std::cout << simplexes[best_idx][n] << "  ";
            }

            //reflected worsed point
            std::vector<double> reflected(n_params, 0.);
            for(uint64_t n = 0; n < n_params; n++)
            {
                reflected[n] = centroid[n] + alpha*(centroid[n]-simplexes[worsed_idx][n]);
                if (reflected[n] < lower_bounds[n]) reflected[n] = lower_bounds[n];
                if (reflected[n] > upper_bounds[n]) reflected[n] = upper_bounds[n];
            }

            double cost_2ndworsed = simplex_cost[secondworsed_idx];
            double cost_reflected = costfunction1(reflected);

            //std::cout << " <- "<< cost_best << std::endl;

            if ((cost_best <= cost_reflected) && (cost_reflected < cost_2ndworsed))
            {
                //reflect
                simplexes[worsed_idx] = reflected;
            }
            else if (cost_reflected < cost_best)
            {
                //expand
                std::vector<double> expansion(n_params, 0.);
                for(uint64_t n = 0; n < n_params; n++)
                {
                    expansion[n] = reflected[n] + beta*(reflected[n]- centroid[n]);
                    if (expansion[n] < lower_bounds[n]) expansion[n] = lower_bounds[n];
                    if (expansion[n] > upper_bounds[n]) expansion[n] = upper_bounds[n];
                }

                double cost_expansion = costfunction1(expansion);
                if (cost_expansion < cost_reflected)
                    simplexes[worsed_idx] = expansion;
                else
                    simplexes[worsed_idx] = reflected;
            }
            else //if (cost_reflected >= cost_2ndworsed)
            {
                double cost_worsed = simplex_cost[worsed_idx];

                //contraction
                std::vector<double> contraction(n_params, 0.);
                for(uint64_t n = 0; n < n_params; n++)
                {
                    contraction[n] = centroid[n] + gamma*(simplexes[worsed_idx][n]-centroid[n]);
                    if (contraction[n] < lower_bounds[n]) contraction[n] = lower_bounds[n];
                    if (contraction[n] > upper_bounds[n]) contraction[n] = upper_bounds[n];
                }


                double cost_contraction = costfunction1(contraction);

                if(cost_contraction < cost_worsed)
                    simplexes[worsed_idx] = contraction;
                else
                {
                    //reduce
                    for(uint64_t j = 0; j < n_simplexes; j++)
                    {
                        if (j == best_idx)
                            continue;
                        else
                        {
                            for(uint64_t n = 0; n < n_params; n++)
                            {
                                simplexes[j][n] = simplexes[best_idx][n] + sigma*(simplexes[j][n]-simplexes[best_idx][n]);
                                if (simplexes[j][n] < lower_bounds[n]) simplexes[j][n] = lower_bounds[n];
                                if (simplexes[j][n] > upper_bounds[n]) simplexes[j][n] = upper_bounds[n];
                            }
                        }
                    }
                }
            }
            //if((iteration+1) % 10 == 0)
            //    std::cin.get();
        }

        result_params.assign(n_params, 0.);
        for(uint64_t n = 0; n < n_params; n++)
            result_params[n] = simplexes[best_idx][n];
        result_iterations = iteration;
        result_cost = cost_best;
        return successfull;
    }

    /********* simplex initialization **********/
    void NelderMead::_initialize_simplexes_type1()
    {
        simplexes.clear();

        //x0:
        std::vector<double> simplex(guess.size(), 0.);
        for (uint64_t j = 0; j < guess.size(); j++)
            simplex[j] = guess[j];
        simplexes.push_back(simplex);

        for(uint64_t j = 0; j < guess.size(); j++)
        {
            if(guess[j] != 0.)
                simplex[j] = guess[j] + 0.05;
            else
                simplex[j] = guess[j] + 0.00025;

            simplexes.push_back(simplex);
            simplex[j] = guess[j];
        }
        return;
    }

    /********* auxiliary ****************/
    uint64_t NelderMead::_getsecondworsedsimplex(std::vector<double> &simplex_cost, uint64_t &worsed_idx)
    {
        uint64_t secondworsed_idx;
        if (worsed_idx == 0)
            secondworsed_idx = std::max_element(simplex_cost.begin()+1, simplex_cost.end())-simplex_cost.begin();
        else if (worsed_idx == (simplexes.size() - 1))
            secondworsed_idx = std::max_element(simplex_cost.begin(), simplex_cost.end()-1)-simplex_cost.begin();
        else
        {
            uint64_t tmp = std::max_element(simplex_cost.begin(), simplex_cost.begin()+worsed_idx)-simplex_cost.begin();
            uint64_t tmp2 =std::max_element(simplex_cost.begin()+worsed_idx+1, simplex_cost.end())-simplex_cost.begin();
            if (simplex_cost[tmp] >= simplex_cost[tmp2])
                secondworsed_idx = tmp;
            else
                secondworsed_idx = tmp2;
        }
        return secondworsed_idx;
    }
}

