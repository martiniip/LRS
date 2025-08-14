/***************************************************************************
                           aco.cpp  -  description
                             -------------------
    begin                : Thu Oct 3 2024
    copyright            : (C) 2024 by Christian Blum
    email                : christian.blum@iiia.csic.es
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <string>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <limits>
#include <random>
#include <chrono>
#include <algorithm>
#include <iomanip>

using namespace std;

struct Run {
    int start;
    int letter;
    int length;
};

struct Solution {
    vector<int> vec;
    int score;
};

struct Option {
    int node;
    double value;
};

// Data structures for the problem data
string orig_input_sequence;
vector<int> input_sequence;
vector<char> alphabet;
int alphabet_size;
int sequence_length;
vector<Run> runs;
int n_of_runs;

// Other data structures
string inputFile;
double computation_time_limit = 100.0;
int trials = 1;

//ACO parameters
int n_of_ants = 50;
double aco_determinism_rate = 0.2;
double tau_max = 0.999;
double tau_min = 0.001;
double learning_rate = 0.1;
bool tuning = false;
bool newheuristic= false;


bool pos_val_compare(const pair<int,double>& i1, const pair<int,double>& i2) {

    return  i1.second > i2.second;
}

inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFile = argv[++iarg];
        else if (strcmp(argv[iarg],"-n_ants")==0) n_of_ants = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-l_rate")==0) learning_rate = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-d_rate")==0) aco_determinism_rate = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-t")==0) computation_time_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-trials")==0) trials = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-tuning")==0) tuning = true;
        else if (strcmp(argv[iarg],"-newheuristic")==0) newheuristic = true;
        iarg++;
    }
}
void generate_aco_solution_new_heuristic(Solution& iSol, vector<double>& pheromone, std::default_random_engine& generator, std::uniform_real_distribution<double>& standard_distribution) {
    vector< pair<int,double> > pos_val_pairs(n_of_runs);
    vector<double> probabilities(n_of_runs);
    for (int i = 0; i < n_of_runs; ++i) {
        pos_val_pairs[i].first = i;
        double score = runs[i].length;
        // Calcular el puntaje adicional por las runs vecinas con la misma letra
        for (int j = 1; j <= 240; ++j) {
            // Verificar vecinos hacia la izquierda (i - j)
            if (i - j >= 0 && runs[i - j].letter == runs[i].letter) {
                score += 0.19 * runs[i - j].length;
            }
            // Verificar vecinos hacia la derecha (i + j)
            if (i + j < n_of_runs && runs[i + j].letter == runs[i].letter) {
                score += 0.19 * runs[i + j].length;
            }
        }
        pos_val_pairs[i].second = score*pheromone[i];
        //pos_val_pairs[i].second = pheromone[i];
    }
    sort(pos_val_pairs.begin(), pos_val_pairs.end(), pos_val_compare);

    for (int i = 0; i < n_of_runs; ++i) probabilities[i] = pos_val_pairs[i].second;

    vector<int> letter_lbs(alphabet_size, -1);
    vector<int> letter_ubs(alphabet_size, -1);
    
    for (int i = 0; i < n_of_runs; ++i) {
        double dec = standard_distribution(generator);
        int current_run;
        if (dec <= aco_determinism_rate) {
            current_run = pos_val_pairs[0].first;
            pos_val_pairs.erase(pos_val_pairs.begin());
            probabilities.erase(probabilities.begin());
        }
        else {
            std::discrete_distribution<> distr(probabilities.begin(), probabilities.end());
            int pos = distr(generator);
            current_run = pos_val_pairs[pos].first;
            pos_val_pairs.erase(pos_val_pairs.begin() + pos);
            probabilities.erase(probabilities.begin() + pos);
  
        }
        bool add = true;
        for (int j = 0; add and (j < alphabet_size); ++j) {
            if (j != runs[current_run].letter) {
                add = false;
                if (not(current_run > letter_lbs[j] and current_run < letter_ubs[j])) {
                    if (not(current_run < letter_lbs[j] and letter_lbs[runs[current_run].letter] > letter_ubs[j])) {
                        if (not(current_run > letter_ubs[j] and letter_ubs[runs[current_run].letter] < letter_lbs[j])) {
                            add = true;
                        }   
                    }
                }
            }
        }
        if (add) {
            (iSol.vec)[current_run] = 1;
            iSol.score += runs[current_run].length;
            if (letter_lbs[runs[current_run].letter] == -1) {
                letter_lbs[runs[current_run].letter] = current_run;
                letter_ubs[runs[current_run].letter] = current_run;
            }
            else {
                if (current_run < letter_lbs[runs[current_run].letter]) letter_lbs[runs[current_run].letter] = current_run;
                else if (current_run > letter_ubs[runs[current_run].letter]) letter_ubs[runs[current_run].letter] = current_run;
            }
        }
    }
}
void generate_aco_solution(Solution& iSol, vector<double>& pheromone, std::default_random_engine& generator, std::uniform_real_distribution<double>& standard_distribution) {

    vector< pair<int,double> > pos_val_pairs(n_of_runs);
    vector<double> probabilities(n_of_runs);
    for (int i = 0; i < n_of_runs; ++i) {
        pos_val_pairs[i].first = i;
        pos_val_pairs[i].second = double(runs[i].length)*pheromone[i];
        //pos_val_pairs[i].second = pheromone[i];
    }
    sort(pos_val_pairs.begin(), pos_val_pairs.end(), pos_val_compare);

    for (int i = 0; i < n_of_runs; ++i) probabilities[i] = pos_val_pairs[i].second;

    vector<int> letter_lbs(alphabet_size, -1);
    vector<int> letter_ubs(alphabet_size, -1);
	
    for (int i = 0; i < n_of_runs; ++i) {
        double dec = standard_distribution(generator);
        int current_run;
        if (dec <= aco_determinism_rate) {
            current_run = pos_val_pairs[0].first;
            pos_val_pairs.erase(pos_val_pairs.begin());
            probabilities.erase(probabilities.begin());
        }
        else {
            std::discrete_distribution<> distr(probabilities.begin(), probabilities.end());
            int pos = distr(generator);
            current_run = pos_val_pairs[pos].first;
            pos_val_pairs.erase(pos_val_pairs.begin() + pos);
            probabilities.erase(probabilities.begin() + pos);
  
        }
        bool add = true;
        for (int j = 0; add and (j < alphabet_size); ++j) {
            if (j != runs[current_run].letter) {
                add = false;
                if (not(current_run > letter_lbs[j] and current_run < letter_ubs[j])) {
                    if (not(current_run < letter_lbs[j] and letter_lbs[runs[current_run].letter] > letter_ubs[j])) {
                        if (not(current_run > letter_ubs[j] and letter_ubs[runs[current_run].letter] < letter_lbs[j])) {
                            add = true;
                        }	
                    }
                }
            }
        }
        if (add) {
            (iSol.vec)[current_run] = 1;
            iSol.score += runs[current_run].length;
            if (letter_lbs[runs[current_run].letter] == -1) {
                letter_lbs[runs[current_run].letter] = current_run;
                letter_ubs[runs[current_run].letter] = current_run;
            }
            else {
                if (current_run < letter_lbs[runs[current_run].letter]) letter_lbs[runs[current_run].letter] = current_run;
                else if (current_run > letter_ubs[runs[current_run].letter]) letter_ubs[runs[current_run].letter] = current_run;
            }
        }
    }
}

double compute_convergence_factor(vector<double>& pheromone) {

    double ret_val = 0.0;
    for (int i = 0; i < n_of_runs; ++i) {
        if ((tau_max - pheromone[i]) > (pheromone[i] - tau_min)) ret_val += tau_max - pheromone[i];
        else ret_val += pheromone[i] - tau_min;
    }
    ret_val = ret_val / (double(n_of_runs) * (tau_max - tau_min));
    ret_val = (ret_val - 0.5) * 2.0;
    return ret_val;
}

/**********
Main function
**********/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);

    std::cout << std::setprecision(3) << std::fixed;

    // initializes the random number generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> standard_distribution(0.0,1.0);

    // reading an instance
    ifstream indata;
    indata.open(inputFile.c_str());
    if(!indata) { // file couldn't be opened
        if (not tuning) cout << "Error: file could not be opened" << endl;
        exit(1);
    }

    map<char,int> rev_mapping;

    indata >> orig_input_sequence;
    sequence_length = int(orig_input_sequence.size());
    //if (not tuning) cout << "Input string: " << orig_input_sequence << " length: " << sequence_length << endl;
    set<char> chars;
    for (int i = 0; i < int(orig_input_sequence.size()); ++i) chars.insert(orig_input_sequence[i]);
    //if (not tuning) cout << "Sigma = {" << *(chars.begin());
    bool startb = true;
    for (set<char>::iterator sit = chars.begin(); sit != chars.end(); ++sit) {
        if (startb) startb = false;
        //else if (not tuning) cout << ", " << *sit;
    }
    //if (not tuning) cout << "}" << endl;
    alphabet_size = int(chars.size());
    //if (not tuning) cout << "|Sigma| = " << alphabet_size << endl;
    alphabet = vector<char>(alphabet_size);
    int j = 0;
    for (set<char>::iterator sit = chars.begin(); sit != chars.end(); ++sit) {
        alphabet[j] = *sit;
        rev_mapping[alphabet[j]] = j;
        ++j;
    }
    for (int j = 0; j < sequence_length; ++j) input_sequence.push_back(rev_mapping[orig_input_sequence[j]]);
    //if (not tuning) cout << "Converted string: ";
    for (int j = 0; j < sequence_length; ++j) //if (not tuning) cout << input_sequence[j];
    //if (not tuning) cout << endl;
    indata.close();
        
    Run r;
    r.letter = -1;
    for (int i = 0; i < sequence_length; ++i) {
        if (input_sequence[i] != r.letter) {
            if (i != 0) {
                r.length = i - r.start;
                runs.push_back(r);
            }
            r.letter = input_sequence[i];
            r.start = i;
        }
    }
    r.length = sequence_length - r.start;
    runs.push_back(r);
    //if (not tuning) cout << "All runs: " << endl;
    for (int i = 0; i < int(runs.size()); ++i) {
        //if (not tuning) cout << "run " << i << " (start: " << runs[i].start << ", letter: " << alphabet[runs[i].letter] << ", length: " << runs[i].length << ")" << endl;
    }
	n_of_runs = int(runs.size());

    vector<double> results(trials, std::numeric_limits<int>::max());
    vector<double> times(trials, 0.0);

    computation_time_limit = sequence_length/5.0;


    // looping over all trials
    for (int trial = 0; trial < trials; trial++) {

        //if (not tuning) cout << "start trial " << trial + 1 << endl;

        // the computation time starts now
        clock_t start = clock();

        // initialization of the pheromone values
        vector<double> pheromone(n_of_runs, 0.5);
        
        Solution bestSol;
        bestSol.score = std::numeric_limits<int>::min();
		bestSol.vec = vector<int>(n_of_runs, 0);

        Solution restartSol;
        restartSol.score = std::numeric_limits<int>::min();
		restartSol.vec = vector<int>(n_of_runs, 0);

        bool global_convergence = false;
        double cf = 0.0;

        clock_t current = clock();
        double ctime = double(current - start) / CLOCKS_PER_SEC;
        bool stop = false;
        while (not stop and ctime < computation_time_limit) {

            // generate n_of_ants solutions
            double iteration_average = 0.0;
            Solution iBestSol;
            iBestSol.score = std::numeric_limits<int>::min();
            for (int na = 0; not stop and na < n_of_ants; ++na) {
                Solution iSol;
				iSol.score = 0;
				iSol.vec = vector<int>(n_of_runs, 0);
                if(newheuristic) generate_aco_solution_new_heuristic(iSol, pheromone, generator, standard_distribution);
                else generate_aco_solution(iSol, pheromone, generator, standard_distribution);
                clock_t current2 = clock();
                ctime = double(current2 - start) / CLOCKS_PER_SEC;
                iteration_average += double(iSol.score);
                if (iSol.score > iBestSol.score) iBestSol = iSol;
                if (iSol.score > restartSol.score) restartSol = iSol;
                if (iSol.score > bestSol.score) {
                    bestSol = iSol;
                    results[trial] = double(bestSol.score);
                    times[trial] = ctime;
                    cout << "score " << bestSol.score << "\ttime " << times[trial] << endl;

                }
                if (ctime >= computation_time_limit) stop = true;
            }
			//cout << "I am here" << endl;
            if (not stop) {
                iteration_average /= double(n_of_ants);

                double i_weight, r_weight, g_weight;
                if (global_convergence) {
                    i_weight = 0.0;
                    r_weight = 0.0;
                    g_weight = 1.0;
                }
                else {
                    if (cf < 0.2) {
                        i_weight = 1.0;
                        r_weight = 0.0;
                        g_weight = 0.0;
                    }
                    else if (cf < 0.5) {
                        i_weight = 2.0/3.0;
                        r_weight = 1.0/3.0;
                        g_weight = 0.0;
                    }
                    else if (cf < 0.8) {
                        i_weight = 1.0/3.0;
                        r_weight = 2.0/3.0;
                        g_weight = 0.0;
                    }
                    else {
                        i_weight = 0.0;
                        r_weight = 1.0;
                        g_weight = 0.0;
                    }
                }

                vector<double> contribution(n_of_runs, 0.0);
				for (int i = 0; i < n_of_runs; ++i) {
					contribution[i] += (i_weight*(iBestSol.vec)[i]);
					contribution[i] += (r_weight*(restartSol.vec)[i]);
					contribution[i] += (g_weight*(bestSol.vec)[i]);
                    pheromone[i] += (learning_rate * (contribution[i] - pheromone[i]));
                    if (pheromone[i] > tau_max) pheromone[i] = tau_max;
                    if (pheromone[i] < tau_min) pheromone[i] = tau_min;
                }

                cf = compute_convergence_factor(pheromone);
                if (cf > 0.99) {
                    cf = 0;
                    if (global_convergence) {
						//cout << "restart" << endl;
                        global_convergence = false;
                        pheromone = vector<double>(n_of_runs, 0.5);
                        restartSol.score = std::numeric_limits<int>::min();
						restartSol.vec = vector<int>(n_of_runs, 0);
                    }
                    else global_convergence = true;
                }
            
                clock_t current_end = clock();
                ctime = double(current_end - start) / CLOCKS_PER_SEC;
            }
        }   

        //if (not tuning) cout << "end trial " << trial + 1 << endl;
    }

    int best_result = std::numeric_limits<int>::min();
    double r_mean = 0.0;
    double g_mean = 0.0;
    for (int i = 0; i < results.size(); i++) {
        r_mean = r_mean + results[i];
        g_mean = g_mean + times[i];
        if (int(results[i]) > best_result) best_result = int(results[i]);
    }
    r_mean = r_mean / ((double)results.size());
    g_mean = g_mean / ((double)times.size());
    double rsd = 0.0;
    double gsd = 0.0;
    for (int i = 0; i < results.size(); i++) {
        rsd = rsd + pow(results[i]-r_mean,2.0);
        gsd = gsd + pow(times[i]-g_mean,2.0);
    }
    rsd = rsd / ((double)(results.size()-1.0));
    if (rsd > 0.0) {
        rsd = sqrt(rsd);
    }
    gsd = gsd / ((double)(times.size()-1.0));
    if (gsd > 0.0) {
        gsd = sqrt(gsd);
    }
    //if (not tuning) cout << best_result << "\t" << r_mean << "\t" << rsd << "\t" << g_mean << "\t" << gsd<< endl;
    //else cout << (-1)*results[0] << endl;
    cout<<best_result <<endl;
}
