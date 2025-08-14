/***************************************************************************
                         brkga.cpp  -  description
                             -------------------
    begin                : Thu Sep 26 2024
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


struct Individual {
    vector<double> vec;
    set<int> selected_runs;
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
vector<double> LLM;



// Other data structures
string inputFile;
string model;
double computation_time_limit = 100.0;
int trials = 1;


// BRKGA parameters
int population_size = 30;
double elite_proportion = 0.15; // normally between 0.1 and 0.25
double mutant_proportion = 0.20; // normally between 0.1 and 0.3
double elite_inheritance_probability = 0.7; // normally greater than 0.5 and <= 0.8
double threshold = 0.7;
bool tuning = false;
int id=1;
bool stnweb = false;




bool individual_compare(const Individual& i1, const Individual& i2) {


    return  i1.score > i2.score;
}


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
        else if (strcmp(argv[iarg],"-p")==0) population_size = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pe")==0) elite_proportion = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pm")==0) mutant_proportion = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-rhoe")==0) elite_inheritance_probability = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-t")==0) computation_time_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-th")==0) threshold = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-trials")==0) trials = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-tuning")==0) tuning = true;
        else if (strcmp(argv[iarg],"-model")==0) model = argv[++iarg];
        else if (strcmp(argv[iarg],"-id")==0) id = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-stnweb")==0) stnweb = true;
        iarg++;
    }
}


int produce_random_integer(int max, default_random_engine& generator, uniform_real_distribution<double>& standard_distribution) {


    double rnd = standard_distribution(generator);
    int num = int(double(max) * rnd);
    if (num == max) num = num - 1;
    return num;
}


int get_random_element(const set<int>& s, default_random_engine& generator, uniform_real_distribution<double>& standard_distribution) {


    double r = produce_random_integer(int(s.size()), generator, standard_distribution);
    set<int>::iterator it = s.begin();
    advance(it, r);
    return *it; 
}


void evaluate(Individual& ind) {


    vector< pair<int,double> > pos_val_pairs(n_of_runs);
    for (int i = 0; i < n_of_runs; ++i) {
        pos_val_pairs[i].first = i;
        pos_val_pairs[i].second = (ind.vec)[i]*LLM[i];
    }
    sort(pos_val_pairs.begin(), pos_val_pairs.end(), pos_val_compare);


    vector<int> letter_lbs(alphabet_size, -1);
    vector<int> letter_ubs(alphabet_size, -1);
    (ind.selected_runs).clear();
    ind.score = 0;
    
    for (int i = 0; i < n_of_runs; ++i) {
        int current_run = pos_val_pairs[i].first;
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
            (ind.selected_runs).insert(current_run);
            ind.score += runs[current_run].length;
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
    //cout << "sol.score: " << ind.score << endl;
}


void generate_random_solution(Individual& ind, default_random_engine& generator, uniform_real_distribution<double>& standard_distribution) {


    ind.vec = vector<double>(n_of_runs);
    for (int i = 0; i < n_of_runs; ++i) (ind.vec)[i] = standard_distribution(generator);
    evaluate(ind);
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
    cout << "Input string: " << orig_input_sequence << " length: " << sequence_length << endl;
    set<char> chars;
    for (int i = 0; i < int(orig_input_sequence.size()); ++i) chars.insert(orig_input_sequence[i]);
    cout << "Sigma = {" << *(chars.begin());
    bool startb = true;
    for (set<char>::iterator sit = chars.begin(); sit != chars.end(); ++sit) {
        if (startb) startb = false;
        else cout << ", " << *sit;
    }
    cout << "}" << endl;
    alphabet_size = int(chars.size());
    cout << "|Sigma| = " << alphabet_size << endl;
    alphabet = vector<char>(alphabet_size);
    int j = 0;
    for (set<char>::iterator sit = chars.begin(); sit != chars.end(); ++sit) {
        alphabet[j] = *sit;
        rev_mapping[alphabet[j]] = j;
        ++j;
    }
    for (int j = 0; j < sequence_length; ++j) input_sequence.push_back(rev_mapping[orig_input_sequence[j]]);
    cout << "Converted string: ";
    for (int j = 0; j < sequence_length; ++j) cout << input_sequence[j];
    cout << endl;
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
    cout << "All runs: " << endl;
    for (int i = 0; i < int(runs.size()); ++i) {
        cout << "run " << i << " (start: " << runs[i].start << ", letter: " << alphabet[runs[i].letter] << ", length: " << runs[i].length << ")" << endl;
    }
    n_of_runs = int(runs.size());


    //Read probabilities from folder LLM_probabilities
    size_t lastSlash = inputFile.find_last_of("/\\");
    string filename = (lastSlash == std::string::npos) ? inputFile : inputFile.substr(lastSlash + 1);

    string probpath = "../LLM_probabilities/" + model + "/" + filename;
    //cout<<probpath<<endl;
    ifstream infile(probpath);
    if (!infile.is_open()) {
        probpath = "../../LLM_probabilities/" + model + "/" + filename;
        //cout<<probpath<<endl;
        ifstream infile(probpath);
        if (!infile.is_open()) {
            cerr << "Error al abrir el archivo de probabilidades: " << probpath << endl;
            return 1;
        }
    }

    double value;
    while (infile >> value) {
        //cout<<value<<endl;

        LLM.push_back(value);
    }

    infile.close();

    cout << "Se leyeron " << LLM.size() << " valores del archivo." <<endl;

    vector<double> results(trials, std::numeric_limits<int>::max());
    vector<double> times(trials, 0.0);


    // looping over all trials
    for (int trial = 0; trial < trials; trial++) {


        if (not tuning) cout << "start trial " << trial + 1 << endl;


        // the computation time starts now
        clock_t start = clock();


        int best_score = -1;


        int n_elites = int(double(population_size)*elite_proportion);
        if (n_elites < 1) n_elites = 1;


        int n_mutants = int(double(population_size)*mutant_proportion);
        if (n_mutants < 1) n_mutants = 1;


        int n_offspring = population_size - n_elites - n_mutants;
        if (n_offspring < 1) {
            if (not tuning) cout << "OHOHOH: wrong parameter settings" << endl;
            exit(0);
        }
        computation_time_limit= sequence_length/5;


        double ctime;
        vector<Individual> population(population_size);
        for (int pi = 0; pi < population_size; ++pi) {
            generate_random_solution(population[pi], generator, standard_distribution);
            if (population[pi].score > best_score) {
                    if(stnweb){

                        //formato para stnweb: id,score,binary_solution
                        cout << id << "," << population[pi].score << ",";
                        
                        string binary_solution(n_of_runs, '0');
                        for (int run_index : population[pi].selected_runs) {
                            if (run_index >= 0 && run_index < n_of_runs) { 
                                binary_solution[run_index] = '1';
                            }
                        }
                        
                        cout << binary_solution << endl;

                    best_score = population[pi].score;
                    clock_t current2 = clock();
                    ctime = double(current2 - start) / CLOCKS_PER_SEC;
                    results[trial] = double(best_score);
                    times[trial] = ctime;


                }else{
                    for (set<int>::iterator sit = (population[pi].selected_runs).begin(); sit != (population[pi].selected_runs).end(); ++sit) {
                        for (int j = 0; j < runs[*sit].length; ++j) cout << alphabet[runs[*sit].letter];
                    }
                    cout << endl;
                    best_score = population[pi].score;
                    clock_t current2 = clock();
                    ctime = double(current2 - start) / CLOCKS_PER_SEC;
                    results[trial] = double(best_score);
                    times[trial] = ctime;
                    if (not tuning) cout << "score " << best_score << "\ttime " << ctime << endl;
                }
            }
        }


        clock_t current3 = clock();
        ctime = double(current3 - start) / CLOCKS_PER_SEC;
        while (ctime < computation_time_limit) {
            sort(population.begin(), population.end(), individual_compare);
            vector<Individual> new_population(population_size);
            for (int ic = 0; ic < n_elites; ++ic) {
                new_population[ic].vec = population[ic].vec;
                new_population[ic].score = population[ic].score;
                new_population[ic].selected_runs = population[ic].selected_runs;
            }
            for (int ic = 0; ic < n_mutants; ++ic) {
                generate_random_solution(new_population[n_elites + ic], generator, standard_distribution);


                if (new_population[n_elites + ic].score > best_score) {
                    // Escribir el formato stnweb: 1,score,solution_binary
                    if (stnweb) {

                            cout << id << "," << new_population[n_elites + ic].score << ",";
                            string binary_solution(n_of_runs, '0');
                            for (int run_index : new_population[n_elites + ic].selected_runs) {
                                if (run_index >= 0 && run_index < n_of_runs) {
                                    binary_solution[run_index] = '1';
                                }
                            }
                            
                            cout << binary_solution << endl;

                        best_score = new_population[n_elites + ic].score;
                        clock_t current4 = clock();
                        ctime = double(current4 - start) / CLOCKS_PER_SEC;
                        results[trial] = double(best_score);
                        times[trial] = ctime;
                    }else{
                        for (set<int>::iterator sit = (new_population[n_elites + ic].selected_runs).begin(); sit != (new_population[n_elites + ic].selected_runs).end(); ++sit) {
                            for (int j = 0; j < runs[*sit].length; ++j) cout << alphabet[runs[*sit].letter];
                        }
                        cout << endl;
                        best_score = new_population[n_elites + ic].score;
                        clock_t current4 = clock();
                        ctime = double(current4 - start) / CLOCKS_PER_SEC;
                        results[trial] = double(best_score);
                        times[trial] = ctime;
                        if (not tuning) cout << "score " << best_score << "\ttime " << ctime << endl;
                    }
                }            
            }
            for (int ic = 0; ic < n_offspring; ++ic) {
                int first_parent = produce_random_integer(n_elites, generator, standard_distribution);
                int second_parent = n_elites + produce_random_integer(population_size - n_elites, generator, standard_distribution);
                new_population[n_elites + n_mutants + ic].vec = vector<double>(n_of_runs);
                for (int i = 0; i < n_of_runs; ++i) {
                    double rnum = standard_distribution(generator);
                    if (rnum <= elite_inheritance_probability) (new_population[n_elites + n_mutants + ic].vec)[i] = (population[first_parent].vec)[i];
                    else (new_population[n_elites + n_mutants + ic].vec)[i] = (population[second_parent].vec)[i];
                }
                evaluate(new_population[n_elites + n_mutants + ic]);


                if (new_population[n_elites + n_mutants + ic].score > best_score) {
                    if (stnweb) {

                            cout << id <<","<< new_population[n_elites + n_mutants + ic].score << ",";
                            string binary_solution(n_of_runs, '0');
                            for (int run_idx : new_population[n_elites + n_mutants + ic].selected_runs) {
                                if (run_idx >= 0 && run_idx < n_of_runs) {
                                    binary_solution[run_idx] = '1';
                                }
                            }
                            
                            cout << binary_solution << endl;

                        best_score = new_population[n_elites + n_mutants + ic].score;
                        clock_t current4 = clock();
                        ctime = double(current4 - start) / CLOCKS_PER_SEC;
                        results[trial] = double(best_score);
                        times[trial] = ctime;
                    } else {
                        for (set<int>::iterator sit = (new_population[n_elites + n_mutants + ic].selected_runs).begin(); sit != (new_population[n_elites + n_mutants + ic].selected_runs).end(); ++sit) {
                            for (int j = 0; j < runs[*sit].length; ++j) cout << alphabet[runs[*sit].letter];
                        }
                        cout << endl;
                        best_score = new_population[n_elites + n_mutants + ic].score;
                        clock_t current4 = clock();
                        ctime = double(current4 - start) / CLOCKS_PER_SEC;
                        results[trial] = double(best_score);
                        times[trial] = ctime;
                        if (not tuning) cout << "score " << best_score << "\ttime " << ctime << endl;
                    }
                }
            }
            population.clear();
            population = new_population;
            clock_t current5 = clock();
            ctime = double(current5 - start) / CLOCKS_PER_SEC;
        }   


        if (not tuning) cout << "end trial " << trial + 1 << endl;
    }


    int best_result = -1;
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
    if (not tuning) cout << best_result << "\t" << r_mean << "\t" << rsd << "\t" << g_mean << "\t" << gsd<< endl;
    else cout << results[0] << endl;
}


