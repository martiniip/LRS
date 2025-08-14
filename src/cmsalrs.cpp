/***************************************************************************
                          cplex.cpp  -  description
                             -------------------
    begin                : Thu Dec 16 2021
    copyright            : (C) 2021 by Christian Blum
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
  // the following "include" is necessary for the correct working/compilation of CPLEX.
#include <ilcplex/ilocplex.h>

using namespace std;

ILOSTLBEGIN

struct Run {
    int start;
    int letter;
    int length;
};



struct sub_instance_LRS{
    vector <int> runs;
    vector <int> runs_age;

};

// Parámetros CMSA
struct CMSAparams {
    double computation_time_limit = 100.0;
    double cplex_time_limit = 10.0;
    double determinism_rate = 0.8;
    int n_solutions = 10;
    int age_limit = 5;
    int candidate_list_size = 10;
    int neighborhood_size = 240;
    double neighborhood_weight = 0.19;
    bool warm_start = false;
    bool cplex_abort = false;
    bool heuristic_emphasis = false;
    int max_iter_local_search= 10;
};

CMSAparams params;

//Data structures for BRKGA

struct Individual {
    vector<double> vec;
    set<int> selected_runs;
    int score;
};

struct Option {
    int node;
    double value;
};


// Other data structures
string inputFile;
double computation_time_limit = 100.0;
//int trials = 1;
clock_t start_time;


// Data structures for the problem data
string orig_input_sequence;
vector<int> input_sequence;
vector<char> alphabet;
int alphabet_size;
int sequence_length;
vector<Run> runs;
int n_of_runs;
int stn_run_number=1; 

// BRKGA parameters
int population_size = 30;
double elite_proportion = 0.15; // normally between 0.1 and 0.25
double mutant_proportion = 0.20; // normally between 0.1 and 0.3
double elite_inheritance_probability = 0.7; // normally greater than 0.5 and <= 0.8
//double threshold = 0.7;
bool tuning = false;

//L_CMSA parameters
double best_BRKGA_fitness=-1;
double best_LCMSA_quality=0;
double best_CPLEX_quality=0;
//double computational_time_limit = 3200.0;
int max_ga_iterations = 1000;
sub_instance_LRS sub_instance;
int individuals_selected_from_LC_for_merge = 2;
int max_age = 5;
int number_CPLEX_solution_copies_injected = 2;
double proportion_CPLEX_solutions_copies_injected = 0.5;

vector<int> best_LCMSA_individual;
vector<int> best_CPLEX_individual;
vector<int> current_CPLEX_individual;
int current_CPLEX_quality=0;

// Variables de estado
int output_type = 0; // Tipo de salida (0=normal, 1=STN)



// vector for keeping all the names of the input files
vector<string> inputFiles;

// time limit for CPLEX (can be supplied to the algorithm via the -t comand line parameter)
double cplex_time_limit = 5;



inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}


bool individual_compare(const Individual& i1, const Individual& i2) {

    return  i1.score > i2.score;
}


bool pos_val_compare(const pair<int,double>& i1, const pair<int,double>& i2) {

    return  i1.second > i2.second;
}






void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFiles.push_back(argv[++iarg]);
        else if (strcmp(argv[iarg],"-n_sol")==0) params.n_solutions = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-d")==0) params.determinism_rate = atof(argv[++iarg]);
        //else if (strcmp(argv[iarg],"-th")==0) threshold = atof(argv[++iarg]);
        //else if (strcmp(argv[iarg],"-trials")==0) trials = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-tuning")==0) tuning = true;
        else if (strcmp(argv[iarg],"-cpl_t")==0) params.cplex_time_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-t")==0) computation_time_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-iter")==0) params.max_iter_local_search = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-na")==0) individuals_selected_from_LC_for_merge = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-m_age")==0) params.age_limit = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-output_type")==0) output_type = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-h")==0) {
            cout << "Usage: " << argv[0] << " -i <input file> -p <population size> -pe <elite proportion> -pm <mutant proportion> -rhoe <elite inheritance probability> -th <threshold> -t <computation time limit> -cpl_t <cplex time limit> -iter <max ga iterations> -na <individuals selected from LC for merge> -m_age <max age> -n_c <number CPLEX solution copies injected>" << endl;
            exit(0);
        }


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


int evaluate(Individual& ind) {

    vector< pair<int,double> > pos_val_pairs(n_of_runs);
    for (int i = 0; i < n_of_runs; ++i) {
        pos_val_pairs[i].first = i;
        pos_val_pairs[i].second = (ind.vec)[i];
    }
    sort(pos_val_pairs.begin(), pos_val_pairs.end(), pos_val_compare);
    /*
    for (int i = 0; i < n_of_runs; ++i) {
        cout << " (" << pos_val_pairs[i].first << "," << pos_val_pairs[i].second << ")";
    }
    cout << endl;
    */

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
    return ind.score;
}

void generate_random_solution(Individual& ind, default_random_engine& generator, uniform_real_distribution<double>& standard_distribution) {

    ind.vec = vector<double>(n_of_runs);
    for (int i = 0; i < n_of_runs; ++i) (ind.vec)[i] = standard_distribution(generator);
    evaluate(ind);
}

void generate_greedy_random(Individual& ind, default_random_engine& generator, 
                          uniform_real_distribution<double>& standard_distribution, 
                          int k, double ponderacion, double determinism_level) {
    ind.vec = vector<double>(n_of_runs);
    uniform_real_distribution<double> dis_random(0.0, 1.0);

    for (int i = 0; i < n_of_runs; ++i) {
        if (dis_random(generator) > determinism_level) {
            ind.vec[i] = standard_distribution(generator);
        } else {
            int current_run = i;
            double score = runs[current_run].length;

            for (int j = 1; j <= k; ++j) {
                if (current_run - j >= 0 && runs[current_run - j].letter == runs[current_run].letter) {
                    score += ponderacion * runs[current_run - j].length;
                }
                if (current_run + j < n_of_runs && runs[current_run + j].letter == runs[current_run].letter) {
                    score += ponderacion * runs[current_run + j].length;
                }
            }

            ind.vec[i] = score;
        }
    }

    evaluate(ind);
}





/*


ILOSOLVECALLBACK4(loggingCallback,
    clock_t&, start,
    double&, time_stamp,
    double&, result,
    double&, gap) {

    if (hasIncumbent()) {
        IloNum nv = getIncumbentObjValue();
        clock_t current = clock();
        double newTime = double(current - start) / CLOCKS_PER_SEC;
        double newGap = 100.0 * getMIPRelativeGap();
        if (result < double(nv)) {
            cout << "value " << nv << "\ttime " << newTime << "\tgap " << newGap << endl;
            result = double(nv);
            time_stamp = newTime;
            gap = newGap;
        }
    }
}

*/

/* 

//ESTE ES EL QUE REALMENTE SE USA
ILOSOLVECALLBACK4(loggingCallback,
    clock_t&, start,
    double&, time_stamp,
    double&, result,
    double&, gap) {

    if (hasIncumbent()) {
        IloNum nv = getIncumbentObjValue();
        clock_t current = clock();
        double newTime = double(current - start) / CLOCKS_PER_SEC;
        double newGap = 100.0 * getMIPRelativeGap();
        if (best_LCMSA_quality < double(nv)) {
            if(output_type==0) cout << "-CP- value " << nv << "\ttime " << newTime << "\tgap " << newGap << endl;
            result = double(nv);
            time_stamp = newTime;
            gap = newGap;
            best_CPLEX_quality = double(nv);
            best_LCMSA_quality = double(nv);
            

        }
    }
}
*/
void print_solution_in_letters(const Individual& ind) {
    for (set<int>::iterator sit = (ind.selected_runs).begin(); sit != (ind.selected_runs).end(); ++sit) {
                    for (int j = 0; j < runs[*sit].length; ++j) cout << alphabet[runs[*sit].letter];
    }
    cout<<endl;
}



//experimental
ILOSOLVECALLBACK4(loggingCallback,
    clock_t&, start,
    IloIntVarArray, x,
    double&, result,
    double&, gap) {

    if (hasIncumbent()) {
        IloNum nv = getIncumbentObjValue();
        clock_t current = clock();
        double newTime = double(current - start) / CLOCKS_PER_SEC;
        double newGap = 100.0 * getMIPRelativeGap();
        if (best_LCMSA_quality < double(nv)) {
            if(tuning == false) cout << "-*CP- value " << nv << "\ttime " << newTime << "\tgap " << newGap << endl;
            else cout << "score " << nv << "\ttime " << newTime << endl;
            result = double(nv);
            //time_stamp = newTime;
            gap = newGap;
            best_CPLEX_quality = double(nv);
            best_LCMSA_quality = double(nv); 
            if(output_type==1) {
                    cout<<"*CP,"<<stn_run_number<<",";
                    for (int i = 0; i < n_of_runs; ++i) {
                        IloNum xval = getIncumbentValue(x[i]); 
                        if (xval > 0.8) cout<<"1"; else cout<<"0";
                        }
                    cout<<","<<-1 * best_LCMSA_quality<<endl;
                    
                }
            

        }
    }
}
class LazyConstraintCallback : public IloCplex::LazyConstraintCallbackI {
    IloIntVarArray x;
    vector<Run> runs;
public:
    LazyConstraintCallback(IloEnv env, IloIntVarArray _x, vector<Run> _runs)
        : IloCplex::LazyConstraintCallbackI(env), x(_x), runs(_runs) {}

    void main() override {
        for (int i = 0; i < runs.size() - 1; ++i) {
            for (int j = i + 1; j < runs.size(); ++j) {
                if (runs[i].letter == runs[j].letter) {
                    double xi = getValue(x[i]);
                    double xj = getValue(x[j]);
                    if (xi > 0.5 && xj > 0.5) {
                        IloExpr expr(getEnv());
                        expr += (j - i) * x[i];
                        expr += (j - i) * x[j];
                        for (int l = i + 1; l < j; ++l) {
                            if (runs[l].letter != runs[i].letter) expr += x[l];
                        }
                        add(expr <= 2 * (j - i)).end();
                    }
                }
            }
        }
    }

    IloCplex::CallbackI* duplicateCallback() const override {
        return new (getEnv()) LazyConstraintCallback(*this);
    }
};



vector<int> run_cplex(clock_t& start, vector<double>& results, vector<double>& times, vector<double>& gaps, int& na) {

    n_of_runs = int(runs.size());

    IloEnv env;
    env.setOut(env.getNullStream());
    try {
        IloModel model(env);

        IloIntVarArray x(env, n_of_runs, 0, 1);
        
        IloExpr obj(env);
        for (int i = 0; i < n_of_runs; ++i) obj += runs[i].length * x[i];
        model.add(IloMaximize(env, obj));

        // Subinstance constraints
        for (int i = 0; i < n_of_runs; i++) {
            if (sub_instance.runs[i] == 0) {
                IloExpr expr(env);
                expr += x[i];
                model.add(expr == 0);
                expr.end();
            }
        }

        IloCplex cpl(model);

        // Parámetros
        cpl.setParam(IloCplex::TiLim, params.cplex_time_limit);
        cpl.setParam(IloCplex::EpGap, 0.01);
        cpl.setParam(IloCplex::EpAGap, 0.1);
        cpl.setParam(IloCplex::Threads, 1);
        cpl.setParam(IloCplex::MIPEmphasis, 1);
        cpl.setParam(IloCplex::PreInd, 0);

        // Callbacks
        cpl.use(loggingCallback(env, start, x, best_LCMSA_quality, gaps[na]));
        cpl.use(new (env) LazyConstraintCallback(env, x, runs)); 

        cpl.setWarning(env.getNullStream());

        cpl.solve();

        if (cpl.getStatus() == IloAlgorithm::Optimal || cpl.getStatus() == IloAlgorithm::Feasible) {
            clock_t current = clock();
            double newTime = double(current - start) / CLOCKS_PER_SEC;
            double lastVal = double(cpl.getObjValue());
            double lastGap = 100.0 * cpl.getMIPRelativeGap();

            if (lastVal > best_LCMSA_quality) {
                gaps[na] = lastGap;
                best_CPLEX_quality = lastVal;
                best_LCMSA_quality = lastVal;

                if (!tuning) cout << "-CP- value " << lastVal << "\ttime " << newTime << "\tgap " << lastGap << endl;
                else cout << "score " << lastVal << "\ttime " << newTime << endl;

                if (output_type == 1) {
                    cout << "CP," << stn_run_number << ",";
                    for (int i = 0; i < n_of_runs; ++i) {
                        IloNum xval = cpl.getValue(x[i]);
                        cout << (xval > 0.8 ? "1" : "0");
                    }
                    cout << "," << -1 * best_LCMSA_quality << endl;
                }
            }

            current_CPLEX_quality = lastVal;
            for (int i = 0; i < n_of_runs; ++i) {
                IloNum xval = cpl.getValue(x[i]);
                current_CPLEX_individual[i] = (xval > 0.8) ? 1 : 0;
            }
        }
    }
    catch (IloException& e) {
        cerr << " ERROR: " << e << endl;
    }
    env.end();
    return current_CPLEX_individual;
}





/*
void run_cplex(clock_t& start, vector<double>& results, vector<double>& times, vector<double>& gaps, int& na) {

    int n_of_runs = int(runs.size());

    IloEnv env;
    env.setOut(env.getNullStream());
    try{
        IloModel model(env);

        IloIntVarArray x(env, n_of_runs, 0, 1);
        
        IloExpr obj(env);
        for (int i = 0; i < n_of_runs; ++i) obj += runs[i].length*x[i];
        model.add(IloMaximize(env, obj));

        for (int i = 0; i < n_of_runs - 1; ++i) {
            for (int j = i + 1; j < n_of_runs; ++j) {
                if (runs[i].letter == runs[j].letter) {
                    IloExpr expr(env);
                    expr += (j - i)*x[i];
                    expr += (j - i)*x[j];
                    for (int l = i + 1; l < j; ++l) {
                        if (runs[l].letter != runs[i].letter) expr += x[l];
                    }
                    model.add(expr <= 2*(j - i));
                }
            }
        }
        
        IloCplex cpl(model);

        cpl.setParam(IloCplex::TiLim, cplex_time_limit);
        cpl.setParam(IloCplex::EpGap, 0.0);
        cpl.setParam(IloCplex::EpAGap, 0.0);
        cpl.setParam(IloCplex::Threads, 1);
        cpl.use(loggingCallback(env, start, times[na], results[na], gaps[na]));
        cpl.setWarning(env.getNullStream());

        cpl.solve();
    
        if (cpl.getStatus() == IloAlgorithm::Optimal || cpl.getStatus() == IloAlgorithm::Feasible) {
            clock_t current = clock();
            double newTime = double(current - start) / CLOCKS_PER_SEC;
            double lastVal = double(cpl.getObjValue());
            double lastGap = 100.0 * cpl.getMIPRelativeGap();
            if (lastGap < 0.0) lastGap *= -1.0;
            if (lastVal > results[na]) {
                cout << "value " << lastVal << "\ttime " << newTime << "\tgap " << lastGap << endl;
                results[na] = lastVal;
                times[na] = newTime;
                gaps[na] = lastGap;
            }
            if (cpl.getStatus() == IloAlgorithm::Optimal) cout << "optimality proven" << endl;
            cout << "Chosen runs:";
            for (int i = 0; i < n_of_runs; ++i) {
                IloNum xval = cpl.getValue(x[i]);
                if (xval > 0.8) {
                    cout << " " << i;
                }
            }
            cout << endl;
            cout << "Solution string: ";
            for (int i = 0; i < n_of_runs; ++i) {
                IloNum xval = cpl.getValue(x[i]);
                if (xval > 0.8) {
                    for (int j = 0; j < runs[i].length; ++j) cout << alphabet[runs[i].letter];
                }
            }
            cout << endl;
        }
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
    }
    env.end();
}

*/
void show_population_score(vector<Individual>& population) {

    for (int i = 0; i < int(population.size()); ++i) {
        cout << "Individual " << i << " score: " << population[i].score << endl;
    }
}

void show_population(vector<Individual>& population) {

    for (int i = 0; i < int(population.size()); ++i) {
        cout << "Individual " << i << ": ";
        for (int j = 0; j < n_of_runs; ++j) {
            cout << population[i].vec[j] << " ";
        }
        cout << endl;
    }
} 

void show_sub_instance() {
    cout << "Sub-instance: ";
    for (int i = 0; i < int(runs.size()); ++i) {
        cout << sub_instance.runs[i]<<"["<<sub_instance.runs_age[i]<<"] ";
    }
    cout << endl;
}

void local_search(Individual& ind, default_random_engine& generator, int max_iterations) {
    int iterations = 0;
    int n = ind.vec.size();
    set<int> runs_in_solution = ind.selected_runs;

    vector<int> runs_not_in_solution;
    for (int i = 0; i < n_of_runs; ++i) {
        if (runs_in_solution.find(i) == runs_in_solution.end()) {
            runs_not_in_solution.push_back(i);
        }
    }

    // Distribuciones para seleccionar elementos aleatorios
    uniform_int_distribution<int> distribution_in_solution(0, runs_in_solution.size() - 1);
    uniform_int_distribution<int> distribution_not_in_solution(0, runs_not_in_solution.size() - 1);

    while (iterations < max_iterations) {
        // Seleccionar un run aleatorio que está en la solución
        int pos_in_solution = distribution_in_solution(generator);
        auto it = runs_in_solution.begin();
        advance(it, pos_in_solution);
        int run_in_solution = *it;

        int pos_not_in_solution = distribution_not_in_solution(generator);
        int run_not_in_solution = runs_not_in_solution[pos_not_in_solution];
        double temp = ind.vec[run_in_solution];
        ind.vec[run_in_solution] = ind.vec[run_not_in_solution];
        ind.vec[run_not_in_solution] = temp;

        // Evaluar la nueva solución
        int old_score = ind.score;
        evaluate(ind);

        // Verificar si la nueva solución es mejor
        if (ind.score > old_score) {
            iterations = 0; // Reiniciar el contador de iteraciones si hay mejora
        } else {
            // Revertir el intercambio si la solución no es mejor
            swap(ind.vec[run_in_solution], ind.vec[run_not_in_solution]);
            evaluate(ind);
        }

        iterations++;
    }
}



/**********
Main function
**********/

int main(int argc, char **argv) {
    read_parameters(argc, argv);

    // Inicialización del generador de números aleatorios
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> standard_distribution(0.0, 1.0);


    // Configuración de formato de salida
    std::cout << std::setprecision(2) << std::fixed;

    // Bucle principal sobre los archivos de entrada
    for (int na = 0; na < inputFiles.size(); ++na) {
        // Lectura del archivo de entrada y preparación de datos
        ifstream indata(inputFiles[na].c_str());
        if (!indata) {
            cout << "Error: file could not be opened" << endl;
            continue;
        }

        // Procesamiento de la secuencia de entrada
        indata >> orig_input_sequence;
        sequence_length = orig_input_sequence.size();
        
        // Preparación del alfabeto y mapeo inverso
        set<char> chars(orig_input_sequence.begin(), orig_input_sequence.end());
        alphabet_size = chars.size();
        alphabet = vector<char>(chars.begin(), chars.end());
        map<char, int> rev_mapping;
        for (int i = 0; i < alphabet_size; ++i) rev_mapping[alphabet[i]] = i;
        
        // Conversión de la secuencia a índices numéricos
        input_sequence.resize(sequence_length);
        for (int j = 0; j < sequence_length; ++j) {
            input_sequence[j] = rev_mapping[orig_input_sequence[j]];
        }
        indata.close();

        // Identificación de runs en la secuencia
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
        n_of_runs = runs.size();

        // Inicialización del tiempo de cómputo
        start_time = clock();

        // 1. Inicialización de CMSA (Algoritmo 1, líneas 2-3)
        Individual Sbsf; // Mejor solución encontrada
        Sbsf.score = std::numeric_limits<int>::min();
        
        sub_instance.runs = vector<int>(n_of_runs, 0);      // C0 en el pseudocódigo
        sub_instance.runs_age = vector<int>(n_of_runs, 0);  // age[c] en el pseudocódigo
        vector<double> results(inputFiles.size(), std::numeric_limits<int>::min());
        vector<double> times(inputFiles.size(), 0.0);
        vector<double> gaps(inputFiles.size(), 0.0);
        current_CPLEX_individual.resize(n_of_runs, 0);
        computation_time_limit= sequence_length/5;
        


        // Bucle principal de CMSA (Algoritmo 1, línea 4)
        while (double(clock() - start_time) / CLOCKS_PER_SEC < computation_time_limit) {
            // Fase Construct (Algoritmo 1, líneas 5-10)
            for (int i = 0; i < params.n_solutions; ++i) {
                // Generación de solución probabilística (Algoritmo 1, línea 7)
                Individual S;
                generate_greedy_random(S, generator, standard_distribution, 
                                     params.neighborhood_size, params.neighborhood_weight, 
                                     params.determinism_rate);
                local_search(S,generator,params.max_iter_local_search);

                if(Sbsf.score<S.score){
                    Sbsf=S;
                    best_LCMSA_quality = Sbsf.score;
                    double elapsed = double(clock() - start_time) / CLOCKS_PER_SEC;
                    cout << "score " << Sbsf.score << "\ttime " << elapsed << endl;
                }
                // Merge: Actualización de C0 y edades (Algoritmo 1, líneas 8-10)
                for (int run_idx : S.selected_runs) {
                    if (sub_instance.runs[run_idx] == 0) {
                        sub_instance.runs[run_idx] = 1;
                        sub_instance.runs_age[run_idx] = 0;
                    }
                }
            }
            //cout<<"antes de entrar a cplex:"<<endl;
            //show_sub_instance();
            // Fase Solve (Algoritmo 1, línea 12)
            vector<int> S_opt = run_cplex(start_time, results, times, gaps, na);
            //cout<<"saliendo de cplex:"<<endl;
            
            // Evaluación de la solución óptima de la subinstancia
            int S_opt_score = 0;
            for (int i = 0; i < n_of_runs; ++i) {
                if (S_opt[i] == 1) S_opt_score += runs[i].length;
            }
            
            // Actualización de la mejor solución (Algoritmo 1, línea 13)
            if (S_opt_score > Sbsf.score) {
                Sbsf.score = S_opt_score;
                Sbsf.selected_runs.clear();
                for (int i = 0; i < n_of_runs; ++i) {
                    if (S_opt[i] == 1) Sbsf.selected_runs.insert(i);
                }
                best_LCMSA_quality = Sbsf.score;
                
 
                print_solution_in_letters(Sbsf);
            }


            // Fase Adapt 
            for (int i = 0; i < n_of_runs; ++i) {
                if (S_opt[i] == 1) {
                    sub_instance.runs_age[i] = 0; // Reinicia edad si se usó
                } else if (sub_instance.runs[i] == 1) {
                    sub_instance.runs_age[i]++;   // Incrementa edad si no se usó
                    
                    // Elimina componentes viejos (Algoritmo 1, parte implícita de Adapt)
                    if (sub_instance.runs_age[i] > params.age_limit) {
                        sub_instance.runs[i] = 0;
                        sub_instance.runs_age[i] = 0;
                    }
                }
            }
        }
        /*

        // Salida final (Algoritmo 1, línea 16)
        if (output_type == 0) {
            cout << "Final best solution for " << inputFiles[na] << ": " << Sbsf.score << endl;
        } else if (output_type == 2) {
            cout << -1 * Sbsf.score << endl;
        }*/
    }

    return 0;
}
