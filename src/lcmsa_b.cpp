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
clock_t start;
int output_type = 0; //0: normal, 1:STN compatible, 2: Tuning 

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
string model;
vector<double> LLM;




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
        else if (strcmp(argv[iarg],"-p")==0) population_size = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pe")==0) elite_proportion = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pm")==0) mutant_proportion = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-rhoe")==0) elite_inheritance_probability = atof(argv[++iarg]);
        //else if (strcmp(argv[iarg],"-th")==0) threshold = atof(argv[++iarg]);
        //else if (strcmp(argv[iarg],"-trials")==0) trials = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-tuning")==0) tuning = true;
        else if (strcmp(argv[iarg],"-cpl_t")==0) cplex_time_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-t")==0) computation_time_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-iter")==0) max_ga_iterations = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-na")==0) individuals_selected_from_LC_for_merge = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-m_age")==0) max_age = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-n_c")==0) proportion_CPLEX_solutions_copies_injected = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-output_type")==0) output_type = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-model")==0) model = argv[++iarg];
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
		pos_val_pairs[i].second = (ind.vec)[i] * LLM[i];
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


void show_sub_instance() {
    cout << "Sub-instance: ";
    for (int i = 0; i < int(runs.size()); ++i) {
        cout << sub_instance.runs[i]<<"["<<sub_instance.runs_age[i]<<"] ";
    }
    cout << endl;
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
            if(output_type==0) cout << "-*CP- value " << nv << "\ttime " << newTime << "\tgap " << newGap << endl;
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
    try{
        IloModel model(env);

        IloIntVarArray x(env, n_of_runs, 0, 1);
        
        IloExpr obj(env);
        for (int i = 0; i < n_of_runs; ++i) obj += runs[i].length*x[i];
        model.add(IloMaximize(env, obj));
        /*
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
                    expr.end();
                }
            }
        }
        */

        //adding subinstance constraints
        for(int i=0; i<n_of_runs; i++) {
            if(sub_instance.runs[i] == 0) {
                IloExpr expr(env);
                expr += x[i];
                model.add(expr == 0);
                expr.end();
            }
        }



        IloCplex cpl(model);

        cpl.setParam(IloCplex::TiLim, cplex_time_limit);
        cpl.setParam(IloCplex::EpGap, 0.0);
        cpl.setParam(IloCplex::EpAGap, 0.0);
        cpl.setParam(IloCplex::Threads, 1);
        cpl.setParam(IloCplex::PreInd, 0); // Recomendado para lazy constraints

        
        //Experimental
        //********************************************************************************* 
        //cpl.setParam(IloCplex::MemoryEmphasis, 1); // Prioriza el uso eficiente de memoria
        //establish a memory limit for the CPLEX solver
        //cpl.setParam(IloCplex::WorkMem, 4096);
        //*********************************************************************************
        //if (output_type!=1) cpl.use(loggingCallback(env, start, times[na], best_LCMSA_quality, gaps[na]));
        cpl.use(loggingCallback(env, start, x, best_LCMSA_quality, gaps[na]));
        cpl.use(new (env) LazyConstraintCallback(env, x, runs)); // ⬅ Lazy constraints
        
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
                if (output_type==0) cout << "-CP- value " << lastVal << "\ttime " << newTime << "\tgap " << lastGap << endl;
                if(output_type==1) {
                    cout<<"CP,"<<stn_run_number<<",";
                    for (int i = 0; i < n_of_runs; ++i) {
                        IloNum xval = cpl.getValue(x[i]);
                        if (xval > 0.8) cout<<"1"; else cout<<"0";
                        }
                    cout<<","<<-1 * best_LCMSA_quality<<endl;
                    
                }

            }
            current_CPLEX_quality = lastVal;
            //cout << "CPLEX value: " << lastVal << endl;
            //if (cpl.getStatus() == IloAlgorithm::Optimal) cout << "optimality proven" << endl;
            //cout << "Chosen runs:";
            for (int i = 0; i < n_of_runs; ++i) {
                IloNum xval = cpl.getValue(x[i]);
                if (xval > 0.8) {
                    //cout << " " << i;
                    current_CPLEX_individual[i] = 1;
                    //recording the best individual
                    

                } else {
                    current_CPLEX_individual[i] = 0;
                }
            }
            //cout << endl;
            //cout << "Solution string: ";
            //for (int i = 0; i < n_of_runs; ++i) {
                //IloNum xval = cpl.getValue(x[i]);
                //if (xval > 0.8) {
                    //for (int j = 0; j < runs[i].length; ++j) cout << alphabet[runs[i].letter];
                //}
            //}
            //cout << endl;
        }
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
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





/**********
Main function
**********/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);

    // initializes the random number generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> standard_distribution(0.0,1.0);


    
    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // number of input files
    int n_files = int(inputFiles.size());

    // vectors for storing the result and the computation time obtained by the applications of CPLEX
    vector<double> results(n_files, std::numeric_limits<int>::min());
    vector<double> times(n_files, 0.0);
    vector<double> gaps(n_files, 0.0);

    // main loop over all input files (problem instances)
    for (int na = 0; na < n_files; ++na) {

        // opening the corresponding input file and reading the problem data
        ifstream indata;
        indata.open(inputFiles[na].c_str());
        if(!indata) { // file couldn't be opened
            cout << "Error: file could not be opened" << endl;
        }

        map<char,int> rev_mapping;

        indata >> orig_input_sequence;
        sequence_length = int(orig_input_sequence.size());
        //cout << "Input string: " << orig_input_sequence << " length: " << sequence_length << endl;
        set<char> chars;
        for (int i = 0; i < int(orig_input_sequence.size()); ++i) chars.insert(orig_input_sequence[i]);
        //cout << "Sigma = {" << *(chars.begin());
        bool startb = true;
        for (set<char>::iterator sit = chars.begin(); sit != chars.end(); ++sit) {
            if (startb) startb = false;
            //else cout << ", " << *sit;
        }
        //cout << "}" << endl;
        alphabet_size = int(chars.size());
        //cout << "|Sigma| = " << alphabet_size << endl;
        alphabet = vector<char>(alphabet_size);
        int j = 0;
        for (set<char>::iterator sit = chars.begin(); sit != chars.end(); ++sit) {
            alphabet[j] = *sit;
            rev_mapping[alphabet[j]] = j;
            ++j;
        }
        for (int j = 0; j < sequence_length; ++j) input_sequence.push_back(rev_mapping[orig_input_sequence[j]]);
        //cout << "Converted string: ";
        //for (int j = 0; j < sequence_length; ++j) cout << input_sequence[j];
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
        //cout << "All runs: " << endl;
        //for (int i = 0; i < int(runs.size()); ++i) {
        //    cout << "run " << i << " (start: " << runs[i].start << ", letter: " << alphabet[runs[i].letter] << ", length: " << runs[i].length << ")" << endl;
        //}
        n_of_runs = int(runs.size());
        // the computation time starts now
        start = clock();

        //cout << "start file " << inputFiles[na] << endl;
                //lectura probabilidades LLM
        size_t lastSlash = inputFiles[0].find_last_of("/\\");
        string filename = (lastSlash == std::string::npos) ? inputFiles[0] : inputFiles[0].substr(lastSlash + 1);
        string probpath = "LLM_probabilities/" + model + "/" + filename;
        cout<<probpath<<endl;
        ifstream infile(probpath);
        if (!infile.is_open()) {
            probpath = "../../LLM_probabilities/" + model + "/" + filename;
            cout<<probpath<<endl;
            infile.open(probpath); // reutilizamos el mismo `infile`
            if (!infile.is_open()) {
                cerr << "Error al abrir el archivo de probabilidades: " << probpath << endl;
                return 1;
            }
        }

        double value;
        while (infile >> value) {
            LLM.push_back(value);
        }

        infile.close();

        cout << "Se leyeron " << LLM.size() << " valores del archivo." <<endl;

        infile.close();

        //********************************************************************************
        //L_CMSA algorithm BEGIN
        //********************************************************************************
        
        //tuning preparation time: execution time is lenght/10
        
            
        computation_time_limit = int(sequence_length/5);

        // time limit for tuning

        number_CPLEX_solution_copies_injected=int(proportion_CPLEX_solutions_copies_injected*population_size*0.5);
            if(number_CPLEX_solution_copies_injected < 1) number_CPLEX_solution_copies_injected = 1;

        //cout<<"number CPLEX solution copies injected: "<<number_CPLEX_solution_copies_injected<<endl;

        //Initialization Phase

        current_CPLEX_individual = vector<int>(int(runs.size()), 0);
        //initialization of the sub_instance
        for(int i=0; i<runs.size(); i++) {
            sub_instance.runs.push_back(0);
            sub_instance.runs_age.push_back(0);

            
        }

        best_LCMSA_quality = -1;

        int n_elites = int(double(population_size)*elite_proportion);
        if (n_elites < 1) n_elites = 1;

        int n_mutants = int(double(population_size)*mutant_proportion);
        if (n_mutants < 1) n_mutants = 1;

        int n_offspring = population_size - n_elites - n_mutants;
        if (n_offspring < 1) {
            if (not tuning) cout << "OHOHOH: wrong parameter settings" << endl;
            exit(0);
        }
        
         double ctime;
        vector<Individual> population(population_size);
        //vector<Individual> new_population(population_size);
        for (int pi = 0; pi < population_size; ++pi) {
            generate_random_solution(population[pi], generator, standard_distribution);
            if (population[pi].score > best_LCMSA_quality) {
                //temporalmente desactivado
				//for (set<int>::iterator sit = (population[pi].selected_runs).begin(); sit != (population[pi].selected_runs).end(); ++sit) {
				//	for (int j = 0; j < runs[*sit].length; ++j) cout << alphabet[runs[*sit].letter];
				//}
				//cout << endl;
                best_LCMSA_quality = population[pi].score;
                clock_t current2 = clock();
                ctime = double(current2 - start) / CLOCKS_PER_SEC;
                best_BRKGA_fitness = best_LCMSA_quality;
                //times[trial] = ctime;
                if (output_type==0) cout << "-BK- value " << best_BRKGA_fitness << "\ttime " << ctime << endl;
                else if(output_type==1){
                    cout<<"BK,"<<stn_run_number<<",";
                    for(int i=0; i<runs.size(); i++){
                        if(population[pi].selected_runs.find(i) != population[pi].selected_runs.end())
                            cout<<"1"; else cout<<"0";
                        
                        }
                    cout<<","<<-1 * best_BRKGA_fitness<<endl;
                    }
            }
        }
        clock_t current3 = clock();

        ctime = double(current3 - start) / CLOCKS_PER_SEC;
        sort(population.begin(), population.end(), individual_compare);
        //MAIN LOOP of the L_CMSA algorithm
        while (ctime < computation_time_limit) {
            for(int it=0; it<max_ga_iterations; it++) {   
                //LC Phase
                //sort(population.begin(), population.end(), individual_compare);
                //show_population_score(population);
                //exit(0);
                vector<Individual> new_population(population_size);
                //new_population.clear();
                for (int ic = 0; ic < n_elites; ++ic) {
                    new_population[ic].vec = population[ic].vec;
                    new_population[ic].score = population[ic].score;
				    new_population[ic].selected_runs = population[ic].selected_runs;
                }
                for (int ic = 0; ic < n_mutants; ++ic) {
                    generate_random_solution(new_population[n_elites + ic], generator, standard_distribution);

                    if (new_population[n_elites + ic].score > best_LCMSA_quality) {
					    //temporalmente desactivado
                        //for (set<int>::iterator sit = (new_population[n_elites + ic].selected_runs).begin(); sit != (new_population[n_elites + ic].selected_runs).end(); ++sit) {
						//    for (int j = 0; j < runs[*sit].length; ++j) cout << alphabet[runs[*sit].letter];
					//}
                    //cout << endl;
                    best_LCMSA_quality = new_population[n_elites + ic].score;
                    clock_t current4 = clock();
					ctime = double(current4 - start) / CLOCKS_PER_SEC;
                    best_BRKGA_fitness = best_LCMSA_quality;
                    //times[trial] = ctime;
                    if (output_type==0) cout << "-BK- value " << best_BRKGA_fitness << "\ttime " << ctime << endl;
                    else if(output_type==1){
                    cout<<"BK,"<<stn_run_number<<",";
                    for(int i=0; i<runs.size(); i++){
                        if(new_population[n_elites + ic].selected_runs.find(i) != new_population[n_elites + ic].selected_runs.end())
                            cout<<"1"; else cout<<"0";
                        
                        }
                    cout<<","<<-1 * best_BRKGA_fitness<<endl;
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

                    if (new_population[n_elites + n_mutants + ic].score > best_LCMSA_quality) {
                        //temporalmente desactivado
					    //for (set<int>::iterator sit = (new_population[n_elites + n_mutants + ic].selected_runs).begin(); sit != (new_population[n_elites + n_mutants + ic].selected_runs).end(); ++sit) {
						//    for (int j = 0; j < runs[*sit].length; ++j) cout << alphabet[runs[*sit].letter];
					    //}
					//cout << endl;
                        best_LCMSA_quality = new_population[n_elites + n_mutants + ic].score;
					    clock_t current4 = clock();
					    ctime = double(current4 - start) / CLOCKS_PER_SEC;
                        best_BRKGA_fitness = best_LCMSA_quality;
                        //times[trial] = ctime;
                        if (output_type==0) cout << "-BK- value " << best_BRKGA_fitness << "\ttime " << ctime << endl;
                        else if(output_type==1){
                    cout<<"BK,"<<stn_run_number<<",";
                    for(int i=0; i<runs.size(); i++){
                        if(new_population[n_elites + n_mutants + ic].selected_runs.find(i) != new_population[n_elites + n_mutants + ic].selected_runs.end())
                            cout<<"1"; else cout<<"0";
                        
                        }
                    cout<<","<<-1 * best_BRKGA_fitness<<endl;
                    }
                    }
                }
            population.clear();
            population = new_population;
            //new_population.clear();
			//clock_t current5 = clock();
			//ctime = double(current5 - start) / CLOCKS_PER_SEC;
            sort(population.begin(), population.end(), individual_compare);

            }
        
            //Merge Phase
            //We extract the best individual of the BRKGA (first elite) and put the components in the subinstance
            for (set<int>::iterator sit = (population[0].selected_runs).begin(); sit != (population[0].selected_runs).end(); ++sit) {
                    if(sub_instance.runs[*sit] == 0) {
                        sub_instance.runs_age[*sit] = 0;
                    }
                        sub_instance.runs[*sit] = 1;
            }
            
                    
            //Now we select a number of individuals_selected_from_LC_for_merge - 1 random individuals from elite individuals of the BRKGA and put the components in the subinstance
            for(int i=1; i<individuals_selected_from_LC_for_merge; i++) {
                int random_individual = produce_random_integer(n_elites, generator, standard_distribution);
                for (set<int>::iterator sit = (population[random_individual].selected_runs).begin(); sit != (population[random_individual].selected_runs).end(); ++sit) {
                    if(sub_instance.runs[*sit] == 0) {
                        sub_instance.runs_age[*sit] = 0;
                    }
                        sub_instance.runs[*sit] = 1;
                }
            }
                    
            //Solve phase
            //show_sub_instance();
            current_CPLEX_individual=run_cplex(start, results, times, gaps, na);



            //ADAPT Phase
            for(int i=0; i<int(runs.size()); i++) {
                if(current_CPLEX_individual[i] == 0) {
                    sub_instance.runs_age[i]++;
                } else {
                    sub_instance.runs_age[i] = 0;
                }
            }
            
            // we remove the runs that have an age greater than max_age
            for(int i=0; i<int(runs.size()); i++) {
                if(sub_instance.runs_age[i] > max_age) {
                    sub_instance.runs[i] = 0;
                    sub_instance.runs_age[i] = 0;
                }
            }
            //ONLY FOR TESTING (verify where return, and how many copies, may be resort ins needed)
            //Reintegration of the solution from the CPLEX into the BRKGA
            /*
            for(int i=0; i<int(runs.size()); i++) {
                if(current_CPLEX_individual[i] == 1) {
                    population[population_size-1].vec[i] = 0.75;
                } else {
                    population[population_size-1].vec[i] = 0.25;
                }
            }
            int eval=evaluate(population[population_size-1]);
            */
            //reintegration of number_CPLEX_solution_copies_injected copies of the solution from the CPLEX into the BRKGA
            //in random positions
            

            for(int i=0; i<number_CPLEX_solution_copies_injected; i++) {
                int random_position = produce_random_integer(population_size, generator, standard_distribution);
                for(int j=0; j<int(runs.size()); j++) {
                    if(current_CPLEX_individual[j] == 1) {
                        population[random_position].vec[j] = 0.75;
                    } else {
                        population[random_position].vec[j] = 0.25;
                    }
                }
                evaluate(population[random_position]);
            }
            
            
            
            
            //cout<<"eval de solución ingresada: "<<eval<<endl; 
            //show_population_score(population);
            //mostrar individuos para verificar inserciòn
            //show_population(population);
            sort(population.begin(), population.end(), individual_compare);
            //show_population(population);
            //exit(0);


            clock_t current5 = clock();
			ctime = double(current5 - start) / CLOCKS_PER_SEC;
            //cout<<"time: "<<ctime<<endl;

        //stn_run_number++;
        }
        


        //run_cplex(start, results, times, gaps, na);

         if (output_type != 2) cout << "end file " << inputFiles[na] << endl;
         if(output_type==2) cout<<-1*best_LCMSA_quality<<endl;

    }

    // calculating the average of the results, computation times, and optimality gaps and write them to the screen
    /*
    double r_mean = 0.0;
    double t_mean = 0.0;
    double g_mean = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        r_mean = r_mean + results[i];
        t_mean = t_mean + times[i];
        g_mean = g_mean + gaps[i];
    }
    r_mean = r_mean/double(results.size());
    t_mean = t_mean/double(times.size());
    g_mean = g_mean/double(gaps.size());
    cout << r_mean << "\t" << t_mean << "\t" << g_mean << endl;
    */
  

}

