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
  // the following "include" is necessary for the correct working/compilation of CPLEX.
#include <ilcplex/ilocplex.h>

using namespace std;

ILOSTLBEGIN

struct Run {
    int start;
    int letter;
    int length;
};


// Data structures for the problem data
string orig_input_sequence;
vector<int> input_sequence;
vector<char> alphabet;
int alphabet_size;
int sequence_length;
vector<Run> runs;

// vector for keeping all the names of the input files
vector<string> inputFiles;

// time limit for CPLEX (can be supplied to the algorithm via the -t comand line parameter)
double time_limit = 3200.0;


inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFiles.push_back(argv[++iarg]);
        else if (strcmp(argv[iarg],"-t")==0) time_limit = atof(argv[++iarg]);
        iarg++;
    }
}

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

        cpl.setParam(IloCplex::TiLim, time_limit);
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



/**********
Main function
**********/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);
    
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

        // the computation time starts now
        clock_t start = clock();

        cout << "start file " << inputFiles[na] << endl;

        run_cplex(start, results, times, gaps, na);

        cout << "end file " << inputFiles[na] << endl;
    }

    // calculating the average of the results, computation times, and optimality gaps and write them to the screen
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
}

