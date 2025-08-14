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
struct Individual {
    vector<double> vec;
    set<int> selected_runs;
    set<int> onlyselected;
    int score;
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
//4length
int k = 240;
double ponderacion=0.19;
double determinism_level=1;



void read_parameters(int argc, char **argv) {
    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg], "-i") == 0) inputFile = argv[++iarg];
        else if (strcmp(argv[iarg], "-pond") == 0) {
            ponderacion = std::stod(argv[++iarg]);  // Cambié std::stoi por std::stod para convertir a double
        }
        else if (strcmp(argv[iarg], "-d") == 0) {
            determinism_level = std::stod(argv[++iarg]);  // Cambié std::stoi por std::stod para convertir a double
        }
        iarg++;
    }
}


bool individual_compare(const Individual& i1, const Individual& i2) {

    return  i1.score > i2.score;
}

bool pos_val_compare(const pair<int,double>& i1, const pair<int,double>& i2) {

    return  i1.second > i2.second;
}

void evaluate(Individual& ind) {

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
}

void generate_random_solution(Individual& ind, default_random_engine& generator, uniform_real_distribution<double>& standard_distribution) {

    ind.vec = vector<double>(n_of_runs);
    for (int i = 0; i < n_of_runs; ++i) (ind.vec)[i] = standard_distribution(generator);
    evaluate(ind);
}


void generate_greedy_random(Individual& ind, default_random_engine& generator, uniform_real_distribution<double>& standard_distribution, int k, double ponderacion, double determinism_level) {
    // Inicializar el vector de scores
    ind.vec = vector<double>(n_of_runs);

    // Generar números aleatorios para decidir si usar greedy o aleatorio
    uniform_real_distribution<double> dis_random(0.0, 0.0);

    for (int i = 0; i < n_of_runs; ++i) {
        if (dis_random(generator) > determinism_level) {
            // Caso aleatorio: asignar un valor aleatorio
            (ind.vec)[i] = standard_distribution(generator);
        } else {
            // Caso greedy: calcular el score basado en la longitud de la run y sus vecinos
            int current_run = i;
            double score = runs[current_run].length; // Puntaje base: Longitud propia

            // Calcular el puntaje adicional por las runs vecinas con la misma letra
            for (int j = 1; j <= k; ++j) {
                // Verificar vecinos hacia la izquierda (i - j)
                if (current_run - j >= 0 && runs[current_run - j].letter == runs[current_run].letter) {
                    score += ponderacion * runs[current_run - j].length;
                }
                // Verificar vecinos hacia la derecha (i + j)
                if (current_run + j < n_of_runs && runs[current_run + j].letter == runs[current_run].letter) {
                    score += ponderacion * runs[current_run + j].length;
                }
            }

            // Asignar el score calculado
            (ind.vec)[i] = score;
        }
    }

    // Evaluar la solución generada
    evaluate(ind);
}


void print_solution_in_letters(const Individual& ind) {
    for (set<int>::iterator sit = (ind.selected_runs).begin(); sit != (ind.selected_runs).end(); ++sit) {
                    for (int j = 0; j < runs[*sit].length; ++j) cout << alphabet[runs[*sit].letter];
    }
    cout<<endl;
}


int main( int argc, char **argv ) {

    read_parameters(argc,argv);

    std::cout << std::setprecision(3) << std::fixed;

    // reading an instance
    ifstream indata;
    indata.open(inputFile.c_str());
    if(!indata) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
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

   
    //neighborhood_size= int(4*alphabet_size);
    // start greedy
    
    cout << "start greedy"<< endl;

    clock_t start = clock();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> standard_distribution(0.0,20.0);

    Individual bestSol;
    bestSol.score = 0;
    bestSol.vec.resize(n_of_runs, 0);    
    //generate_greedy_solution1(bestSol,neighborhood_size,ponderacion);
    generate_greedy_random(bestSol, generator, standard_distribution, 240, 0.19, 1);
    //generate_greedy_solution(bestSol,neighborhood_size);

    clock_t current = clock();
    double ctime = double(current - start) / CLOCKS_PER_SEC;

    // Imprimir la solución en forma de letras
    print_solution_in_letters(bestSol);

    cout << "score " << bestSol.score << "\ttime " << ctime << endl;

}