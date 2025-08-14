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
#include <unordered_set>

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

// Data structures for the problem data
string orig_input_sequence;
vector<int> input_sequence;
vector<char> alphabet;
vector<Run> sol;
int alphabet_size;
int sequence_length;
vector<Run> runs;
int n_of_runs;

// Other data structures
string inputFile;
int k = 240;
double ponderacion = 0.19;
double determinism_level=0.42;
double computation_time_limit = 5.0;
int max_iterations = 1000;
double time_threshold = 7.0; // Umbral de tiempo sin mejoras (en segundos)
bool newheuristic = false;

void read_parameters(int argc, char **argv) {
    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg], "-i") == 0) inputFile = argv[++iarg];
        else if (strcmp(argv[iarg], "-max_iter") == 0) {
            max_iterations = atoi(argv[++iarg]);
        }else if (strcmp(argv[iarg], "-d") == 0) {
            determinism_level = std::stod(argv[++iarg]);
        }else if (strcmp(argv[iarg], "-time_limit") == 0) {
            computation_time_limit = std::stod(argv[++iarg]);
        }else if (strcmp(argv[iarg], "-time_th") == 0) {
            time_threshold = std::stod(argv[++iarg]);
        }        
        else if (strcmp(argv[iarg],"-newheuristic")==0) newheuristic = true;
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
    uniform_real_distribution<double> dis_random(0.0, 1.0);

    for (int i = 0; i < n_of_runs; ++i) {
        if (dis_random(generator) > determinism_level) {
            // Caso aleatorio: asignar un valor aleatorio
            (ind.vec)[i] = standard_distribution(generator);
        } else {
            // Caso greedy: calcular el score basado en la longitud de la run y sus vecinos
            int current_run = i;

            //HEURISTICA BASE 
            
            double score = runs[current_run].length; // Puntaje base: Longitud propia

            
            //HEURISTICA NUEVA
            if(newheuristic){
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


int main(int argc, char **argv) {
    read_parameters(argc, argv);

    std::cout << std::setprecision(3) << std::fixed;

    // reading an instance
    ifstream indata;
    indata.open(inputFile.c_str());
    if (!indata) {
        cout << "Error: file could not be opened" << endl;
        exit(1);
    }

    map<char, int> rev_mapping;

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

    // start greedy
    cout << "start greedy" << endl;

    clock_t start = clock();

    Individual bestInd;
    bestInd.score = 0;
    bestInd.vec.resize(n_of_runs, 0.0);
    computation_time_limit = sequence_length / 5;
    //computation_time_limit = 20;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> standard_distribution(0.0,20.0);
    clock_t current = clock();
    double ctime = double(current - start) / CLOCKS_PER_SEC;
    // Variables para controlar el tiempo sin mejoras
    double time_without_improvement = 0.0;
    bool use_local_search = true; // Bandera para activar la búsqueda local

    //clock_t last_improvement_time = start; // Tiempo de la última mejora

    while (true) {
        clock_t current = clock();
        double ctime = double(current - start) / CLOCKS_PER_SEC;
        //cout<<ctime<<endl;

        // Si se supera el límite de tiempo, salir del bucle
        if (ctime >= computation_time_limit) {
            break;
        }

        // Generar una solución greedy aleatoria
        Individual randInd;
        randInd.score = 0;
        randInd.vec.resize(n_of_runs, 0.0);
        generate_greedy_random(randInd, generator, standard_distribution, k, ponderacion, determinism_level);

        // Si se activó la búsqueda local, aplicarla a la solución generada
        //if (use_local_search) {
        local_search(randInd, generator, max_iterations);
        //}

        // Verificar si la nueva solución es mejor
        if (randInd.score > bestInd.score) {
            bestInd = randInd;
            print_solution_in_letters(bestInd);
            cout << "score " << bestInd.score << "\ttime " << ctime << endl;

            // Reiniciar el contador de tiempo sin mejoras
            //last_improvement_time = current;
            //time_without_improvement = 0.0;
        } 
/*
        //else {
            // Actualizar el tiempo sin mejoras
            //time_without_improvement = double(current - last_improvement_time) / CLOCKS_PER_SEC;
        //}

        // Si el tiempo sin mejoras supera el umbral, activar la búsqueda local
        if (!use_local_search && time_without_improvement >= time_threshold) {
            use_local_search = true;
            //cout << "Activando búsqueda local..." << endl;
        }
    */
    }
    cout<<bestInd.score<<endl;

    return 0;
}