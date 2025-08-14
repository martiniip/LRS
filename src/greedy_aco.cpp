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


void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFile = argv[++iarg];
        iarg++;
    }
}

void generate_greedy_solution(Solution& iSol) {
    // Crear un vector de pares (índice, longitud de la run)
    vector< pair<int, double> > run_length_pairs(n_of_runs);

    // Llenar el vector con las longitudes de las runs
    for (int i = 0; i < n_of_runs; ++i) {
        run_length_pairs[i].first = i; // Índice de la run
        run_length_pairs[i].second = double(runs[i].length); // Longitud de la run
    }

    // Ordenar las runs de mayor a menor según su longitud
    sort(run_length_pairs.begin(), run_length_pairs.end(), [](const pair<int, double>& a, const pair<int, double>& b) {
        return a.second > b.second; // Orden descendente
    });


    // Inicializar los límites de las letras
    vector<int> letter_lbs(alphabet_size, -1);
    vector<int> letter_ubs(alphabet_size, -1);

    // Seleccionar las runs de manera greedy
    for (int i = 0; i < n_of_runs; ++i) {
        int current_run = run_length_pairs[i].first;

        // Verificar si la run puede ser añadida sin conflictos
        bool add = true;
        for (int j = 0; add and (j < alphabet_size); ++j) {
            if (j != runs[current_run].letter) {
                add = false;
                if (!(current_run > letter_lbs[j] && current_run < letter_ubs[j])) {
                    if (!(current_run < letter_lbs[j] && letter_lbs[runs[current_run].letter] > letter_ubs[j])) {
                        if (!(current_run > letter_ubs[j] && letter_ubs[runs[current_run].letter] < letter_lbs[j])) {
                            add = true;
                        }
                    }
                }
            }
        }

        // Si no hay conflictos, añadir la run a la solución
        if (add) {
            iSol.vec[current_run] = 1; // Marcar la run como seleccionada
            iSol.score += runs[current_run].length; // Actualizar el puntaje

            // Actualizar los límites de la letra correspondiente
            if (letter_lbs[runs[current_run].letter] == -1) {
                letter_lbs[runs[current_run].letter] = current_run;
                letter_ubs[runs[current_run].letter] = current_run;
            } else {
                if (current_run < letter_lbs[runs[current_run].letter]) {
                    letter_lbs[runs[current_run].letter] = current_run;
                } else if (current_run > letter_ubs[runs[current_run].letter]) {
                    letter_ubs[runs[current_run].letter] = current_run;
                }
            }
        }
    }
}
void print_solution_in_letters(const Solution& iSol) {
    // Crear una copia de la secuencia original
    string solution_sequence(orig_input_sequence);

    // Primero, imprimir las runs seleccionadas en el orden deseado
    cout << "Solution sequence (selected runs in uppercase):" << endl;
    for (int i = 0; i < n_of_runs; ++i) {
        if (iSol.vec[i] == 1) { // Si la run está seleccionada
            for (int j = runs[i].start; j < runs[i].start + runs[i].length; ++j) {
                cout << (char)toupper(solution_sequence[j]); // Convertir a mayúsculas e imprimir
            }
        }
    }

    // Luego, imprimir el resto de la secuencia en minúsculas
    for (int i = 0; i < solution_sequence.length(); ++i) {
        bool is_part_of_selected_run = false;
        for (int j = 0; j < n_of_runs; ++j) {
            if (iSol.vec[j] == 1 && i >= runs[j].start && i < runs[j].start + runs[j].length) {
                is_part_of_selected_run = true;
                break;
            }
        }
        if (!is_part_of_selected_run) {
            cout << (char)tolower(solution_sequence[i]); // Imprimir en minúsculas
        }
    }

    cout << endl;
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

   

    // start greedy
    
    cout << "start greedy"<< endl;

    clock_t start = clock();

    Solution bestSol;
    bestSol.score = 0;
    bestSol.vec.resize(n_of_runs, 0);    

    generate_greedy_solution(bestSol);

    clock_t current = clock();
    double ctime = double(current - start) / CLOCKS_PER_SEC;

    // Imprimir la solución en forma de letras
    print_solution_in_letters(bestSol);

    cout << "score " << bestSol.score << "\ttime " << ctime << endl;

}