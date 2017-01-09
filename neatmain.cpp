#include <iostream>
#include "neat.h"
#include "population.h"
#include "experiments.h"

using namespace std;

int main(int argc, char *argv[]) {
    srand((unsigned) time(NULL));
    if (argc != 2) {
        cerr << "A NEAT parameters file (.ne file) is required to run the experiments!" << endl;
        return -1;
    }

    NEAT::load_neat_params(argv[1], true);
    cout << "loaded" << endl;
    cout << "Press any key: ";
    getchar();

    NEAT::Population *p = 0;
    p = pole2TestRealTime();
    if (p)
        delete p;

    return 0;
}

