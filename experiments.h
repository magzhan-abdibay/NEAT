#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "neat.h"
#include "network.h"
#include "population.h"
#include "organism.h"
#include "genome.h"
#include "species.h"

using namespace std;

using namespace NEAT;

class CartPole;

Population *pole2_test_realtime();

int pole2_realtime_loop(Population *pop, CartPole *thecart);

bool pole2_evaluate(Organism *org, CartPole *thecart);

class CartPole {
public:
    CartPole(bool velocity);

    virtual double evalNet(Network *net);

    double maxFitness;
    bool MARKOV;

    bool nmarkov_long;  //Flag that we are looking at the champ
    bool generalization_test;  //Flag we are testing champ's generalization

    double state[6];

    double jigglestep[1000];

protected:
    virtual void init();

private:

    void performAction(double output, int stepnum);

    void step(double action, double *state, double *derivs);

    void rk4(double f, double y[], double dydx[], double yout[]);

    bool outsideBounds();

    const static int NUM_INPUTS = 7;
    const static double MUP = 0.000002;
    const static double GRAVITY = -9.8;
    const static double MASSCART = 1.0;
    const static double MASSPOLE_1 = 0.1;

    const static double LENGTH_1 = 0.5;          /* actually half the pole's length */

    const static double FORCE_MAG = 10.0;
    const static double TAU = 0.01;          //seconds between state updates

    const static double thirty_six_degrees = 0.628329;

    double LENGTH_2;
    double MASSPOLE_2;

    //Queues used for Gruau's fitness which damps oscillations
    int balanced_sum;
    double cartpos_sum;
    double cartv_sum;
    double polepos_sum;
    double polev_sum;


};

#endif
