#include "experiments.h"

CartPole::CartPole(bool velocity) {
    maxFitness = 100000;

    MARKOV = velocity;

    LENGTH_2 = 0.05;
    MASSPOLE_2 = 0.01;

}

double CartPole::evalNet(Network *net) {
    int steps = 0;
    double input[NUM_INPUTS];
    double output;

    int nmarkovmax;

    double nmarkov_fitness;

    double jiggletotal; //total jiggle in last 100
    int count;  //step counter

    //init(randomize);		// restart at some point

    if (nmarkov_long) nmarkovmax = 100000;
    else if (generalization_test) nmarkovmax = 1000;
    else nmarkovmax = 1000;


    init();

    if (MARKOV) {
        while (steps++ < maxFitness) {


            input[0] = state[0] / 4.8;
            input[1] = state[1] / 2;
            input[2] = state[2] / 0.52;
            input[3] = state[3] / 2;
            input[4] = state[4] / 0.52;
            input[5] = state[5] / 2;
            input[6] = .5;

            net->load_sensors(input);

            //Activate the net
            //If it loops, exit returning only fitness of 1 step
            if (!(net->activate())) return 1.0;

            output = (*(net->outputs.begin()))->activation;

            performAction(output, steps);

            if (outsideBounds())    // if failure
                break;            // stop it now
        }
        return (double) steps;
    } else {  //NON MARKOV CASE

        while (steps++ < nmarkovmax) {

            input[0] = state[0] / 4.8;
            input[1] = state[2] / 0.52;
            input[2] = state[4] / 0.52;
            input[3] = .5;

            net->load_sensors(input);

            //Activate the net
            //If it loops, exit returning only fitness of 1 step
            if (!(net->activate())) return 0.0001;

            output = (*(net->outputs.begin()))->activation;

            performAction(output, steps);

            if (outsideBounds())    // if failure
                break;            // stop it now

            if (nmarkov_long && (outsideBounds()))    // if failure
                break;            // stop it now
        }

        //If we are generalizing we just need to balance it a while
        if (generalization_test)
            return (double) balanced_sum;

        //Sum last 100
        if ((steps > 100) && (!nmarkov_long)) {

            jiggletotal = 0;
            cout << "step " << steps - 99 - 2 << " to step " << steps - 2 << endl;
            //Adjust for array bounds and count
            for (count = steps - 99 - 2; count <= steps - 2; count++)
                jiggletotal += jigglestep[count];
        }

        if (!nmarkov_long) {
            if (balanced_sum > 100)
                nmarkov_fitness = ((0.1 * (((double) balanced_sum) / 1000.0)) +
                                   (0.9 * (0.75 / (jiggletotal))));
            else nmarkov_fitness = (0.1 * (((double) balanced_sum) / 1000.0));

#ifndef NO_SCREEN_OUTR
            cout << "Balanced:  " << balanced_sum << " jiggle: " << jiggletotal << " ***" << endl;
#endif

            return nmarkov_fitness;
        } else return (double) steps;

    }

}

void CartPole::init() {
    static int first_time = 1;

    if (!MARKOV) {
        //Clear all fitness records
        cartpos_sum = 0.0;
        cartv_sum = 0.0;
        polepos_sum = 0.0;
        polev_sum = 0.0;
    }

    balanced_sum = 0; //Always count # balanced

    /*if (randomize) {
      state[0] = (lrand48()%4800)/1000.0 - 2.4;
      state[1] = (lrand48()%2000)/1000.0 - 1;
      state[2] = (lrand48()%400)/1000.0 - 0.2;
      state[3] = (lrand48()%400)/1000.0 - 0.2;
      state[4] = (lrand48()%3000)/1000.0 - 1.5;
      state[5] = (lrand48()%3000)/1000.0 - 1.5;
    }
    else {*/


    if (!generalization_test) {
        state[0] = state[1] = state[3] = state[4] = state[5] = 0;
        state[2] = 0.07; // one_degree;
    } else {
        state[4] = state[5] = 0;
    }

    //}
    if (first_time) {
        cout << "Initial Long pole angle = %f\n" << state[2] << endl;;
        cout << "Initial Short pole length = %f\n" << LENGTH_2 << endl;
        first_time = 0;
    }
}

void CartPole::performAction(double output, int stepnum) {

    int i;
    double dydx[6];

    const bool RK4 = true; //Set to Runge-Kutta 4th order integration method
    const double EULER_TAU = TAU / 4;

    /*random start state for long pole*/
    /*state[2]= drand48();   */

    /*--- Apply action to the simulated cart-pole ---*/

    if (RK4) {
        for (i = 0; i < 2; ++i) {
            dydx[0] = state[1];
            dydx[2] = state[3];
            dydx[4] = state[5];
            step(output, state, dydx);
            rk4(output, state, dydx, state);
        }
    }

    //Record this state
    cartpos_sum += fabs(state[0]);
    cartv_sum += fabs(state[1]);
    polepos_sum += fabs(state[2]);
    polev_sum += fabs(state[3]);
    if (stepnum <= 1000)
        jigglestep[stepnum - 1] = fabs(state[0]) + fabs(state[1]) + fabs(state[2]) + fabs(state[3]);

    if (!(outsideBounds()))
        ++balanced_sum;

}

void CartPole::step(double action, double *st, double *derivs) {
    double force, costheta_1, costheta_2, sintheta_1, sintheta_2,
            gsintheta_1, gsintheta_2, temp_1, temp_2, ml_1, ml_2, fi_1, fi_2, mi_1, mi_2;

    force = (action - 0.5) * FORCE_MAG * 2;
    costheta_1 = cos(st[2]);
    sintheta_1 = sin(st[2]);
    gsintheta_1 = GRAVITY * sintheta_1;
    costheta_2 = cos(st[4]);
    sintheta_2 = sin(st[4]);
    gsintheta_2 = GRAVITY * sintheta_2;

    ml_1 = LENGTH_1 * MASSPOLE_1;
    ml_2 = LENGTH_2 * MASSPOLE_2;
    temp_1 = MUP * st[3] / ml_1;
    temp_2 = MUP * st[5] / ml_2;
    fi_1 = (ml_1 * st[3] * st[3] * sintheta_1) +
           (0.75 * MASSPOLE_1 * costheta_1 * (temp_1 + gsintheta_1));
    fi_2 = (ml_2 * st[5] * st[5] * sintheta_2) +
           (0.75 * MASSPOLE_2 * costheta_2 * (temp_2 + gsintheta_2));
    mi_1 = MASSPOLE_1 * (1 - (0.75 * costheta_1 * costheta_1));
    mi_2 = MASSPOLE_2 * (1 - (0.75 * costheta_2 * costheta_2));

    derivs[1] = (force + fi_1 + fi_2)
                / (mi_1 + mi_2 + MASSCART);

    derivs[3] = -0.75 * (derivs[1] * costheta_1 + gsintheta_1 + temp_1)
                / LENGTH_1;
    derivs[5] = -0.75 * (derivs[1] * costheta_2 + gsintheta_2 + temp_2)
                / LENGTH_2;

}

void CartPole::rk4(double f, double y[], double dydx[], double yout[]) {

    int i;

    double hh, h6, dym[6], dyt[6], yt[6];


    hh = TAU * 0.5;
    h6 = TAU / 6.0;
    for (i = 0; i <= 5; i++) yt[i] = y[i] + hh * dydx[i];
    step(f, yt, dyt);
    dyt[0] = yt[1];
    dyt[2] = yt[3];
    dyt[4] = yt[5];
    for (i = 0; i <= 5; i++) yt[i] = y[i] + hh * dyt[i];
    step(f, yt, dym);
    dym[0] = yt[1];
    dym[2] = yt[3];
    dym[4] = yt[5];
    for (i = 0; i <= 5; i++) {
        yt[i] = y[i] + TAU * dym[i];
        dym[i] += dyt[i];
    }
    step(f, yt, dyt);
    dyt[0] = yt[1];
    dyt[2] = yt[3];
    dyt[4] = yt[5];
    for (i = 0; i <= 5; i++)
        yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
}

bool CartPole::outsideBounds() {
    const double failureAngle = thirty_six_degrees;

    return
            state[0] < -2.4 ||
            state[0] > 2.4 ||
            state[2] < -failureAngle ||
            state[2] > failureAngle ||
            state[4] < -failureAngle ||
            state[4] > failureAngle;
}

/* ------------------------------------------------------------------ */
/* Real-time NEAT Validation                                          */
/* ------------------------------------------------------------------ */

//Perform evolution on double pole balacing using rtNEAT methods calls
//Always uses Markov case (i.e. velocities provided)
//This test is meant to validate the rtNEAT methods and show how they can be used instead
// of the usual generational NEAT
Population *pole2_test_realtime() {
    Population *pop;
    Genome *start_genome;
    char curword[20];
    int id;

    CartPole *thecart;

    ifstream iFile("pole2startgenes1", ios::in);

    cout << "START DOUBLE POLE BALANCING REAL-TIME EVOLUTION VALIDATION" << endl;

    cout << "Reading in the start genome" << endl;

    //Read in the start Genome
    iFile >> curword;
    iFile >> id;
    cout << "Reading in Genome id " << id << endl;
    start_genome = new Genome(id, iFile);
    iFile.close();

    cout << "Start Genome: " << start_genome << endl;

    //Spawn the Population from starter gene
    cout << "Spawning Population off Genome" << endl;
    pop = new Population(start_genome, NEAT::pop_size);

    cout << "Verifying Spawned Pop" << endl;
    pop->verify();

    //Create the Cart
    thecart = new CartPole(1);

    //Start the evolution loop using rtNEAT method calls 
    pole2_realtime_loop(pop, thecart);

    return pop;

}

int pole2_realtime_loop(Population *pop, CartPole *thecart) {
    vector<Organism *>::iterator curorg;
    vector<Species *>::iterator curspec; //used in printing out debug info

    int pause;
    bool win = false;

    //Real-time evolution variables
    int offspring_count;
    Organism *new_org;

    thecart->nmarkov_long = false;
    thecart->generalization_test = false;

    //We try to keep the number of species constant at this number
    int num_species_target = NEAT::pop_size / 15;

    //This is where we determine the frequency of compatibility threshold adjustment
    int compat_adjust_frequency = NEAT::pop_size / 10;
    if (compat_adjust_frequency < 1)
        compat_adjust_frequency = 1;

    //Initially, we evaluate the whole population
    //Evaluate each organism on a test
    for (curorg = (pop->organisms).begin(); curorg != (pop->organisms).end(); ++curorg) {

        //shouldn't happen
        if (((*curorg)->gnome) == 0) {
            cout << "ERROR EMPTY GEMOME!" << endl;
            cin >> pause;
        }

        if (pole2_evaluate((*curorg), thecart)) win = true;

    }

    //Get ready for real-time loop

    //Rank all the organisms from best to worst in each species
    pop->rank_within_species();

    //Assign each species an average fitness
    //This average must be kept up-to-date by rtNEAT in order to select species probabailistically for reproduction
    pop->estimate_all_averages();


    //Now create offspring one at a time, testing each offspring,
    // and replacing the worst with the new offspring if its better
    for (offspring_count = 0; offspring_count < 20000; offspring_count++) {


        //Every pop_size reproductions, adjust the compat_thresh to better match the num_species_targer
        //and reassign the population to new species
        if (offspring_count % compat_adjust_frequency == 0) {

            int num_species = pop->species.size();
            double compat_mod = 0.1;  //Modify compat thresh to control speciation

            // This tinkers with the compatibility threshold
            if (num_species < num_species_target) {
                NEAT::compat_threshold -= compat_mod;
            } else if (num_species > num_species_target)
                NEAT::compat_threshold += compat_mod;

            if (NEAT::compat_threshold < 0.3)
                NEAT::compat_threshold = 0.3;

            cout << "compat_thresh = " << NEAT::compat_threshold << endl;

            //Go through entire population, reassigning organisms to new species
            for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
                pop->reassign_species(*curorg);
            }
        }


        //For printing only
        for (curspec = (pop->species).begin(); curspec != (pop->species).end(); curspec++) {
            cout << "Species " << (*curspec)->id << " size" << (*curspec)->organisms.size() << " average= "
                 << (*curspec)->average_est << endl;
        }

        cout << "Pop size: " << pop->organisms.size() << endl;

        //Here we call two rtNEAT calls:
        //1) choose_parent_species() decides which species should produce the next offspring
        //2) reproduce_one(...) creates a single offspring fromt the chosen species
        new_org = (pop->choose_parent_species())->reproduce_one(offspring_count, pop, pop->species);

        //Now we evaluate the new individual
        //Note that in a true real-time simulation, evaluation would be happening to all individuals at all times.
        //That is, this call would not appear here in a true online simulation.
        cout << "Evaluating new baby: " << endl;
        if (pole2_evaluate(new_org, thecart)) win = true;

        if (win) {
            cout << "WINNER" << endl;
            pop->print_to_file_by_species("rt_winpop");
            break;
        }

        //Now we reestimate the baby's species' fitness
        new_org->species->estimate_average();

        //Remove the worst organism
        pop->remove_worst();

    }
    return 0;
}

bool pole2_evaluate(Organism *org, CartPole *thecart) {
    Network *net;

    int pause;

    net = org->net;

    //Try to balance a pole now
    org->fitness = thecart->evalNet(net);

#ifndef NO_SCREEN_OUT
    if (org->pop_champ_child)
        cout << " <<DUPLICATE OF CHAMPION>> ";

    //Output to screen
    cout << "Org " << (org->gnome)->genome_id << " fitness: " << org->fitness;
    cout << " (" << (org->gnome)->genes.size();
    cout << " / " << (org->gnome)->nodes.size() << ")";
    cout << "   ";
    if (org->mut_struct_baby) cout << " [struct]";
    if (org->mate_baby) cout << " [mate]";
    cout << endl;
#endif

    if ((!(thecart->generalization_test)) && (!(thecart->nmarkov_long)))
        if (org->pop_champ_child) {
            cout << org->gnome << endl;
            //DEBUG CHECK
            if (org->high_fit > org->fitness) {
                cout << "ALERT: ORGANISM DAMAGED" << endl;
                print_Genome_tofile(org->gnome, "failure_champ_genome");
                cin >> pause;
            }
        }

    //Decide if its a winner, in Markov Case
    if (thecart->MARKOV) {
        if (org->fitness >= (thecart->maxFitness - 1)) {
            org->winner = true;
            return true;
        } else {
            org->winner = false;
            return false;
        }
    }
        //if doing the long test non-markov
    else if (thecart->nmarkov_long) {
        if (org->fitness >= 99999) {
            //if (org->fitness>=9000) {
            org->winner = true;
            return true;
        } else {
            org->winner = false;
            return false;
        }
    } else if (thecart->generalization_test) {
        if (org->fitness >= 999) {
            org->winner = true;
            return true;
        } else {
            org->winner = false;
            return false;
        }
    } else {
        org->winner = false;
        return false;  //Winners not decided here in non-Markov
    }
}
