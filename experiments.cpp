#include "experiments.h"

/* ------------------------------------------------------------------ */
/* Real-time NEAT Validation                                          */
/* ------------------------------------------------------------------ */

//Perform evolution on double pole balacing using rtNEAT methods calls
//Always uses Markov case (i.e. velocities provided)
//This test is meant to validate the rtNEAT methods and show how they can be used instead
// of the usual generational NEAT
Population *pole2TestRealTime() {
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
    pole2RealTimeLoop(pop, thecart);

    return pop;

}

int pole2RealTimeLoop(Population *pop, CartPole *cart) {
    vector<Organism *>::iterator curorg;
    vector<Species *>::iterator curspec; //used in printing out debug info

    int pause;
    bool win = false;

    //Real-time evolution variables
    int offspring_count;
    Organism *new_org;

    cart->N_MARKOV_LONG = false;
    cart->GENERALIZATION_TEST = false;

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

        if (pole2Evaluate((*curorg), cart)) win = true;

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

            int num_species = (int) pop->species.size();
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
        if (pole2Evaluate(new_org, cart)) win = true;

        if (win) {
            cout << "WINNER" << endl;
            pop->print_to_file_by_species((char *) "rt_winpop");
            break;
        }

        //Now we reestimate the baby's species' fitness
        new_org->species->estimate_average();

        //Remove the worst organism
        pop->remove_worst();

    }
    return 0;
}

bool pole2Evaluate(Organism *org, CartPole *cart) {
    Network *net;

    int pause;

    net = org->net;

    //Try to balance a pole now
    org->fitness = cart->evalNet(net);

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

    if ((!(cart->GENERALIZATION_TEST)) && (!(cart->N_MARKOV_LONG)))
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
    if (cart->MARKOV) {
        if (org->fitness >= (cart->MAX_FITNESS - 1)) {
            org->winner = true;
            return true;
        } else {
            org->winner = false;
            return false;
        }
    }
        //if doing the long test non-markov
    else if (cart->N_MARKOV_LONG) {
        if (org->fitness >= 99999) {
            org->winner = true;
            return true;
        } else {
            org->winner = false;
            return false;
        }
    } else if (cart->GENERALIZATION_TEST) {
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
