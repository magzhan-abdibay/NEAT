#ifndef _POPULATION_H_
#define _POPULATION_H_

#include <cmath>
#include <vector>
#include "Innovation.h"
#include "Genome.h"
#include "Species.h"
#include "organism.h"

namespace NEAT {

    class Species;

    class Organism;

    // ---------------------------------------------
    // POPULATION CLASS:
    //   A Population is a group of Organisms
    //   including their species
    // ---------------------------------------------
    class Population {

    protected:

        // A Population can be spawned off of a single Genome
        // There will be size Genomes added to the Population
        // The Population does not have to be empty to add Genomes
        bool spawn(Genome *g, int size);

    public:

        std::vector<Organism *> organisms; //The organisms in the Population

        std::vector<Species *> species;  // Species in the Population. Note that the species should comprise all the genomes

        // ******* Member variables used during reproduction *******
        std::vector<Innovation *> innovations;  // For holding the genetic innovations of the newest generation
        int cur_node_id;  //Current label number available
        double cur_innov_num;

        int last_species;  //The highest species number

        // Separate the Organisms into species
        bool speciate();

        bool print_to_file_by_species(char *filename);

        // Run verify on all Genomes in this Population (Debugging)
        bool verify();

        // *** Real-time methods ***

        // Places the organisms in species in order from best to worst fitness
        bool rank_within_species();

        // Estimates average fitness for all existing species
        void estimate_all_averages();

        // Probabilistically choose a species to reproduce
        // Note that this method is effectively real-time fitness sharing in that the
        // species will tend to produce offspring in an amount proportional
        // to their average fitness, which approximates the generational
        // method of producing the next generation of the species en masse
        // based on its average (shared) fitness.
        Species *choose_parent_species();

        //Remove a species from the species list (sometimes called by remove_worst when a species becomes empty)
        bool remove_species(Species *spec);

        // Removes worst member of population that has been around for a minimum amount of time and returns
        // a pointer to the Organism that was removed (note that the pointer will not point to anything at all,
        // since the Organism it was pointing to has been deleted from memory)
        Organism *remove_worst();

        //KEN: New 2/17/04
        //This method takes an Organism and reassigns what Species it belongs to
        //It is meant to be used so that we can reasses where Organisms should belong
        //as the speciation threshold changes.
        void reassign_species(Organism *org);

        //Move an Organism from one Species to another (called by reassign_species)
        void switch_species(Organism *org, Species *orig_species, Species *new_species);

        // Construct off of a single spawning Genome
        Population(Genome *g, int size);

        // Construct off of a single spawning Genome without mutation
        Population(Genome *g, int size, float power);

        //MSC Addition
        // Construct off of a vector of genomes with a mutation rate of "power"
        Population(std::vector<Genome *> genomeList, float power);

        bool clone(Genome *g, int size, float power);

        //// Special constructor to create a population of random topologies
        //// uses Genome(int i, int o, int n,int nmax, bool r, double linkprob)
        //// See the Genome constructor for the argument specifications
        //Population(int size,int i,int o, int nmax, bool r, double linkprob);

        // Construct off of a file of Genomes
        Population(const char *filename);

        // It can delete a Population in two ways:
        //    -delete by killing off the species
        //    -delete by killing off the organisms themselves (if not speciated)
        // It does the latter if it sees the species list is empty
        ~Population();


    };

}
#endif
