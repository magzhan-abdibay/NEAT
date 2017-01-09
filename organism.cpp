#include "organism.h"

using namespace NEAT;

Organism::Organism(double fit, Genome *g, int gen, const char *md) {
    fitness = fit;
    orig_fitness = fitness;
    gnome = g;
    net = gnome->genesis(gnome->genome_id);
    species = 0;  //Start it in no Species
    expected_offspring = 0;
    generation = gen;
    eliminate = false;
    error = 0;
    winner = false;
    champion = false;
    super_champ_offspring = 0;

    // If md is null, then we don't have metadata, otherwise we do have metadata so copy it over
    if (md == 0) {
        strcpy(metadata, "");
    } else {
        strncpy(metadata, md, 128);
    }

    time_alive = 0;

    //DEBUG vars
    pop_champ = false;
    pop_champ_child = false;
    high_fit = 0;
    mut_struct_baby = 0;
    mate_baby = 0;

    modified = true;
}

Organism::Organism(const Organism &org) {
    fitness = org.fitness;
    orig_fitness = org.orig_fitness;
    gnome = new Genome(*(org.gnome));    // Associative relationship
    //gnome = org.gnome->duplicate(org.gnome->genome_id);
    net = new Network(*(org.net)); // Associative relationship
    species = org.species;    // Delegation relationship
    expected_offspring = org.expected_offspring;
    generation = org.generation;
    eliminate = org.eliminate;
    error = org.error;
    winner = org.winner;
    champion = org.champion;
    super_champ_offspring = org.super_champ_offspring;

    strcpy(metadata, org.metadata);
    //printf("copying %s did it work? %s", org.metadata, metadata);

    time_alive = org.time_alive;
    pop_champ = org.pop_champ;
    pop_champ_child = org.pop_champ_child;
    high_fit = org.high_fit;
    mut_struct_baby = org.mut_struct_baby;
    mate_baby = org.mate_baby;

    modified = false;
}

Organism::~Organism() {
    delete net;
    delete gnome;
}

bool Organism::print_to_file(char *filename) {

    std::ofstream oFile(filename);

    return write_to_file(oFile);
}

bool Organism::write_to_file(std::ostream &outFile) {

    char tempbuf2[1024];
    if (modified) {
        sprintf(tempbuf2, "/* Organism #%d Fitness: %f Time: %d */\n", (gnome)->genome_id, fitness, time_alive);
    } else {
        sprintf(tempbuf2, "/* %s */\n", metadata);
    }
    outFile << tempbuf2;
    gnome->print_to_file(outFile);
    return 1;
}

bool NEAT::order_orgs(Organism *x, Organism *y) {
    return (x)->fitness > (y)->fitness;
}

