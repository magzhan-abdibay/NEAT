#include "population.h"
#include <iostream>
#include <sstream>

using namespace NEAT;

extern int NEAT::time_alive_minimum;

Population::Population(Genome *g, int size) {
    spawn(g, size);
}

Population::Population(Genome *g, int size, float power) {
    clone(g, size, power);
}

//MSC Addition
//Added the ability for a population to be spawned
//off of a vector of Genomes.  Useful when converging.
Population::Population(std::vector<Genome *> genomeList, float power) {

    Genome *new_genome;
    Organism *new_organism;

    //Create size copies of the Genome
    //Start with perturbed linkweights
    for (std::vector<Genome *>::iterator iter = genomeList.begin(); iter != genomeList.end(); ++iter) {

        new_genome = (*iter);
        if (power > 0)
            new_genome->mutate_link_weights(power, 1.0, GAUSSIAN);
        //new_genome->mutate_link_weights(1.0,1.0,COLDGAUSSIAN);
        new_genome->randomize_traits();
        new_organism = new Organism(0.0, new_genome, 1);
        organisms.push_back(new_organism);
    }

    //Keep a record of the innovation and node number we are on
    cur_node_id = new_genome->get_last_node_id();
    cur_innov_num = new_genome->get_last_gene_innovnum();

    //Separate the new Population into species
    speciate();
}

Population::Population(const char *filename) {

    char curword[128];  //max word size of 128 characters
    char curline[1024]; //max line size of 1024 characters

    Genome *new_genome;

    cur_node_id = 0;
    cur_innov_num = 0.0;

    std::ifstream iFile(filename);
    if (!iFile) {
        printf("Can't open genomes file for input");
        return;
    } else {
        bool md = false;
        char metadata[128];
        //Loop until file is finished, parsing each line
        while (!iFile.eof()) {
            iFile.getline(curline, sizeof(curline));
            std::stringstream ss(curline);
            //strcpy(curword, NEAT::getUnit(curline, 0, delimiters));
            ss >> curword;
            //std::cout << curline << std::endl;

            //Check for next
            if (strcmp(curword, "genomestart") == 0) {
                //strcpy(curword, NEAT::getUnit(curline, 1, delimiters));
                //int idcheck = atoi(curword);

                int idcheck;
                ss >> idcheck;

                // If there isn't metadata, set metadata to ""
                if (!md) {
                    strcpy(metadata, "");
                }
                md = false;

                new_genome = new Genome(idcheck, iFile);
                organisms.push_back(new Organism(0, new_genome, 1, metadata));
                if (cur_node_id < (new_genome->get_last_node_id()))
                    cur_node_id = new_genome->get_last_node_id();

                if (cur_innov_num < (new_genome->get_last_gene_innovnum()))
                    cur_innov_num = new_genome->get_last_gene_innovnum();
            } else if (strcmp(curword, "/*") == 0) {
                // New metadata possibly, so clear out the metadata
                strcpy(metadata, "");
                //strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
                ss >> curword;

                while (strcmp(curword, "*/") != 0) {
                    // If we've started to form the metadata, put a space in the front
                    if (md) {
                        strncat(metadata, " ", 128 - strlen(metadata));
                    }

                    // Append the next word to the metadata, and say that there is metadata
                    strncat(metadata, curword, 128 - strlen(metadata));
                    md = true;

                    //strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
                    ss >> curword;
                }
            }
            //Ignore comments - they get printed to screen
            //else if (strcmp(curword,"/*")==0) {
            //	iFile>>curword;
            //	while (strcmp(curword,"*/")!=0) {
            //cout<<curword<<" ";
            //		iFile>>curword;
            //	}
            //	cout<<endl;

            //}
            //Ignore comments surrounded by - they get printed to screen
        }

        iFile.close();

        speciate();

    }
}


Population::~Population() {

    std::vector<Species *>::iterator curspec;
    std::vector<Organism *>::iterator curorg;
    //std::vector<Generation_viz*>::iterator cursnap;

    if (species.begin() != species.end()) {
        for (curspec = species.begin(); curspec != species.end(); ++curspec) {
            delete (*curspec);
        }
    } else {
        for (curorg = organisms.begin(); curorg != organisms.end(); ++curorg) {
            delete (*curorg);
        }
    }

    for (std::vector<Innovation *>::iterator iter = innovations.begin(); iter != innovations.end(); ++iter)
        delete *iter;

    //Delete the snapshots
    //		for(cursnap=generation_snapshots.begin();cursnap!=generation_snapshots.end();++cursnap) {
    //			delete (*cursnap);
    //		}
}

bool Population::verify() {
    std::vector<Organism *>::iterator curorg;

    bool verification;

    for (curorg = organisms.begin(); curorg != organisms.end(); ++curorg) {
        verification = ((*curorg)->gnome)->verify();
    }

    return verification;
}

bool Population::clone(Genome *g, int size, float power) {
    int count;
    Genome *new_genome;
    Organism *new_organism;

    new_genome = g->duplicate(1);
    new_organism = new Organism(0.0, new_genome, 1);
    organisms.push_back(new_organism);

    //Create size copies of the Genome
    //Start with perturbed linkweights
    for (count = 2; count <= size; count++) {
        new_genome = g->duplicate(count);
        if (power > 0)
            new_genome->mutate_link_weights(power, 1.0, GAUSSIAN);

        new_genome->randomize_traits();
        new_organism = new Organism(0.0, new_genome, 1);
        organisms.push_back(new_organism);
    }

    //Keep a record of the innovation and node number we are on
    cur_node_id = new_genome->get_last_node_id();
    cur_innov_num = new_genome->get_last_gene_innovnum();

    //Separate the new Population into species
    speciate();

    return true;
}

bool Population::spawn(Genome *g, int size) {
    int count;
    Genome *new_genome;
    Organism *new_organism;

    //Create size copies of the Genome
    //Start with perturbed linkweights
    for (count = 1; count <= size; count++) {
        //cout<<"CREATING ORGANISM "<<count<<endl;

        new_genome = g->duplicate(count);
        //new_genome->mutate_link_weights(1.0,1.0,GAUSSIAN);
        new_genome->mutate_link_weights(1.0, 1.0, COLDGAUSSIAN);
        new_genome->randomize_traits();
        new_organism = new Organism(0.0, new_genome, 1);
        organisms.push_back(new_organism);
    }

    //Keep a record of the innovation and node number we are on
    cur_node_id = new_genome->get_last_node_id();
    cur_innov_num = new_genome->get_last_gene_innovnum();

    //Separate the new Population into species
    speciate();

    return true;
}

bool Population::speciate() {
    std::vector<Organism *>::iterator curorg;  //For stepping through Population
    std::vector<Species *>::iterator curspecies; //Steps through species
    Organism *comporg = 0;  //Organism for comparison
    Species *newspecies; //For adding a new species

    int counter = 0; //Species counter

    //Step through all existing organisms
    for (curorg = organisms.begin(); curorg != organisms.end(); ++curorg) {

        //For each organism, search for a species it is compatible to
        curspecies = species.begin();
        if (curspecies == species.end()) {
            //Create the first species
            newspecies = new Species(++counter);
            species.push_back(newspecies);
            newspecies->add_Organism(*curorg);  //Add the current organism
            (*curorg)->species = newspecies;  //Point organism to its species
        } else {
            comporg = (*curspecies)->first();
            while ((comporg != 0) &&
                   (curspecies != species.end())) {

                if ((((*curorg)->gnome)->compatibility(comporg->gnome)) < NEAT::compat_threshold) {

                    //Found compatible species, so add this organism to it
                    (*curspecies)->add_Organism(*curorg);
                    (*curorg)->species = (*curspecies);  //Point organism to its species
                    comporg = 0;  //Note the search is over
                } else {

                    //Keep searching for a matching species
                    ++curspecies;
                    if (curspecies != species.end())
                        comporg = (*curspecies)->first();
                }
            }

            //If we didn't find a match, create a new species
            if (comporg != 0) {
                newspecies = new Species(++counter);
                species.push_back(newspecies);
                newspecies->add_Organism(*curorg);  //Add the current organism
                (*curorg)->species = newspecies;  //Point organism to its species
            }

        } //end else

    } //end for

    last_species = counter;  //Keep track of highest species

    return true;
}

bool Population::print_to_file_by_species(char *filename) {

    std::vector<Species *>::iterator curspecies;

    std::ofstream outFile(filename, std::ios::out);

    //Make sure it worked
    if (!outFile) {
        std::cerr << "Can't open " << filename << " for output" << std::endl;
        return false;
    }


    //Step through the Species and print them to the file
    for (curspecies = species.begin(); curspecies != species.end(); ++curspecies) {
        (*curspecies)->print_to_file(outFile);
    }

    outFile.close();

    return true;

}


bool Population::rank_within_species() {
    std::vector<Species *>::iterator curspecies;

    //Add each Species in this generation to the snapshot
    for (curspecies = species.begin(); curspecies != species.end(); ++curspecies) {
        (*curspecies)->rank();
    }

    return true;
}

void Population::estimate_all_averages() {
    std::vector<Species *>::iterator curspecies;

    for (curspecies = species.begin(); curspecies != species.end(); ++curspecies) {
        (*curspecies)->estimate_average();
    }

}

Species *Population::choose_parent_species() {

    double total_fitness = 0;
    std::vector<Species *>::iterator curspecies;
    double marble; //The roulette marble
    double spin; //Spins until the marble reaches its chosen point

    //Keeping species diverse
    //This commented out code forces the system to aim for
    // num_species species at all times, enforcing diversity
    //This tinkers with the compatibility threshold, which
    // normally would be held constant

    //Use the roulette method to choose the species

    //Sum all the average fitness estimates of the different species
    //for the purposes of the roulette
    for (curspecies = species.begin(); curspecies != species.end(); ++curspecies) {
        total_fitness += (*curspecies)->average_est;
    }

    marble = randfloat() * total_fitness;
    curspecies = species.begin();
    spin = (*curspecies)->average_est;
    while (spin < marble) {
        ++curspecies;

        //Keep the wheel spinning
        spin += (*curspecies)->average_est;
    }
    //Finished roulette

    //Return the chosen species
    return (*curspecies);
}

bool Population::remove_species(Species *spec) {
    std::vector<Species *>::iterator curspec;

    curspec = species.begin();
    while ((curspec != species.end()) &&
           ((*curspec) != spec))
        ++curspec;

    if (curspec == species.end()) {
        //   cout<<"ALERT: Attempt to remove nonexistent Species from Population"<<endl;
        return false;
    } else {
        species.erase(curspec);
        return true;
    }
}

Organism *Population::remove_worst() {

    double adjusted_fitness;
    double min_fitness = 999999;
    std::vector<Organism *>::iterator curorg;
    Organism *org_to_kill = 0;
    std::vector<Organism *>::iterator deadorg;
    Species *orgs_species; //The species of the dead organism

    //Make sure the organism is deleted from its species and the population

    //First find the organism with minimum *adjusted* fitness
    for (curorg = organisms.begin(); curorg != organisms.end(); ++curorg) {
        adjusted_fitness = ((*curorg)->fitness) / ((*curorg)->species->organisms.size());
        if ((adjusted_fitness < min_fitness) &&
            (((*curorg)->time_alive) >= NEAT::time_alive_minimum)) {
            min_fitness = adjusted_fitness;
            org_to_kill = (*curorg);
            deadorg = curorg;
            orgs_species = (*curorg)->species;
        }
    }

    if (org_to_kill) {

        //Remove the organism from its species and the population
        orgs_species->remove_org(org_to_kill);  //Remove from species
        delete org_to_kill;  //Delete the organism itself
        organisms.erase(deadorg); //Remove from population list

        //Did the species become empty?
        if ((orgs_species->organisms.size()) == 0) {

            remove_species(orgs_species);
            delete orgs_species;
        }
            //If not, re-estimate the species average after removing the organism
        else {
            orgs_species->estimate_average();
        }
    }

    return org_to_kill;
}

//KEN: New 2/17/04
//This method takes an Organism and reassigns what Species it belongs to
//It is meant to be used so that we can reasses where Organisms should belong
//as the speciation threshold changes.
void Population::reassign_species(Organism *org) {

    Organism *comporg;
    bool found = false;  //Note we don't really need this flag but it
    //might be useful if we change how this function works
    std::vector<Species *>::iterator curspecies;
    Species *newspecies;

    curspecies = (species).begin();

    comporg = (*curspecies)->first();
    while ((curspecies != (species).end()) && (!found)) {
        if (comporg == 0) {
            //Keep searching for a matching species
            ++curspecies;
            if (curspecies != (species).end())
                comporg = (*curspecies)->first();
        } else if (((org->gnome)->compatibility(comporg->gnome)) < NEAT::compat_threshold) {
            //If we found the same species it's already in, return 0
            if ((*curspecies) == org->species) {
                found = true;
                break;
            }
                //Found compatible species
            else {
                switch_species(org, org->species, (*curspecies));
                found = true;  //Note the search is over
            }
        } else {
            //Keep searching for a matching species
            ++curspecies;
            if (curspecies != (species).end())
                comporg = (*curspecies)->first();
        }
    }

    //If we didn't find a match, create a new species, move the org to
    // that species, check if the old species is empty,
    //re-estimate averages, and return 0
    if (!found) {

        //Create a new species for the org
        newspecies = new Species(++(last_species), true);
        (species).push_back(newspecies);

        switch_species(org, org->species, newspecies);
    }

}

//Move an Organism from one Species to another
void Population::switch_species(Organism *org, Species *orig_species, Species *new_species) {

    //Remove organism from the species we want to remove it from
    orig_species->remove_org(org);

    //Add the organism to the new species it is being moved to
    new_species->add_Organism(org);
    org->species = new_species;

    //KEN: Delete orig_species if empty, and remove it from pop
    if ((orig_species->organisms.size()) == 0) {

        remove_species(orig_species);
        delete orig_species;

        //Re-estimate the average of the species that now has a new member
        new_species->estimate_average();
    }
        //If not, re-estimate the species average after removing the organism
        // AND the new species with the new member
    else {
        orig_species->estimate_average();
        new_species->estimate_average();
    }
}