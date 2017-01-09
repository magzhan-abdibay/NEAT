#include "gene.h"

#include <iostream>
#include <sstream>

using namespace NEAT;

Gene::Gene(double w, NNode *inode, NNode *onode, bool recur, double innov, double mnum) {
    lnk = new Link(w, inode, onode, recur);
    innovation_num = innov;
    mutation_num = mnum;

    enable = true;

    frozen = false;
}


Gene::Gene(Trait *tp, double w, NNode *inode, NNode *onode, bool recur, double innov, double mnum) {
    lnk = new Link(tp, w, inode, onode, recur);
    innovation_num = innov;
    mutation_num = mnum;

    enable = true;

    frozen = false;
}

Gene::Gene(Gene *g, Trait *tp, NNode *inode, NNode *onode) {
    lnk = new Link(tp, (g->lnk)->weight, inode, onode, (g->lnk)->isRecurrent);
    innovation_num = g->innovation_num;
    mutation_num = g->mutation_num;
    enable = g->enable;

    frozen = g->frozen;
}

Gene::Gene(const char *argline, std::vector<Trait *> &traits, std::vector<NNode *> &nodes) {
    //Gene parameter holders
    int traitnum;
    int inodenum;
    int onodenum;
    NNode *inode;
    NNode *onode;
    double weight;
    int recur;
    Trait *traitptr;

    std::vector<Trait *>::iterator curtrait;
    std::vector<NNode *>::iterator curnode;

    //Get the gene parameters

    std::stringstream ss(argline);
    ss >> traitnum >> inodenum >> onodenum >> weight >> recur >> innovation_num >> mutation_num >> enable;

    frozen = false; //TODO: MAYBE CHANGE

    //Get a pointer to the linktrait
    if (traitnum == 0) traitptr = 0;
    else {
        curtrait = traits.begin();
        while (((*curtrait)->trait_id) != traitnum)
            ++curtrait;
        traitptr = (*curtrait);
    }

    //Get a pointer to the input node
    curnode = nodes.begin();
    while (((*curnode)->node_id) != inodenum)
        ++curnode;
    inode = (*curnode);

    //Get a pointer to the output node
    curnode = nodes.begin();
    while (((*curnode)->node_id) != onodenum)
        ++curnode;
    onode = (*curnode);

    lnk = new Link(traitptr, weight, inode, onode, recur);
}

Gene::Gene(const Gene &gene) {
    innovation_num = gene.innovation_num;
    mutation_num = gene.mutation_num;
    enable = gene.enable;
    frozen = gene.frozen;

    lnk = new Link(*gene.lnk);
}

Gene::~Gene() {
    delete lnk;
}


void Gene::print_to_file(std::ofstream &outFile) {
    outFile << "gene ";
    //Start off with the trait number for this gene
    if ((lnk->linkTrait) == 0) outFile << "0 ";
    else outFile << ((lnk->linkTrait)->trait_id) << " ";
    outFile << (lnk->inNode)->node_id << " ";
    outFile << (lnk->outNode)->node_id << " ";
    outFile << (lnk->weight) << " ";
    outFile << (lnk->isRecurrent) << " ";
    outFile << innovation_num << " ";
    outFile << mutation_num << " ";
    outFile << enable << std::endl;
}


void Gene::print_to_file(std::ostream &outFile) {
    outFile << "gene ";
    //outFile.write(5, "gene ");

    //Start off with the trait number for this gene
    if ((lnk->linkTrait) == 0) {
        outFile << "0 ";
        //outFile.write(2, "0 ");
    } else {
        outFile << ((lnk->linkTrait)->trait_id) << " ";
    }
    outFile << (lnk->inNode)->node_id << " ";
    outFile << (lnk->outNode)->node_id << " ";
    outFile << (lnk->weight) << " ";
    outFile << (lnk->isRecurrent) << " ";
    outFile << innovation_num << " ";
    outFile << mutation_num << " ";
    outFile << enable << std::endl;
}
