#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include "Polymer.h"
#include "System.h"
#include "Tests.h"
#include "settings.h"

using namespace std;

bool set_template_indestructible;
bool set_monomer_count_is_constant; //Only useful if no_rebinding = true
bool set_no_rebinding;
double set_seed;
double set_k;//Polymerisation rate
double set_k0;//Binding rate
double set_G_bb;//Backbone forming free energy
double set_G_spec;//Specific bond forming free energy
double set_G_gen;//Generic bond forming free energy
double set_M_eff;//Effective concentration of monomers in zipping
int set_monomers_family_zero;
int set_monomers_family_one;
int set_template_length;
int set_transition_limit;

int main() {
/*
    Tests tests;
    tests.run();
*/

/*
 * TODO
 * Make an input csv set initialisers
 * Make outputs relative
 * Merge the python and c into one seamless package
 */
    set_template_indestructible = true;
    set_monomer_count_is_constant = false;
    set_no_rebinding = false;
    set_seed = 200;
    set_k = 1;
    set_k0 = 1;
    set_G_bb = -15;
    set_G_spec = -4;
    set_G_gen = -5;
    set_M_eff = 100;
    set_monomers_family_zero = 0;
    set_monomers_family_one = 100;
    set_template_length = 3;
    set_transition_limit = 1000;

    double seed = set_seed;
    mt19937 gen(seed);

    System * system = new System();

    int count = 0;
    int transition_limit = set_transition_limit;
    bool transitions_possible = true;
    ofstream f_hist("/Users/nicksimpson/PycharmProjects/MyProject/histogram.txt", ofstream::out);

    while(count<transition_limit && transitions_possible){
        transitions_possible = system->chooseTransition(gen());
        count ++;
        if(f_hist.is_open()){
            f_hist << '[';
            for(int i = 0; i< system->lengths.size()-1; i++){
                f_hist << system->lengths[i] << ", ";
            }
            f_hist << system->lengths[system->lengths.size()-1] << ']' << "\n";
        }

    }

    f_hist.close();

    ofstream fw("/Users/nicksimpson/PycharmProjects/MyProject/input.txt", ofstream::out);

    if(fw.is_open()){
        //Creating all the nodes and joining polymers
        for(int cong = 0; cong<system->conglomerates.size(); cong ++) {
            //We only want to view the things that are interesting
            //If there is a single monomer or a dimer then we ignore it
            if (system->conglomerates[cong]->polymers.size() == 1 &&
                system->conglomerates[cong]->polymers[0]->length == 1) {
                //TODO make a monomer image?
                //Monomer
            } else if (system->conglomerates[cong]->polymers.size() == 2 &&
                       system->conglomerates[cong]->polymers[0]->length == 1 &&
                       system->conglomerates[cong]->polymers[1]->length == 1) {
                //Dimer
                //TODO make a dimer image?
            } else {

                int mon_index = -1;
                int pol_index = -1;
                bool is_first = true;
                for (int i = 0; i < system->conglomerates[cong]->polymers.size(); i++) {
                    system->conglomerates[cong]->polymers[i]->index = ++pol_index;
                    if (is_first) {
                        fw << "([";
                        is_first = false;
                    } else {
                        fw << ", [";
                    }
                    for (int j = 0; j < system->conglomerates[cong]->polymers[i]->length - 1; j++) {
                        fw << ++mon_index << ',';
                    }
                    fw << ++mon_index << ']';
                }
                fw << ')' << "\n";

                is_first = true;
                fw << '(';
                for (auto &poly: system->conglomerates[cong]->polymers) {
                    if (is_first) {
                        fw << poly->family;
                        is_first = false;
                    } else {
                        fw << ", " << poly->family;
                    }
                }
                fw << ')' << "\n";



                //Creating inter-polymer connection polymers
                is_first = true;
                fw << '[';
                for (int poly = 0; poly < system->conglomerates[cong]->polymer_connections.size(); poly++) {
                    for (int mon = 0; mon < system->conglomerates[cong]->polymer_connections[poly].size(); mon++) {
                        if (!system->conglomerates[cong]->polymer_connections[poly][mon].empty()) {
                            if (is_first) {
                                //We need to not have a comma before the first connection
                                fw << "(Polymers["
                                   << system->conglomerates[cong]->polymer_connections[poly][mon][0]->polymers_in_connection[0]->index;
                                fw << "], Polymers["
                                   << system->conglomerates[cong]->polymer_connections[poly][mon][0]->polymers_in_connection[1]->index;
                                fw << "])";
                                is_first = false;
                            } else {
                                fw << ", (Polymers["
                                   << system->conglomerates[cong]->polymer_connections[poly][mon][0]->polymers_in_connection[0]->index;
                                fw << "], Polymers["
                                   << system->conglomerates[cong]->polymer_connections[poly][mon][0]->polymers_in_connection[1]->index;
                                fw << "])";
                            }
                        }
                    }
                }
                fw << ']' << "\n";
                //Creating inter-polymer connection indexes
                is_first = true;
                fw << '[';
                for (int poly = 0; poly < system->conglomerates[cong]->polymer_connections.size(); poly++) {
                    for (int mon = 0; mon < system->conglomerates[cong]->polymer_connections[poly].size(); mon++) {
                        if (!system->conglomerates[cong]->polymer_connections[poly][mon].empty()) {
                            if (is_first) {
                                //We need to not have a comma before the first connection
                                fw << '(' << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[0];
                                fw << ',' << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[1]
                                   << ")";
                                is_first = false;
                            } else {
                                fw << ", ("
                                   << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[0];
                                fw << ',' << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[1]
                                   << ")";
                            }
                        }
                    }
                }
                fw << ']' << endl;
            }
        }

        fw.close();
    }

    delete system;
    return 0;
}
