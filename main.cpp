#include <iostream>
#include <random>
#include <iostream>
#include <string>
#include <fstream>
#include "Polymer.h"
#include "System.h"
#include "Tests.h"

using namespace std;

int main() {

    Tests t;
    t.run();

    double seed = 200;
    mt19937 gen(seed);

    //Initialisations
    double k = 1;//Polymerisation rate
    double k0 = 1;//Binding rate
    double G_bb = -16;//Backbone forming free energy
    double G_spec = 2;//Specific bond forming free energy
    double G_gen = -2;//Generic bond forming free energy
    double M_eff = 100;//Effective concentration of monomers in zipping

    vector<double> rates({k, k0});
    vector<double> energies({G_bb, G_spec, G_gen, M_eff});

    int monomers_family_zero = 1000;
    int monomers_family_one = 1000;

    vector<int> free_monomers({monomers_family_zero, monomers_family_one});

    int template_length = 6;

    Polymer * template_polymer = new Polymer(-1, template_length, 0);

    System * system = new System(rates, energies, free_monomers, template_polymer);

    int count = 0;
    int transition_limit = 10000;

    while(count<transition_limit){
        system->chooseTransition(gen());
        count ++;
    }
/*
    for(auto & cong : system->conglomerates){
        cout << endl << "Conglomerate: " << cong->index << endl;
        cout << "Polymers:" << endl;
        for(auto & pol : cong->polymers) {
            cout << "Index:" << pol->index << " Length: " << pol->length << endl;
        }
        cout << endl << "Connections:" << endl;
        for(auto & con : cong->connections){
            cout << "Polymer: " << con->polymers_in_connection[0]->index << " Index: " << con->indexes[0] << endl;
            cout << "Polymer: " << con->polymers_in_connection[1]->index << " Index: " << con->indexes[1] << endl << endl;
        }
    }*/

    ofstream fw("/Users/nicksimpson/PycharmProjects/MyProject/input.txt", ofstream::out);

    if(fw.is_open()){
        //Creating all the nodes and joining polymers
        for(int cong = 0; cong<system->conglomerates.size(); cong ++) {
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
            for(int poly = 0; poly<system->conglomerates[cong]->polymer_connections.size(); poly++){
                for(int mon=0; mon<system->conglomerates[cong]->polymer_connections[poly].size(); mon++){
                    if(!system->conglomerates[cong]->polymer_connections[poly][mon].empty()){
                        if(is_first){
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
            for(int poly = 0; poly<system->conglomerates[cong]->polymer_connections.size(); poly++){
                for(int mon=0; mon<system->conglomerates[cong]->polymer_connections[poly].size(); mon++){
                    if(!system->conglomerates[cong]->polymer_connections[poly][mon].empty()){
                        if(is_first){
                            //We need to not have a comma before the first connection
                            fw << '(' << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[0];
                            fw << ',' << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[1] << ")";
                            is_first = false;
                        } else {
                            fw << ", (" << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[0];
                            fw << ',' << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[1] << ")";
                        }
                    }
                }
            }
            fw << ']' << endl;
        }

        fw.close();
    }

    delete system;
    return 0;
}
