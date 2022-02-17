#include <iostream>
#include <random>
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
    double G_spec = -4;//Specific bond forming free energy
    double G_gen = -8;//Generic bond forming free energy
    double M_eff = 100;//Effective concentration of monomers in zipping

    vector<double> rates({k, k0});
    vector<double> energies({G_bb, G_spec, G_gen, M_eff});

    int monomers_family_zero = 100;
    int monomers_family_one = 100;

    vector<int> free_monomers({monomers_family_zero, monomers_family_one});

    int template_length = 6;

    Polymer * template_polymer = new Polymer(-1, template_length, 0);

    System * system = new System(rates, energies, free_monomers, template_polymer);

    int count = 0;
    int transition_limit = 100000;

    while(count<transition_limit){
        cout << endl << "-----------------------------" << endl;
        system->chooseTransition(gen());
        count ++;
    }

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
    }

    delete system;
    return 0;
}
