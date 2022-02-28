//
// Created by Nicholas Simpson on 07/02/2022.
//

#ifndef MODELTWO_SYSTEM_H
#define MODELTWO_SYSTEM_H


#include "Conglomerate.h"

class System {
public:
    int conglomerate_index;
    int polymer_index;

    vector<Conglomerate *> conglomerates;

    double total_rate;

    vector<vector<int>> external_sites; //Family<Conglomerate<number of sites>>
    vector<int> total_external_sites; //Family<number of sites>
    double external_connection_rate;

    /*Transition numbers:
     * 0 = head unbinding
     * 1 = head binding
     * 2 = tail unbinding
     * 3 = tail binding
     * 4 = neighbours unbind
     * 5 = neighbours bind
     */
    vector<double> transition_rates; //transition<total_rate>
    vector<vector<int>> conglomerate_rates; //transition<conglomerate<number of bonds>>
    vector<int> total_connections; //transition<total number of bonds>

    double k;//Polymerisation rate
    double k0;//Binding rate
    double G_bb;//Backbone forming free energy
    double G_spec;//Specific bond forming free energy
    double G_gen;//Generic bond forming free energy
    double M_eff;//Effective concentration of monomers in zipping


    System(vector<double> rates, vector<double> energies, vector<int> free_monomers, Polymer * template_polymer);

    void updateRates(int cong);
    void chooseTransition(double seed);
    void removeConglomerate(int cong);
    void addConglomerate(Conglomerate * new_cong);


};


#endif //MODELTWO_SYSTEM_H
