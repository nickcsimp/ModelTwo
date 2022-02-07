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
    mt19937 gen();

    vector<Conglomerate *> conglomerates;

    vector<vector<int>> external_sites; //Family<Conglomerate>
    vector<int> head_unbinding_connections;
    vector<int> head_binding_connections;
    vector<int> tail_unbinding_connections;
    vector<int> tail_binding_connections;
    vector<int> connected_neighbours;
    vector<int> unconnected_neighbours;

    vector<double> external_sites_rate;
    vector<double> head_unbinding_connections_rate;
    vector<double> head_binding_connections_rate;
    vector<double> tail_unbinding_connections_rate;
    vector<double> tail_binding_connections_rate;
    vector<double> connected_neighbours_rate;
    vector<double> unconnected_neighbours_rate;

    void updateRates(int cong);
    void updateBonds(int cong);
    void chooseTransition(double seed);
    void mergeConglomerates();


};


#endif //MODELTWO_SYSTEM_H
