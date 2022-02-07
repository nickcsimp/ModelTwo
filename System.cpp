//
// Created by Nicholas Simpson on 07/02/2022.
//

#include "System.h"


void System::updateBonds(int cong){
    external_sites[0][cong] = conglomerates[cong]->available_free_sites_list[0].size();
    external_sites[1][cong] = conglomerates[cong]->available_free_sites_list[1].size();

    head_unbinding_connections[cong] = conglomerates[cong]->head_unbinding_list.size();
    head_binding_connections[cong] = conglomerates[cong]->head_binding_list.size();

    tail_unbinding_connections[cong] = conglomerates[cong]->tail_unbinding_list.size();
    tail_binding_connections[cong] = conglomerates[cong]->tail_binding_list.size();

    connected_neighbours[cong] = conglomerates[cong]->connected_neighbours_list.size();
    unconnected_neighbours[cong] = conglomerates[cong]->unconnected_neighbours_list.size();
}

void System::updateRates(int cong){
    /*TODO:
     * Think really hard about external connections
     * Do we need rate vectors or just do the maths on the fly?
     * Need to store rates, initialise in system initialisation
     * Make sure polymers not being saved in system doesnt fuck us over
     */

}

void System::chooseTransition(double seed){
    mt19937 gen(seed);

    (gen()/mt19937::max())
}