//
// Created by Nicholas Simpson on 07/02/2022.
//

#include "System.h"


void System::updateRates(int cong){
    total_rate = 0;

    total_external_sites[0] = total_external_sites[0] - external_sites[0][cong];
    total_external_sites[1] = total_external_sites[1] - external_sites[1][cong];
    external_sites[0][cong] = conglomerates[cong]->available_free_sites_list[0].size();
    external_sites[1][cong] = conglomerates[cong]->available_free_sites_list[1].size();
    total_external_sites[0] = total_external_sites[0] + external_sites[0][cong];
    total_external_sites[1] = total_external_sites[1] + external_sites[1][cong];

    //TODO: is there a better way to update this rather than starting from nought
    external_connection_rate = 0;
    //For every conglomerate, we look at the number of family 0 free sites
    //These can bind with any other conglomerate's family 1 free sites
    //Therefore two loops to cover entire system
    for(int i=0; i<conglomerates.size(); i++) { //Loop conglomerates looking at family 0
        for(int j=0; j<conglomerates.size(); j++) { //Loop conglomerates looking at family 1
            if(i!=j){ //Conglomerate can't bind within itself here
                external_connection_rate = external_connection_rate + external_sites[0][i]*external_sites[1][j]*k0;
            }
        }
    }
    total_rate = total_rate + external_connection_rate;

    transition_rates[0] = transition_rates[0] - conglomerate_rates[0][cong]*(k0*exp(G_spec+G_gen));
    conglomerate_rates[0][cong] = conglomerates[cong]->head_unbinding_list.size();
    transition_rates[0] = transition_rates[0] + conglomerate_rates[0][cong]*(k0*exp(G_spec+G_gen));
    total_rate = total_rate + transition_rates[0];

    transition_rates[1] = transition_rates[1] - conglomerate_rates[1][cong]*(k0*M_eff);
    conglomerate_rates[1][cong] = conglomerates[cong]->head_binding_list.size();
    transition_rates[1] = transition_rates[1] + conglomerate_rates[1][cong]*(k0*M_eff);
    total_rate = total_rate + transition_rates[1];

    transition_rates[2] = transition_rates[2] - conglomerate_rates[2][cong]*(k0*exp(G_spec));
    conglomerate_rates[2][cong] = conglomerates[cong]->tail_unbinding_list.size();
    transition_rates[2] = transition_rates[2] + conglomerate_rates[2][cong]*(k0*exp(G_spec));
    total_rate = total_rate + transition_rates[2];

    transition_rates[3] = transition_rates[3] - conglomerate_rates[3][cong]*(k0*M_eff);
    conglomerate_rates[3][cong] = conglomerates[cong]->tail_binding_list.size();
    transition_rates[3] = transition_rates[3] + conglomerate_rates[3][cong]*(k0*M_eff);
    total_rate = total_rate + transition_rates[3];

    transition_rates[4] = transition_rates[4] - conglomerate_rates[4][cong]*(k*exp(G_bb-G_gen));
    conglomerate_rates[4][cong] = conglomerates[cong]->connected_neighbours_list.size();
    transition_rates[4] = transition_rates[4] + conglomerate_rates[4][cong]*(k*exp(G_bb-G_gen));
    total_rate = total_rate + transition_rates[4];

    transition_rates[5] = transition_rates[5] - conglomerate_rates[5][cong]*k;
    conglomerate_rates[5][cong] = conglomerates[cong]->unconnected_neighbours_list.size();
    transition_rates[5] = transition_rates[5] + conglomerate_rates[5][cong]*k;
    total_rate = total_rate + transition_rates[5];
}

    /*TODO:
     * Think really hard about external connections
     * Make sure polymers not being saved in system doesnt fuck us over
     */


void System::chooseTransition(double seed){
    mt19937 gen(seed);

    //First choose a transition
    //External rates is weird, so we see if that is chosen first

    double random_number_transition = gen();
    double current_rate = external_connection_rate;

    if((current_rate/total_rate)>=(random_number_transition/mt19937::max())){
        //external transition chosen
        //We need to find two conglomerates: the first will be family 0, the second family 1

        double random_number_conglomerate_zero = gen();
        double current_conglomerate = 0;
        int chosen_conglomerate_zero = -1;

        for(int cong_zero = 0; cong_zero < conglomerates.size(); cong_zero ++){
            current_conglomerate = current_conglomerate + external_sites[0][cong_zero];

            if((current_conglomerate/total_external_sites[0])>=(random_number_conglomerate_zero/mt19937::max())){
                chosen_conglomerate_zero = cong_zero;
            }
        }

        //Remove chosen cong zero from total list
        total_external_sites[1] = total_external_sites[1] - external_sites[1][chosen_conglomerate_zero];
        double random_number_conglomerate_one = gen();
        current_conglomerate = 0;
        int chosen_conglomerate_one = -1;

        for(int cong_one = 0; cong_one < conglomerates.size(); cong_one ++){
            if(cong_one != chosen_conglomerate_zero) {
                current_conglomerate = current_conglomerate + external_sites[1][cong_one];

                if ((current_conglomerate / total_external_sites[1]) >=
                    (random_number_conglomerate_one / mt19937::max())) {
                    chosen_conglomerate_one = cong_one;
                }
            }
        }

        //We need to find which sites we are going to
        double random_number_site_zero = gen();
        int chosen_site_zero = -1;

        for(int site_zero = 1; site_zero <= external_sites[0][chosen_conglomerate_zero]; site_zero ++){

            if((site_zero/external_sites[0][chosen_conglomerate_zero])>=(random_number_site_zero/mt19937::max())){
                chosen_site_zero = site_zero;
            }
        }

        //We need to find which sites we are going to
        double random_number_site_one = gen();
        int chosen_site_one = -1;

        for(int site_one = 1; site_one <= external_sites[1][chosen_conglomerate_one]; site_one ++){

            if((site_one/external_sites[1][chosen_conglomerate_one])>=(random_number_site_one/mt19937::max())){
                chosen_site_one = site_one;
            }
        }

        //TODO we have two conglomerates and their sites
        //We need to join them
        //Conglomerate zero -> add connection(Conglomerate one, connection)
        //Add the connection and all other connections to conglomerate zero
        //Update conglomerate zero
        //Remove conglomerate one
        //update rates cong zero
        //Remove all rates associated with cong one

    } else {
        int chosen_transition =-1;
        for(int transition=0; transition<=transition_rates.size(); transition++){
            current_rate = current_rate + transition_rates[transition];
            if((current_rate/total_rate)>=(random_number_transition/mt19937::max())){
                chosen_transition = transition;
            }
        }

        //Choose conglomerate
        int chosen_conglomerate =-1;
        double random_number_conglomerate = gen();
        int current_conglomerate = 0;

        for(int conglomerate=0; conglomerate<conglomerate_rates[chosen_transition].size(); conglomerate++){
            current_conglomerate = current_conglomerate + conglomerate_rates[chosen_transition][conglomerate];
            if((current_conglomerate/total_connections[chosen_transition])>=(random_number_conglomerate/mt19937::max())){
                chosen_conglomerate = conglomerate;
            }
        }
        
        //Choose bond
        int chosen_bond =-1;
        double random_number_bond = gen();

        for(int current_bond = 1; current_bond<=conglomerate_rates[chosen_transition][chosen_conglomerate]; current_bond++){
            if((current_bond/conglomerate_rates[chosen_transition][chosen_conglomerate])>=(random_number_bond/mt19937::max())){
                chosen_bond = current_bond;
            }
        }
        
        if(chosen_transition == 0){
            //Head unbinding
            vector<Conglomerate *> output = conglomerates[chosen_conglomerate]->chooseHeadUnbinding(chosen_bond);
            if(!output.empty()){
                //If there is something in the output, the unbinding has split a conglomerate
                conglomerates.push_back(output[0]);
                output[0]->index = ++conglomerate_index;
            }
        } else if (chosen_transition == 1) {
            // Head binding
            conglomerates[chosen_conglomerate]->chooseHeadBinding(chosen_bond);
        } else if (chosen_transition == 2) {
            // Tail unbinding
            vector<Conglomerate *> output = conglomerates[chosen_conglomerate]->chooseTailUnbinding(chosen_bond);
            if(!output.empty()){
                //If there is something in the output, the unbinding has split a conglomerate
                conglomerates.push_back(output[0]);
                output[0]->index = ++conglomerate_index;
                //TODO will this affect other lists? yes
            }
        } else if (chosen_transition == 3) {
            // Tail binding
            conglomerates[chosen_conglomerate]->chooseTailBinding(chosen_bond);
        } else if (chosen_transition == 4) {
            // Unbind some connected neighbours
            Polymer * new_polymer = conglomerates[chosen_conglomerate]->chooseNeighboursUnbind(chosen_bond);
            new_polymer->index = ++polymer_index;
        } else if(chosen_transition == 5){
            //Bind some unconnected neighbours
            Polymer * removed_polymer = conglomerates[chosen_conglomerate]->chooseNeighboursBind(chosen_bond);
        }
        updateRates(chosen_conglomerate);
    }

    
}