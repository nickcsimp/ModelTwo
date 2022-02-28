//
// Created by Nicholas Simpson on 07/02/2022.
//

#include "System.h"


System::System(settings set){

    sett = set;

    //Initialise index count
    conglomerate_index = -1;
    polymer_index = -1;

    //Give vectors correct sizes
    vector<int> empty_vector;
    vector<vector<int>> vect(6, empty_vector);
    conglomerate_rates = vect;

    vector<int> v(6, 0);
    total_connections = v;
    vector<double> ve(6, 0);
    transition_rates = ve;

    vector<int> vec(2, 0);
    total_external_sites = vec;

    vector<vector<int>> vecto(2, empty_vector);
    external_sites = vecto;


    //Initialise system parameters
    k = set.k;
    k0 = set.k0;
    G_bb = set.G_bb;
    G_spec = set.G_spec;
    G_gen = set.G_gen;
    M_eff = set.M_eff;

    //Make template polymer into a conglomerate
    Polymer * template_polymer = new Polymer(++polymer_index, set.template_length, 0);
    Conglomerate * template_conglomerate = new Conglomerate(template_polymer, set);
    template_conglomerate->index = ++conglomerate_index;
    addConglomerate(template_conglomerate);

    //Add to lengths list
    vector<int> victor(set.template_length, 0);
    victor[set.template_length-1]++;
    lengths = victor;

    //Add all free monomers to the system as polymers within conglomerates
    for(int i=0; i<set.monomers_family_zero; i++){
        Polymer * new_poly = new Polymer(++polymer_index, 1, 0);
        Conglomerate * new_cong = new Conglomerate(new_poly, set);
        new_cong->index = ++conglomerate_index;
        addConglomerate(new_cong);
        lengths[0]++;
    }
    for(int i=0; i<set.monomers_family_one; i++){
        Polymer * new_poly = new Polymer(++polymer_index, 1, 1);
        Conglomerate * new_cong = new Conglomerate(new_poly, set);
        new_cong->index = ++conglomerate_index;
        addConglomerate(new_cong);
        lengths[0]++;
    }

    template_indestructible = set.template_indestructible;
    monomer_count_is_constant = set.monomer_count_is_constant;
    no_rebinding = set.no_rebinding;

}


void System::updateRates(int cong){
    total_rate = 0;

    //Remove ext rate from chosen conglomerate
    for(int i=0; i<conglomerates.size(); i++) { //Loop conglomerates
        if (i != cong) { //Conglomerate can't bind within itself here
            external_connection_rate = external_connection_rate - external_sites[0][i] * external_sites[1][cong] * k0;
            external_connection_rate = external_connection_rate - external_sites[1][i] * external_sites[0][cong] * k0;
        }
    }
    //Update external rates
    total_external_sites[0] = total_external_sites[0] - external_sites[0][cong];
    total_external_sites[1] = total_external_sites[1] - external_sites[1][cong];
    external_sites[0][cong] = conglomerates[cong]->available_free_sites_list[0].size();
    external_sites[1][cong] = conglomerates[cong]->available_free_sites_list[1].size();
    total_external_sites[0] = total_external_sites[0] + external_sites[0][cong];
    total_external_sites[1] = total_external_sites[1] + external_sites[1][cong];



    //Add new to total
    for(int i=0; i<conglomerates.size(); i++) { //Loop conglomerates
        if (i != cong) { //Conglomerate can't bind within itself here
            external_connection_rate = external_connection_rate + external_sites[0][i] * external_sites[1][cong] * k0;
            external_connection_rate = external_connection_rate + external_sites[1][i] * external_sites[0][cong] * k0;
        }
    }
    total_rate = total_rate + external_connection_rate;
    transition_rates[0] = transition_rates[0] - conglomerate_rates[0][cong]*(k0*exp(G_spec+G_gen));
    total_connections[0] = total_connections[0] - conglomerate_rates[0][cong];
    conglomerate_rates[0][cong] = conglomerates[cong]->head_unbinding_list.size();
    transition_rates[0] = transition_rates[0] + conglomerate_rates[0][cong]*(k0*exp(G_spec+G_gen));
    total_connections[0] = total_connections[0] + conglomerate_rates[0][cong];
    total_rate = total_rate + transition_rates[0];

    transition_rates[1] = transition_rates[1] - conglomerate_rates[1][cong]*(k0*M_eff);
    total_connections[1] = total_connections[1] - conglomerate_rates[1][cong];
    conglomerate_rates[1][cong] = conglomerates[cong]->head_binding_list.size();
    transition_rates[1] = transition_rates[1] + conglomerate_rates[1][cong]*(k0*M_eff);
    total_connections[1] = total_connections[1] + conglomerate_rates[1][cong];
    total_rate = total_rate + transition_rates[1];

    transition_rates[2] = transition_rates[2] - conglomerate_rates[2][cong]*(k0*exp(G_spec));
    total_connections[2] = total_connections[2] - conglomerate_rates[2][cong];
    conglomerate_rates[2][cong] = conglomerates[cong]->tail_unbinding_list.size();
    transition_rates[2] = transition_rates[2] + conglomerate_rates[2][cong]*(k0*exp(G_spec));
    total_connections[2] = total_connections[2] + conglomerate_rates[2][cong];

    total_rate = total_rate + transition_rates[2];

    transition_rates[3] = transition_rates[3] - conglomerate_rates[3][cong]*(k0*M_eff);
    total_connections[3] = total_connections[3] - conglomerate_rates[3][cong];
    conglomerate_rates[3][cong] = conglomerates[cong]->tail_binding_list.size();
    transition_rates[3] = transition_rates[3] + conglomerate_rates[3][cong]*(k0*M_eff);
    total_connections[3] = total_connections[3] + conglomerate_rates[3][cong];
    total_rate = total_rate + transition_rates[3];

    transition_rates[4] = transition_rates[4] - conglomerate_rates[4][cong]*(k*exp(G_bb-G_gen));
    total_connections[4] = total_connections[4] - conglomerate_rates[4][cong];
    conglomerate_rates[4][cong] = conglomerates[cong]->connected_neighbours_list.size();
    transition_rates[4] = transition_rates[4] + conglomerate_rates[4][cong]*(k*exp(G_bb-G_gen));
    total_connections[4] = total_connections[4] + conglomerate_rates[4][cong];
    total_rate = total_rate + transition_rates[4];

    transition_rates[5] = transition_rates[5] - conglomerate_rates[5][cong]*k;
    total_connections[5] = total_connections[5] - conglomerate_rates[5][cong];
    conglomerate_rates[5][cong] = conglomerates[cong]->unconnected_neighbours_list.size();
    transition_rates[5] = transition_rates[5] + conglomerate_rates[5][cong]*k;
    total_connections[5] = total_connections[5] + conglomerate_rates[5][cong];
    total_rate = total_rate + transition_rates[5];
}


bool System::chooseTransition(double seed){
    mt19937 gen(seed);

    //First choose a transition
    //External rates is weird, so we see if that is chosen first

    double random_number_transition = gen();
    double current_rate = external_connection_rate;

    int chosen_conglomerate_zero = -1;
    int chosen_conglomerate_one = -1;
    int chosen_site_zero = -1;
    int chosen_site_one = -1;
    int first = 0;
    int second = 1;

    if(total_rate == 0){
        cout << "No transitions possible" << endl;
        return false;
    }

    if((current_rate/total_rate)>=(random_number_transition/mt19937::max())){
        //external transition chosen
        //We need to find two conglomerates: the first will be family 0, the second family 1

        //If either family has just one free site, we need to make sure that it is picked
        //Family 0 will be chosen automatically first so just need to check family 1

        int family_one_has_sites_in_this_many_conglomerates = 0;
        for(int i=0; i<external_sites[1].size(); i++){
            if(external_sites[1][i]!=0){
                family_one_has_sites_in_this_many_conglomerates++;
            }
        }

        if(family_one_has_sites_in_this_many_conglomerates==1){
            //We need the chosen_conglomerate_one to be the conglomerate with the one external site
            first = 1;
            second = 0;
        }

        double random_number_conglomerate_zero = gen();
        double current_conglomerate = 0;

        for(int cong_zero = 0; cong_zero < conglomerates.size(); cong_zero ++){
            current_conglomerate = current_conglomerate + external_sites[first][cong_zero];

            if((current_conglomerate/total_external_sites[first])>=(random_number_conglomerate_zero/mt19937::max())){
                chosen_conglomerate_zero = cong_zero;
                break;
            }
        }

        //Remove chosen cong zero from total list
        total_external_sites[second] = total_external_sites[second] - external_sites[second][chosen_conglomerate_zero];
        double random_number_conglomerate_one = gen();
        current_conglomerate = 0;

        for(int cong_one = 0; cong_one < conglomerates.size(); cong_one ++){
            if(cong_one != chosen_conglomerate_zero) {
                current_conglomerate = current_conglomerate + external_sites[second][cong_one];

                if ((current_conglomerate / total_external_sites[second]) >=
                    (random_number_conglomerate_one / mt19937::max())) {
                    chosen_conglomerate_one = cong_one;
                    break;
                }
            }
        }
        total_external_sites[second] = total_external_sites[second] + external_sites[second][chosen_conglomerate_zero];

        //We need to find which sites we are going to use
        double random_number_site_zero = gen();

        for(int site_zero = 1; site_zero <= external_sites[first][chosen_conglomerate_zero]; site_zero ++){
            if((site_zero/external_sites[first][chosen_conglomerate_zero])>=(random_number_site_zero/mt19937::max())){
                chosen_site_zero = site_zero - 1;
                break;
            }
        }

        //We need to find which sites we are going to use
        double random_number_site_one = gen();

        for(int site_one = 1; site_one <= external_sites[second][chosen_conglomerate_one]; site_one ++){
            if((site_one/external_sites[second][chosen_conglomerate_one])>=(random_number_site_one/mt19937::max())){
                chosen_site_one = site_one - 1;
                break;
            }
        }

        //Make new connection between the two sites
        FreeSite * free_site_zero = conglomerates[chosen_conglomerate_zero]->available_free_sites_list[first][chosen_site_zero];
        FreeSite * free_site_one = conglomerates[chosen_conglomerate_one]->available_free_sites_list[second][chosen_site_one];

        Connection * new_connection = new Connection(free_site_zero->polymer, free_site_zero->index, free_site_one->polymer, free_site_one->index);

        //Add the connection to conglomerate zero and transfer all connections from one to zero
        conglomerates[chosen_conglomerate_zero]->addConnection(conglomerates[chosen_conglomerate_one], new_connection);
        //Update rates of conglomerate zero
        updateRates(chosen_conglomerate_zero);
        //Once the rates have been updated, we remove everything related to conglomerate one and then the system should run as normal
        //Remove conglomerate one from the system
        removeConglomerate(chosen_conglomerate_one);
        //Now all lists have one less thing and the total rates should have been updated
    } else {
        int chosen_transition =-1;
        for(int transition=0; transition<=transition_rates.size(); transition++){
            current_rate = current_rate + transition_rates[transition];
            if((current_rate/total_rate)>=(random_number_transition/mt19937::max())){
                chosen_transition = transition;
                break;
            }
        }

        //Choose conglomerate
        int chosen_conglomerate =-1;
        double random_number_conglomerate = gen();
        double current_conglomerate = 0;
        for(int conglomerate=0; conglomerate<conglomerate_rates[chosen_transition].size(); conglomerate++){
            current_conglomerate = current_conglomerate + conglomerate_rates[chosen_transition][conglomerate];
            if((current_conglomerate/total_connections[chosen_transition])>=(random_number_conglomerate/mt19937::max())){
                chosen_conglomerate = conglomerate;
                break;
            }
        }

        //Choose bond
        int chosen_bond =-1;
        double random_number_bond = gen();
        for(double current_bond = 1; current_bond<=conglomerate_rates[chosen_transition][chosen_conglomerate]; current_bond++){
            if((current_bond/conglomerate_rates[chosen_transition][chosen_conglomerate])>=(random_number_bond/mt19937::max())){
                chosen_bond = current_bond - 1;
                break;
            }
        }

        if(chosen_transition == 0){
            //Head unbinding
            vector<Conglomerate *> output = conglomerates[chosen_conglomerate]->chooseHeadUnbinding(chosen_bond);
            if(!output.empty()){
                //If there is something in the output, the unbinding has split a conglomerate
                if(no_rebinding){
                    int monomers_to_add;
                    //if there is no rebinding, we need to find which conglomerate has the template
                    //we delete the one which hasn't
                    output[0]->index = ++conglomerate_index;
                    bool has_template = false;
                    for(auto & pol: output[0]->polymers){
                        if(pol->family == 0){
                            has_template = true;
                        }
                    }

                    if(has_template){
                        //add new
                        addConglomerate(output[0]);
                        monomers_to_add = conglomerates[chosen_conglomerate]->polymers[0]->length;
                        //delete old
                        removeConglomerate(chosen_conglomerate);
                    } else {
                        monomers_to_add = output[0]->polymers[0]->length;
                        delete output[0];
                    }
                    if(monomer_count_is_constant) {
                        for (int i = 0; i < monomers_to_add; i++) {
                            //Create new monomers (polymers and conglomerates)
                            Polymer * p = new Polymer(++polymer_index, 1, 1);
                            Conglomerate * c = new Conglomerate(p, sett);
                            c->index = ++conglomerate_index;
                            addConglomerate(c);
                            lengths[0]++;
                        }
                    }
                } else {
                    //If there is rebinding we can just allow things to continue
                    output[0]->index = ++conglomerate_index;
                    addConglomerate(output[0]);
                }
            }
        } else if (chosen_transition == 1) {
            //Head Binding
            conglomerates[chosen_conglomerate]->chooseHeadBinding(chosen_bond);
        } else if (chosen_transition == 2) {
            // Tail unbinding

            vector<Conglomerate *> output = conglomerates[chosen_conglomerate]->chooseTailUnbinding(chosen_bond);
            if(!output.empty()){
                //If there is something in the output, the unbinding has split a conglomerate
                if(no_rebinding){
                    int monomers_to_add;
                    //if there is no rebinding, we need to find which conglomerate has the template
                    //we delete the one which hasn't
                    output[0]->index = ++conglomerate_index;
                    bool has_template = false;
                    for(auto & pol: output[0]->polymers){
                        if(pol->family == 0){
                            has_template = true;
                        }
                    }

                    if(has_template){
                        //add new
                        addConglomerate(output[0]);
                        monomers_to_add = conglomerates[chosen_conglomerate]->polymers[0]->length;
                        //delete old
                        removeConglomerate(chosen_conglomerate);
                    } else {
                        monomers_to_add = output[0]->polymers[0]->length;
                        delete output[0];
                    }
                    if(monomer_count_is_constant) {
                        for (int i = 0; i < monomers_to_add; i++) {
                            //Create new monomers (polymers and conglomerates)
                            Polymer * p = new Polymer(++polymer_index, 1, 1);
                            Conglomerate * c = new Conglomerate(p, sett);
                            c->index = ++conglomerate_index;
                            addConglomerate(c);
                            lengths[0]++;
                        }
                    }
                } else {
                    //If there is rebinding we can just allow things to continue
                    output[0]->index = ++conglomerate_index;
                    addConglomerate(output[0]);
                }
            }
        } else if (chosen_transition == 3) {
            // Tail binding
            conglomerates[chosen_conglomerate]->chooseTailBinding(chosen_bond);
        } else if (chosen_transition == 4) {
            // Unbind some connected neighbours
            //Remove lengths of polymer from histogram
            Polymer * old_polymer = conglomerates[chosen_conglomerate]->connected_neighbours_list[chosen_bond]->polymer;
            lengths[old_polymer->length - 1]--;

            //Do transition and get new polymer
            Polymer * new_polymer = conglomerates[chosen_conglomerate]->chooseNeighboursUnbind(chosen_bond);
            new_polymer->index = ++polymer_index;

            //Add lengths to histogram
            lengths[old_polymer->length-1]++;
            lengths[new_polymer->length-1]++;

        } else if(chosen_transition == 5){
            //Bind some unconnected neighbours

            //Get polymers involved to remove from length histogram and join again
            Polymer * p_one = conglomerates[chosen_conglomerate]->unconnected_neighbours_list[chosen_bond]->polymer_one;
            Polymer * p_two = conglomerates[chosen_conglomerate]->unconnected_neighbours_list[chosen_bond]->polymer_two;
            lengths[p_one->length-1]--;
            lengths[p_two->length-1]--;
            if(lengths.size()<p_one->length+p_two->length){
                for(int i=0; i<p_one->length+p_two->length-lengths.size(); i++){
                    lengths.push_back(0);
                }
            }
            lengths[p_one->length+p_two->length-1]++;

            Polymer * removed_polymer = conglomerates[chosen_conglomerate]->chooseNeighboursBind(chosen_bond);
            delete removed_polymer;
        }
        updateRates(chosen_conglomerate);
    }
    return true;
}

void System::removeConglomerate(int cong){
    total_rate = 0;
    for(int i=0; i<conglomerates.size(); i++) { //Loop conglomerates
        if (i != cong) { //Conglomerate can't bind within itself here
            external_connection_rate = external_connection_rate - external_sites[0][i] * external_sites[1][cong] * k0;
            external_connection_rate = external_connection_rate - external_sites[1][i] * external_sites[0][cong] * k0;
        }
    }
    delete conglomerates[cong];
    conglomerates.erase(conglomerates.begin()+cong);

    total_rate = total_rate + external_connection_rate;

    total_external_sites[0] = total_external_sites[0] - external_sites[0][cong];
    total_external_sites[1] = total_external_sites[1] - external_sites[1][cong];
    external_sites[0].erase(external_sites[0].begin()+cong);
    external_sites[1].erase(external_sites[1].begin()+cong);


    transition_rates[0] = transition_rates[0] - conglomerate_rates[0][cong]*(k0*exp(G_spec+G_gen));
    total_connections[0] = total_connections[0] - conglomerate_rates[0][cong];
    conglomerate_rates[0].erase(conglomerate_rates[0].begin()+cong);
    total_rate = total_rate + transition_rates[0];

    transition_rates[1] = transition_rates[1] - conglomerate_rates[1][cong]*(k0*M_eff);
    total_connections[1] = total_connections[1] - conglomerate_rates[1][cong];
    conglomerate_rates[1].erase(conglomerate_rates[1].begin()+cong);
    total_rate = total_rate + transition_rates[1];

    transition_rates[2] = transition_rates[2] - conglomerate_rates[2][cong]*(k0*exp(G_spec));
    total_connections[2] = total_connections[2] - conglomerate_rates[2][cong];
    conglomerate_rates[2].erase(conglomerate_rates[2].begin()+cong);
    total_rate = total_rate + transition_rates[2];

    transition_rates[3] = transition_rates[3] - conglomerate_rates[3][cong]*(k0*M_eff);
    total_connections[3] = total_connections[3] - conglomerate_rates[3][cong];
    conglomerate_rates[3].erase(conglomerate_rates[3].begin()+cong);
    total_rate = total_rate + transition_rates[3];

    transition_rates[4] = transition_rates[4] - conglomerate_rates[4][cong]*(k*exp(G_bb-G_gen));
    total_connections[4] = total_connections[4] - conglomerate_rates[4][cong];
    conglomerate_rates[4].erase(conglomerate_rates[4].begin()+cong);
    total_rate = total_rate + transition_rates[4];

    transition_rates[5] = transition_rates[5] - conglomerate_rates[5][cong]*k;
    total_connections[5] = total_connections[5] - conglomerate_rates[5][cong];
    conglomerate_rates[5].erase(conglomerate_rates[5].begin()+cong);
    total_rate = total_rate + transition_rates[5];
}

void System::addConglomerate(Conglomerate * new_cong){

    conglomerates.push_back(new_cong);

    int cong = conglomerates.size()-1;

    total_rate = 0;

    external_sites[0].push_back(conglomerates[cong]->available_free_sites_list[0].size());
    external_sites[1].push_back(conglomerates[cong]->available_free_sites_list[1].size());
    total_external_sites[0] = total_external_sites[0] + external_sites[0][cong];
    total_external_sites[1] = total_external_sites[1] + external_sites[1][cong];

    for(int i=0; i<conglomerates.size(); i++) { //Loop conglomerates
        if (i != cong) { //Conglomerate can't bind within itself here
            external_connection_rate = external_connection_rate + external_sites[0][i] * external_sites[1][cong] * k0;
            external_connection_rate = external_connection_rate + external_sites[1][i] * external_sites[0][cong] * k0;
        }
    }
    total_rate = total_rate + external_connection_rate;

    conglomerate_rates[0].push_back(conglomerates[cong]->head_unbinding_list.size());
    transition_rates[0] = transition_rates[0] + conglomerate_rates[0][cong]*(k0*exp(G_spec+G_gen));
    total_connections[0] = total_connections[0] + conglomerate_rates[0][cong];
    total_rate = total_rate + transition_rates[0];

    conglomerate_rates[1].push_back(conglomerates[cong]->head_binding_list.size());
    transition_rates[1] = transition_rates[1] + conglomerate_rates[1][cong]*(k0*M_eff);
    total_connections[1] = total_connections[1] + conglomerate_rates[1][cong];
    total_rate = total_rate + transition_rates[1];

    conglomerate_rates[2].push_back(conglomerates[cong]->tail_unbinding_list.size());
    transition_rates[2] = transition_rates[2] + conglomerate_rates[2][cong]*(k0*exp(G_spec));
    total_connections[2] = total_connections[2] + conglomerate_rates[2][cong];
    total_rate = total_rate + transition_rates[2];

    conglomerate_rates[3].push_back(conglomerates[cong]->tail_binding_list.size());
    transition_rates[3] = transition_rates[3] + conglomerate_rates[3][cong]*(k0*M_eff);
    total_connections[3] = total_connections[3] + conglomerate_rates[3][cong];
    total_rate = total_rate + transition_rates[3];

    conglomerate_rates[4].push_back(conglomerates[cong]->connected_neighbours_list.size());
    transition_rates[4] = transition_rates[4] + conglomerate_rates[4][cong]*(k*exp(G_bb-G_gen));
    total_connections[4] = total_connections[4] + conglomerate_rates[4][cong];
    total_rate = total_rate + transition_rates[4];

    conglomerate_rates[5].push_back(conglomerates[cong]->unconnected_neighbours_list.size());
    transition_rates[5] = transition_rates[5] + conglomerate_rates[5][cong]*k;
    total_connections[5] = total_connections[5] + conglomerate_rates[5][cong];
    total_rate = total_rate + transition_rates[5];
}