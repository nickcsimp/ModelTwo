//
// Created by Nicholas Simpson on 03/02/2022.
//

#include "Conglomerate.h"

#include <utility>

Conglomerate::Conglomerate(vector<Connection *> con){
    connections = move(con);
    updateConglomerate();
}

Conglomerate::Conglomerate(Polymer * pol){
    polymers.push_back(pol);
    updateConglomerate();
    //cout << "Available free sites: " << available_free_sites_list[0].size() << ' ' << available_free_sites_list[1].size() << ' ' << available_free_sites_list[2].size() << ' ' << available_free_sites_list[3].size() << endl;
    //cout << "End unbinding: " << end_unbinding_list.size() << endl;
    //cout << "Head unbinding: " << head_unbinding_list.size() << endl;
    //cout << "Tail unbinding: " << tail_unbinding_list.size() << endl;
    //cout << "End binding: " << end_binding_list.size() << endl;
    //cout << "Tail binding: " << tail_binding_list.size() << endl;
    //cout << "Head binding: " << head_binding_list.size() << endl;
    //cout << "Polymerisation: " << unconnected_neighbours_list.size() << endl;
    //cout << "Depolymerisation: " << connected_neighbours_list.size() << endl;
}

//Calls all update functions
void Conglomerate::updateConglomerate() {
    updatePolymersInConglomerate();
    updateAvailableFreeSites();
    updateUnbindingLists();
    updateBindingLists();
    updateNeighboursBindingList();
    updateNeighboursUnbindingList();
}

//Finds all polymers used in the connections of the conglomerate
//If there are no connections - leave the polymer that is in there already
//Then find the connections on the polymer
void Conglomerate::updatePolymersInConglomerate(){
    if(connections.size()==0){
        //Check for errors
        if(polymers.size()!=1){
            cout << "ERROR: no connections in conglomerate but polymer count does not equal one." << endl;
            cout << "Conglomerate class. Update polymers in conglomerate";
        }
        //Initialise polymer_connections
        polymer_connections.clear();

        vector<Connection *> vecto; //Empty vector
        //Create a vector of empty vectors of the same size as the polymer
        vector<vector<Connection *>> vect(polymers[0]->length, vecto);
        polymer_connections.push_back(vect);
    } else {
        polymers.clear(); //Clear the lists
        polymer_connections.clear();
        for(auto & con : connections){ //Loop all connections
            for(int i=0; i<2; i++){ //Loop the polymers in connection
                bool in_list = false; //Used to determine if the polymer is already in the list
                for(int j=0; j<polymers.size(); j++){ //Loop all polymers in the list
                    if(*polymers[j] == *con->polymers_in_connection[i]){ //If the polymer is already in the list:
                        in_list=true; //Update the boolean
                        polymer_connections[j][con->indexes[i]].push_back(con); //Add connection to polymer connection list
                    }
                }
                if(!in_list){ //If the polymer is not in the list but it is in a connection:
                    polymers.push_back(con->polymers_in_connection[i]); //Add polymer to the list
                    vector<Connection *> vecto; //Empty vector
                    //Create a vector of empty vectors of the same size as the polymer
                    vector<vector<Connection *>> vect(con->polymers_in_connection[i]->length, vecto);
                    //Add connection to the correct monomer
                    vect[con->indexes[i]].push_back(con);
                    //Length 1, index 1
                    polymer_connections.push_back(vect);
                }
            }
        }
    }


}

//Loop all polymers in the conglomerate
//Available free sites are created where there is no connection found
void Conglomerate::updateAvailableFreeSites(){
    //Clear and delete free sites in list already
    for(int i=0; i<available_free_sites_list.size(); i++){
        for(int j = 0; j<available_free_sites_list[i].size(); j++){
            delete available_free_sites_list[i][j];
        }
    }
    available_free_sites_list.clear(); //Clear vector
    available_free_sites_list.emplace_back(); //Create a vector for family 0
    available_free_sites_list.emplace_back(); //Create a vector for family 1
    if(set_end_monomer){
        available_free_sites_list.emplace_back(); //Create a vector for family 0 ends
        available_free_sites_list.emplace_back(); //Create a vector for family 1 ends
    }
    for(int i=0; i<polymer_connections.size(); i++){ //Loop polymers
        int family = polymers[i]->family;
        for(int j=0; j<polymer_connections[i].size(); j++){ //Loop monomers
            if(polymer_connections[i][j].size()==0){ //If there is not a connection, we create a free site
                //If there are end monomers and this monomer is an end then we add to family+2 which holds the end monomers
                if(set_parallel_growth){
                    if(set_end_monomer && polymers[i]->has_end_monomer && j==0){
                        available_free_sites_list[family+2].push_back(new FreeSite(polymers[i], j));
                    } else {
                        available_free_sites_list[family].push_back(new FreeSite(polymers[i], j));
                    }
                } else {
                    if(set_end_monomer && polymers[i]->has_end_monomer && j==polymers[i]->length-1){
                        available_free_sites_list[family+2].push_back(new FreeSite(polymers[i], j));
                    } else {
                        available_free_sites_list[family].push_back(new FreeSite(polymers[i], j));
                    }
                }
            }
        }
    }
}

//Finds all head connections
//Finds all tail connections that can unzip
void Conglomerate::updateUnbindingLists(){

    tail_unbinding_list.clear();
    head_unbinding_list.clear();
    end_unbinding_list.clear();

    //Loop Connections, find valid unbinding
    for(auto & con : connections){

        //Is the connection on the end of a template?
        bool end_monomer = false;
        if(set_weakened_template_end) {
            if (set_parallel_growth) {
                if (((con->polymers_in_connection[0]->is_template && con->indexes[0] == 0) || (con->polymers_in_connection[1]->is_template && con->indexes[1] == 0))) {
                    end_monomer = true;
                }
            } else {
                if (((con->polymers_in_connection[0]->is_template && con->indexes[0] == con->polymers_in_connection[0]->length -1) || (con->polymers_in_connection[1]->is_template && con->indexes[1] == con->polymers_in_connection[1]->length -1))) {
                    end_monomer = true;
                }
            }
        }
        //Is the connection part of an end monomer connection?
        if(set_end_monomer) {
            if (set_parallel_growth) {
                if ((con->polymers_in_connection[0]->has_end_monomer && con->indexes[0] == 0) || (con->polymers_in_connection[1]->has_end_monomer && con->indexes[1] == 0)) {
                    end_monomer = true;
                }
            } else {
                if ((con->polymers_in_connection[0]->has_end_monomer && con->indexes[0] == con->polymers_in_connection[0]->length -1) || (con->polymers_in_connection[1]->has_end_monomer && con->indexes[1] == con->polymers_in_connection[1]->length -1)) {
                    end_monomer = true;
                }
            }
        }

        //There is a weakened end on the template polymer with index 0 on the length - 1 monomer
        if(end_monomer){
            end_unbinding_list.push_back(con);
        } else if(con->indexes[0] == 0 || con->indexes[1] == 0){ //If either are 0, then the connection is a head
            head_unbinding_list.push_back(con); //Add connection to list
        } else if(con->indexes[0] == con->polymers_in_connection[0]->length-1 || con->indexes[1] == con->polymers_in_connection[1]->length-1){
            //End of tail, definitely tail unbinding
            tail_unbinding_list.push_back(con); //Add connection to list

        } else {
            //Need to check if it is in the middle of a zip or able to be zipped

            //Find the polymer in the conglomerate saving list
            //Only need to look at one as they need to be both bound for there to be a middle monomer
            int polymer_one = -1;
            for(int i=0; i<polymers.size(); i++){
                if(*polymers[i]==*con->polymers_in_connection[0]){
                    polymer_one = i;
                    break;
                }
            }

            bool middle_connection = true;
            //If the connection is surrounded by connections to the same polymers, then no unzipping can occur
            for(int i=-1; i<2; i+=2){ //Looking at the previous(-1) and the next (+1)
                if(polymer_connections[polymer_one][con->indexes[0]+i].size()==0){ //If there is no connection,
                    middle_connection = false; //The tail can unbind
                    break;
                }

                //If the first neighbour is not the same as either polymer then it's not central
                if(!(*polymer_connections[polymer_one][con->indexes[0]+i][0]->polymers_in_connection[0] == *con->polymers_in_connection[0] || *polymer_connections[polymer_one][con->indexes[0]+i][0]->polymers_in_connection[0] == *con->polymers_in_connection[1])){
                    middle_connection = false;
                    break;
                }
                //If the second neighbour is not the same as either polymer then it's not central
                if(!(*polymer_connections[polymer_one][con->indexes[0]+i][0]->polymers_in_connection[1] == *con->polymers_in_connection[0] || *polymer_connections[polymer_one][con->indexes[0]+i][0]->polymers_in_connection[1] == *con->polymers_in_connection[1])){
                    middle_connection = false;
                    break;
                }
            }
            if(!middle_connection){
                tail_unbinding_list.push_back(con);
            }
        }
    }
}


//Loop all monomers on all polymers
//If the monomer has no connection AND a neighbouring monomer has a connection then we need to check it
void Conglomerate::updateBindingLists(){
    for(int i=0; i<head_binding_list.size(); i++){
        delete head_binding_list[i];
    }
    for(int i=0; i<tail_binding_list.size(); i++){
        delete tail_binding_list[i];
    }
    for(int i=0; i<end_binding_list.size(); i++){
        delete end_binding_list[i];
    }
    head_binding_list.clear();
    tail_binding_list.clear();
    end_binding_list.clear();
    for(int pol=0; pol<polymer_connections.size(); pol++) { //Loop polymers
        for (int mon = 0; mon < polymer_connections[pol].size(); mon++) { //Loop monomers
            if (polymer_connections[pol][mon].empty()) { //No connection on this site
                //Check the possibility of zipping on both sides
                //We can check toward the head if mon is not 0
                //We can check toward the tail if mon is not length-1
                if (mon != polymers[pol]->length -1) {
                    checkBindingValidity(pol, mon, 1);
                }
                if (mon != 0) {
                    checkBindingValidity(pol, mon, -1);
                }
            }
        }
    }
}

//Checks if there is a possible zipping binding at position pol mon in the direction indicated
//Direction +1 will check lagging edge binding
//Direction -1 will check leading edge binding
void Conglomerate::checkBindingValidity(int polymer, int monomer, int direction){
    //If the vector is empty there is no connection and therefore there can be no zipping
    if (!polymer_connections[polymer][monomer + direction].empty()) {
        //Need to know which polymer is which in the connection
        int connected_polymer;
        if (*polymer_connections[polymer][monomer + direction][0]->polymers_in_connection[0] == *polymers[polymer]) {
            connected_polymer = 1;
        } else {
            connected_polymer = 0;
        }


        //Check if there is a tail on the connected polymer
        bool is_tail = false;

        if((direction == 1 && set_parallel_growth == true) || (direction == -1 && set_parallel_growth == false)){
            //Need to check there is a leading edge (we are not at the beginning of the polymer)
            if(polymer_connections[polymer][monomer + direction][0]->indexes[connected_polymer] != 0){
                //Binding possible
                is_tail = true;
            }
        } else {
            //Need to check there is a lagging edge (we are not at the end of the polymer)
            if(polymer_connections[polymer][monomer + direction][0]->indexes[connected_polymer] != polymer_connections[polymer][monomer + direction][0]->polymers_in_connection[connected_polymer]->length-1){
                //Binding possible
                is_tail = true;
            }
        }

        if(is_tail){
            //Check that there is no connection on the neighbouring monomer
            //First need to find the polymer in polymer_connections
            int con_pol_ind = -1;
            for(int i=0; i<polymers.size(); i++){
                if(*polymers[i] == *polymer_connections[polymer][monomer + direction][0]->polymers_in_connection[connected_polymer]){
                    con_pol_ind = i;
                    break;
                }
            }

            //If there is no connection on this index then there can be binding
            int growth = 1;
            if(set_parallel_growth){
                growth = -1;
            }
            if(polymer_connections[con_pol_ind][polymer_connections[polymer][monomer + direction][0]->indexes[connected_polymer]+growth*direction].empty()){
                //See if either point is an end monomer
                bool p1_end = false;
                bool p2_end = false;
                if(set_end_monomer && polymers[polymer]->has_end_monomer){
                    if(set_parallel_growth && monomer == 0){
                        p1_end = true;
                    } else if(!set_parallel_growth && monomer == polymers[polymer]->length -1){
                        p1_end = true;
                    }
                }
                if(set_end_monomer && polymers[con_pol_ind]->has_end_monomer){
                    if(set_parallel_growth && polymer_connections[polymer][monomer + direction][0]->indexes[connected_polymer]+growth*direction == 0){
                        p2_end = true;
                    } else if(!set_parallel_growth && polymer_connections[polymer][monomer + direction][0]->indexes[connected_polymer]+growth*direction == polymers[con_pol_ind]->length -1){
                        p2_end = true;
                    }
                }

                if(!(p1_end ^ p2_end)) { //Exclusive or so that if one is an end and the other is not then we can't bind
                    //Binding possible!!
                    Connection *connection_found = new Connection(polymers[polymer], monomer, polymers[con_pol_ind],
                                                                  polymer_connections[polymer][monomer +
                                                                                               direction][0]->indexes[connected_polymer] +
                                                                  growth * direction);
                    if(p1_end){ //if p1_end is true then p2_end is true and this is an end connection
                        //End connection
                        bool already_exists = false;
                        for (auto &con: end_binding_list) {
                            if (*connection_found == *con) {
                                already_exists = true;
                            }
                        }
                        if (!already_exists) {
                            //If it has not been added, then please add
                            end_binding_list.push_back(connection_found);
                        }
                    } else if (monomer == 0 ||
                               polymer_connections[polymer][monomer + direction][0]->indexes[connected_polymer] +
                               growth * direction == 0) { //Check if there is a head involved
                        //Check this connection has not already been found
                        bool already_exists = false;
                        for (auto &con: head_binding_list) {
                            if (*connection_found == *con) {
                                already_exists = true;
                            }
                        }
                        if (!already_exists) {
                            //If it has not been added, then please add
                            head_binding_list.push_back(connection_found);
                        }
                    } else {
                        //Check this connection has not already been found
                        bool already_exists = false;
                        for (auto &con: tail_binding_list) {
                            if (*connection_found == *con) {
                                already_exists = true;
                            }
                        }
                        if (!already_exists) {
                            //If it has not been added, then please add
                            tail_binding_list.push_back(connection_found);
                        }
                    }
                }
            }

        }

    }
}


//Loop polymers, check if head is neighbouring a tail
void Conglomerate::updateNeighboursBindingList(){
    for(int i=0; i<unconnected_neighbours_list.size(); i++){
        delete unconnected_neighbours_list[i];
    }
    unconnected_neighbours_list.clear();
    for(int pol = 0; pol<polymer_connections.size(); pol++) { //Loop polymers
        //If we have parallel growth, end monomers, and this polymer has an end monomer, then we cant grow in this direction
        if (!(set_end_monomer && polymers[pol]->has_end_monomer && set_parallel_growth)) {
            if (!polymer_connections[pol][0].empty()) { //if the first monomer (head) has a connection, then we need to check the neighbour
                //Find polymer connected to *the* polymer
                int connected_polymer;
                if (*polymer_connections[pol][0][0]->polymers_in_connection[0] == *polymers[pol]) {
                    connected_polymer = 1;
                } else {
                    connected_polymer = 0;
                }

                bool neighbour_possible = false;
                if (set_parallel_growth && polymer_connections[pol][0][0]->indexes[connected_polymer] != 0) {
                    //If the connection is on the first monomer of the connected polymer and the growth is parallel then no neighbour is possible
                    neighbour_possible = true;
                } else if (!set_parallel_growth && (polymer_connections[pol][0][0]->indexes[connected_polymer] !=
                                                    (polymer_connections[pol][0][0]->polymers_in_connection[connected_polymer]->length -
                                                     1))) {
                    //If the connection is on the final monomer of the connected polymer and the growth is antiparallel then no neighbour is possible
                    neighbour_possible = true;
                }

                if (neighbour_possible) {
                    //Find the polymer in polymers list
                    int con_pol_ind = -1;
                    for (int i = 0; i < polymers.size(); i++) {
                        if (*polymers[i] ==
                            *polymer_connections[pol][0][0]->polymers_in_connection[connected_polymer]) {
                            con_pol_ind = i;
                        }
                    }

                    int growth = 1;
                    if (set_parallel_growth) {
                        growth = -1;
                    }
                    //If there is not a connection on the neighbouring monomer of the connected polymer, there will be no neighbours binding
                    if (!polymer_connections[con_pol_ind][polymer_connections[pol][0][0]->indexes[connected_polymer] +
                                                          growth].empty()) {
                        //Find polymer connected to the connected polymer - the neighbour
                        int neighbour;
                        if (*polymer_connections[con_pol_ind][
                                polymer_connections[pol][0][0]->indexes[connected_polymer] +
                                growth][0]->polymers_in_connection[0] == *polymers[con_pol_ind]) {
                            neighbour = 1;
                        } else {
                            neighbour = 0;
                        }
                        Polymer *neighbour_polymer = polymer_connections[con_pol_ind][
                                polymer_connections[pol][0][0]->indexes[connected_polymer] +
                                growth][0]->polymers_in_connection[neighbour];
                        //Now if the neighbour is connected with its final monomer, then we have neighbours!!!!
                        if (polymer_connections[con_pol_ind][
                                    polymer_connections[pol][0][0]->indexes[connected_polymer] +
                                    growth][0]->indexes[neighbour] == neighbour_polymer->length - 1) {
                            //If there is antiparallel growth, end monomers, and the neighbour polymer has an end monomer, then we cant grow in this direction
                            if (!(set_end_monomer && neighbour_polymer->has_end_monomer && !set_parallel_growth)) {
                                unconnected_neighbours_list.push_back(
                                        new UnconnectedNeighbours(neighbour_polymer, polymers[pol]));
                            }
                        }
                    }
                }

            }
        }
    }
}

//Loop polymers, Loop monomers, if two consecutive monomers have connections we check if they are connected to the same polymer

void Conglomerate::updateNeighboursUnbindingList(){
    for(int i=0; i<connected_neighbours_list.size(); i++){
        delete connected_neighbours_list[i];
    }
    connected_neighbours_list.clear();
    for(int pol=0; pol<polymer_connections.size(); pol++){ //Loop polymers
        if(!(set_template_indestructible && polymers[pol]->is_template)) { //If template is indestructible we don't allow the template to break
            for (int mon = 0; mon < polymer_connections[pol].size() -
                                    1; mon++) { //Loop monomers, but not the last as we are checking mon+1
                //We need mon and mon+1 to both have a connection for them to potentially be unbinding neighbours
                if (!polymer_connections[pol][mon].empty() && !polymer_connections[pol][mon + 1].empty()) {
                    //Find which polymers they are connected to
                    int polymer_connected_to_mon = -1;
                    if (*polymer_connections[pol][mon][0]->polymers_in_connection[0] == *polymers[pol]) {
                        polymer_connected_to_mon = 1;
                    } else {
                        polymer_connected_to_mon = 0;
                    }

                    int polymer_connected_to_mon_plus_one = -1;
                    if (*polymer_connections[pol][mon + 1][0]->polymers_in_connection[0] == *polymers[pol]) {
                        polymer_connected_to_mon_plus_one = 1;
                    } else {
                        polymer_connected_to_mon_plus_one = 0;
                    }

                    //If the connected polymers are the same then we have unbindable neighbours
                    if (*polymer_connections[pol][mon +
                                                  1][0]->polymers_in_connection[polymer_connected_to_mon_plus_one] ==
                        *polymer_connections[pol][mon][0]->polymers_in_connection[polymer_connected_to_mon]) {
                        //They can bind!!
                        connected_neighbours_list.push_back(new ConnectedNeighbours(polymers[pol], mon));
                    }
                }
            }
        }
    }
}

void Conglomerate::addConnection(Conglomerate * conglomerate_being_joined, Connection * new_connection){
    //Add the new connection
    connections.push_back(new_connection);
    //Add all connections from the other conglomerate
    for(auto & con : conglomerate_being_joined->connections){
        connections.push_back(con);
    }
    //update this conglomerate
    updateConglomerate();
}

//Remove chosen bond from connection list
//Update conglomerate
vector<Conglomerate *> Conglomerate::chooseHeadUnbinding(int chosen_bond){
    vector<Conglomerate *> output;
    //Find which connection is to be removed
    int connection_index = -1;
    for(int i=0; i<connections.size(); i++){
        if(*connections[i]==*head_unbinding_list[chosen_bond]){
            //This is the connection to be removed
            connection_index = i;
            break;
        }
    }
    //Remove connection
    Connection * delete_me = connections[connection_index];
    connections.erase(connections.begin()+connection_index);

    //Create new conglomerate if needed
    output = checkSeparation(head_unbinding_list[chosen_bond]);

    delete delete_me;
    //Update conglomerate
    updateConglomerate();

    return output;
}

//Remove chosen bond from connection list
//Update conglomerate
vector<Conglomerate *> Conglomerate::chooseEndUnbinding(int chosen_bond){
    vector<Conglomerate *> output;
    //Find which connection is to be removed
    int connection_index = -1;
    for(int i=0; i<connections.size(); i++){
        if(*connections[i]==*end_unbinding_list[chosen_bond]){
            //This is the connection to be removed
            connection_index = i;
            break;
        }
    }
    //Remove connection
    Connection * delete_me = connections[connection_index];
    connections.erase(connections.begin()+connection_index);

    //Create new conglomerate if needed
    output = checkSeparation(end_unbinding_list[chosen_bond]);

    delete delete_me;
    //Update conglomerate
    updateConglomerate();
    return output;
}

//Remove chosen bond from connection list
//Update conglomerate
vector<Conglomerate *> Conglomerate::chooseTailUnbinding(int chosen_bond){
    vector<Conglomerate *> output;
    //Find which connection is to be removed
    int connection_index = -1;
    for(int i=0; i<connections.size(); i++){
        if(*connections[i]==*tail_unbinding_list[chosen_bond]){
            //This is the connection to be removed
            connection_index = i;
            break;
        }
    }
    //Remove connection
    Connection * delete_me = connections[connection_index];
    connections.erase(connections.begin()+connection_index);

    output = checkSeparation(tail_unbinding_list[chosen_bond]);

    delete delete_me;
    //Update conglomerate
    updateConglomerate();

    return output;
}

//Check for conglomerate separation
//If there are no connections remaining, we need to figure this out
//Loop connections, if the two polymers are still connected in some way then dont worry
//If they aren't then the conglomerate hath split
//We return a vector of new conglomerates
vector<Conglomerate *> Conglomerate::checkSeparation(Connection * removed_connection){
    vector<Conglomerate *> output;

    //If there are no connections, the two polymers have split and are now individuals
    //Create a new conglomerate made of one polymer, remove the polymer from this conglomerate
    if(connections.empty()){
        int removed_poly = -1;
        for(int i=0; i<polymers.size(); i++){
            if(*removed_connection->polymers_in_connection[1] == *polymers[i]){
                removed_poly = i;
            }
        }
        polymers.erase(polymers.begin()+removed_poly);
        output.push_back(new Conglomerate(removed_connection->polymers_in_connection[1]));
        return output;
    }

    //If there are still connections, we need to see if the conglom has separated

    bool separated = true;

    for(auto & con : connections){ //Loop all connections
        //If the connection has both polymers, then the conglomerate has not separated
        if((*con->polymers_in_connection[0] == *removed_connection->polymers_in_connection[0]) &&
           (*con->polymers_in_connection[1] == *removed_connection->polymers_in_connection[1])){
            separated = false;
        } else if((*con->polymers_in_connection[1] == *removed_connection->polymers_in_connection[0]) &&
                  (*con->polymers_in_connection[0] == *removed_connection->polymers_in_connection[1])){
            separated = false;
        }
        //If no separation, we leave output empty and no new conglomerate will be added to the system
    }

    if(separated){
        //Using the polymer_connection list will make life easier so need to update it
        int p_one_index = -1;
        int p_two_index = -1;
        for(int pol = 0; pol<polymers.size(); pol++){ //Loop polymers
            if(*removed_connection->polymers_in_connection[0] == *polymers[pol]){
                //Find the 1st polymer in the connection
                p_one_index = pol;
            }
            if(*removed_connection->polymers_in_connection[1] == *polymers[pol]){
                //Find the second polymer in the connection
                p_two_index = pol;
            }
        }

        //Remove the connections from the polymer_connections list
        polymer_connections[p_one_index][removed_connection->indexes[0]].clear();
        polymer_connections[p_two_index][removed_connection->indexes[1]].clear();

        //Find the tree for one of the polymers
        vector<Polymer *> empty_vector;
        vector<Polymer *> connected_polymers = getTree(polymers[p_one_index], empty_vector);

        //If the vector is just one polymer long, there will be no connections
        //In this case, we make the new cong from a polymer rather than a list of connections
        if(connected_polymers.size()==1){
            //Make new conglomerate and add to output
            output.push_back(new Conglomerate(connected_polymers[0]));
            //Remove the polymer from this conglomerate
            //We need to do this so that the update polymers in conglomerate is not confused
            int removed_poly = -1;
            for(int i=0; i<polymers.size(); i++){
                if(*removed_connection->polymers_in_connection[1] == *polymers[i]){
                    removed_poly = i;
                }
            }
            polymers.erase(polymers.begin()+removed_poly);
        } else {
            //We can use the connections here
            //Loop connections, find polymers included in connected polymers
            //Remove the connection from this polymer but add to the new polymer
            vector<int> removed_connections;
            vector<Connection *> new_cong_conns;
            for(int i=0; i<connections.size(); i++){//Loop all connections
                for(auto & pol : connected_polymers){ //Loop all connected polymers
                    //If they equal, we do our jiggerypuff
                    //We only check one polymer because that should be enough
                    if(*pol == *connections[i]->polymers_in_connection[0]){
                        removed_connections.push_back(i);
                        new_cong_conns.push_back(connections[i]);
                        break;
                    }
                }
            }

            //Need to erase the connections found
            //Do in reverse so the pointers don't point to the wrong thing
            for(int i=removed_connections.size()-1; i>-1; i--){
                connections.erase(connections.begin()+removed_connections[i]);
            }

            //Need to erase polymers in connected polymers
            vector<int> polymers_removal_indexes;
            for(int j=0; j<polymers.size(); j++){
                for(auto & pol : connected_polymers){
                    if(*polymers[j] == *pol){
                        polymers_removal_indexes.push_back(j);
                    }
                }
            }

            //Need to erase the polymers found
            //Do in reverse so the pointers don't point to the wrong thing
            for(int i=polymers_removal_indexes.size()-1; i>-1; i--){
                polymers.erase(polymers.begin()+polymers_removal_indexes[i]);
            }

            //Need to add the connections to the new conglomerate
            output.push_back(new Conglomerate(new_cong_conns));
        }
    }
    return output;
}

//Takes a polymer and finds which other polymers are connected to it but are not in the list already provided
vector<Polymer *> Conglomerate::getTree(Polymer * p, vector<Polymer *> connected_polymers){
    //First we add our polymer to the connected polymer list so that it is not investigated in the future
    connected_polymers.push_back(p);
    //Find polymer in polymers list
    int ind = -1;
    for(int pol = 0; pol<polymers.size(); pol++){ //Loop polymers
        if(*p == *polymers[pol]){
            ind = pol;
        }
    }
    //Find the tree for any connected polymer found - ensuring up-to-date connected polymers is input

    //Loop monomers on polymer
    for(int i=0; i<polymer_connections[ind].size(); i++){
        //If the vector is not empty then we need to investigate the connection
        if(!polymer_connections[ind][i].empty()){
            //Find which polymer is new
            int new_polymer = -1;
            if(*polymer_connections[ind][i][0]->polymers_in_connection[0] == *p){
                new_polymer = 1;
            } else {
                new_polymer = 0;
            }
            //Need to see if polymer has already been investigated
            bool in_list = false;
            for(auto & poly : connected_polymers){ //Loop investigated polymers
                if(*poly == *polymer_connections[ind][i][0]->polymers_in_connection[new_polymer]){
                    //If it is found then it has already been investigated
                    in_list = true;
                }
            }
            //If the polymer is not in the list
            //We want to get the polymers that it is connected to
            if(!in_list){
                vector<Polymer *> new_tree = getTree(polymer_connections[ind][i][0]->polymers_in_connection[new_polymer], connected_polymers);
                for(int i=0; i<new_tree.size(); i++){
                    bool in_list = false;
                    for(int j= 0 ; j<connected_polymers.size(); j++){
                        if(*new_tree[i]==*connected_polymers[j]){
                            in_list = true;
                        }
                    }
                    if(!in_list){
                        connected_polymers.push_back(new_tree[i]);
                    }
                }
            }
        }
    }
    return connected_polymers;
}

//Add chosen bond to connection list
//Update conglomerate
void Conglomerate::chooseHeadBinding(int chosen_bond){

    //Add connection
    connections.push_back(new Connection(head_binding_list[chosen_bond]->polymers_in_connection[0], head_binding_list[chosen_bond]->indexes[0], head_binding_list[chosen_bond]->polymers_in_connection[1], head_binding_list[chosen_bond]->indexes[1]));
    //Update conglomerate
    updateConglomerate();
}

//Add chosen bond to connection list
//Update conglomerate
void Conglomerate::chooseTailBinding(int chosen_bond){
    //Add connection
    connections.push_back(new Connection(tail_binding_list[chosen_bond]->polymers_in_connection[0], tail_binding_list[chosen_bond]->indexes[0], tail_binding_list[chosen_bond]->polymers_in_connection[1], tail_binding_list[chosen_bond]->indexes[1]));

    //Update conglomerate
    updateConglomerate();
}

//Add chosen bond to connection list
//Update conglomerate
void Conglomerate::chooseEndBinding(int chosen_bond){
    //Add connection
    connections.push_back(new Connection(end_binding_list[chosen_bond]->polymers_in_connection[0], end_binding_list[chosen_bond]->indexes[0], end_binding_list[chosen_bond]->polymers_in_connection[1], end_binding_list[chosen_bond]->indexes[1]));
    //Update conglomerate
    updateConglomerate();
}


//Join polymers
//Update connections
//Update conglomerate
//Return removed polymer to remove from the system
Polymer * Conglomerate::chooseNeighboursBind(int chosen_bond){
    //Get polymers for easy access
    Polymer * polymer_being_extended = unconnected_neighbours_list[chosen_bond]->polymer_one;
    Polymer * polymer_to_be_removed = unconnected_neighbours_list[chosen_bond]->polymer_two;

    //Extend polymer being extended
    int original_length = polymer_being_extended->length;
    polymer_being_extended->length = polymer_being_extended->length + polymer_to_be_removed->length;

    //Find polymer_to_be_removed in conglomerate
    int ptbr_index = -1;
    for(int i=0; i<polymers.size(); i++){
        if(*polymers[i] == *polymer_to_be_removed){
            ptbr_index = i;
            break;
        }
    }

    //Transfer connections
    for(int i=0; i<polymer_connections[ptbr_index].size(); i++){ //Loop monomers on polymer being removed
        for(auto & con : polymer_connections[ptbr_index][i]){ //For every connection on the monomer (really not needed)
            for(int j=0; j<2; j++){ //Loop polymers in connection
                if(*con->polymers_in_connection[j] == *polymer_to_be_removed){ //If the polymer is the polymer being removed
                    //Any connection involving removed polymer needs to be changed to the extended polymer
                    con->polymers_in_connection[j] = polymer_being_extended; //Change the polymer in the connection
                    con->indexes[j] = con->indexes[j] + original_length;
                }
            }
        }
    }
    if(!set_parallel_growth && set_end_monomer && polymer_to_be_removed->has_end_monomer){
        polymer_being_extended->has_end_monomer = true;
    }

    //Remove polymer
    polymers.erase(polymers.begin()+ptbr_index);
    //Update conglomerate
    updateConglomerate();

    return polymer_to_be_removed;
}

//Return new polymer to add to the system
Polymer * Conglomerate::chooseNeighboursUnbind(int chosen_bond){
    //Get polymer for easy access
    Polymer * polymer_being_separated = connected_neighbours_list[chosen_bond]->polymer;

    //We make a new polymer
    //It will be the same family as polymer being separated
    //Index starts from 0 and is the smaller of the two in the interaction
    //Therefore the con_neighbour switching a polymer length 8 into a 2 and 6 will have index 1
    //This polymer has length 6=8-1-1
    Polymer * new_polymer;
    if(set_end_monomer && !set_parallel_growth && polymer_being_separated->has_end_monomer){
        //If antiparallel growth, the spearated monomer will lose its end and the new polymer will gain it
        new_polymer = new Polymer(-1, polymer_being_separated->length - connected_neighbours_list[chosen_bond]->index -1, polymer_being_separated->family, true);
        polymer_being_separated->has_end_monomer = false;
    } else {
        //With parallel growth (or no set end monomer or pbs does not have a head) the new polymer will have no end
        new_polymer = new Polymer(-1, polymer_being_separated->length - connected_neighbours_list[chosen_bond]->index -1, polymer_being_separated->family, false);
    }

    //Shorten polymer being separated
    //The polymer being separated keeps its head, so it has a length index+1 (1+1=2)
    polymer_being_separated->length = connected_neighbours_list[chosen_bond]->index + 1;


    //Find polymer_being_separated in conglomerate
    int pbs_index = -1;
    for(int i=0; i<polymers.size(); i++){
        if(*polymers[i] == *polymer_being_separated){
            pbs_index = i;
            break;
        }
    }


    //Transfer connections
    for(int i=0; i<polymer_connections[pbs_index].size(); i++){
        for(auto & con : polymer_connections[pbs_index][i]){
            for(int j=0; j<2; j++){
                if(*con->polymers_in_connection[j] == *polymer_being_separated && con->indexes[j]>connected_neighbours_list[chosen_bond]->index){
                    //Any connection involving separated polymer above the index of separation needs to be transferred
                    con->polymers_in_connection[j] = new_polymer;
                    con->indexes[j] = con->indexes[j] - connected_neighbours_list[chosen_bond]->index - 1;
                }
            }
        }
    }


    //Add to conglomerate
    polymers.push_back(new_polymer);
    //Update conglomerate
    updateConglomerate();
    return new_polymer;
}

void Conglomerate::deleteConglomerate(){
    for(int i=tail_binding_list.size()-1; i>0; i--){
        delete tail_binding_list[i];
    }
    for(int i=head_binding_list.size()-1; i>0; i--){
        delete head_binding_list[i];
    }
    for(int i=connected_neighbours_list.size()-1; i>0; i--){
        delete connected_neighbours_list[i];
    }
    for(int i=unconnected_neighbours_list.size()-1; i>0; i--){
        delete unconnected_neighbours_list[i];
    }
    for(int i=available_free_sites_list.size()-1; i>0; i--){
        for(int j=available_free_sites_list[i].size()-1; j>0; j--){
            delete available_free_sites_list[i][j];
        }
    }
}
