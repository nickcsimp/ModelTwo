//
// Created by Nicholas Simpson on 16/02/2022.
//

#include "Tests.h"

void Tests::run(){

    if(!testConglomerateInitialissation()){
        cout << "Error: Test conglomerate initialisation" << endl;
    }
    if(!testConglomerateUpdate()){
        cout << "Error: Test conglomerate update" << endl;
    }
    if(!testConglomerateAddConnection()){
        cout << "Error: Test conglomerate add connection" << endl;
    }
    if(!testConglomerate()){
        cout << "Error: Test conglomerate" << endl;
    }

    if(!testMiddleTailUnbinding()){
        cout << "Error: Test middle tail unbinding" << endl;
    }

    if(!testEndTailUnbinding()) {
        cout << "Error: Test End Tail Unbinding()" << endl;
    }


}

bool Tests::testConglomerateInitialissation() {

    Polymer * polymer = new Polymer(1, 6, 0);
    Conglomerate * conglomerate = new Conglomerate(polymer);

    if(conglomerate->polymers.size()!=1){
        delete conglomerate;
        delete polymer;
        return false;
    }

    if(!(*conglomerate->polymers[0] == *polymer)){
        delete conglomerate;
        delete polymer;
        return false;
    }

    if(conglomerate->polymer_connections[0].size()!=6){
        delete conglomerate;
        delete polymer;
        return false;
    }

    for(int i=0; i<6; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            delete conglomerate;
            delete polymer;
            return false;
        }
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(conglomerate->available_free_sites_list[0].size()!=6){
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=0){
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(!conglomerate->head_unbinding_list.empty() || !conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        delete conglomerate;
        delete polymer;
        return false;
    }

    delete conglomerate;

    Polymer * poly = new Polymer(2, 1, 1);
    Connection * con = new Connection(polymer, 5, poly, 0);
    vector<Connection *> connections;
    connections.push_back(con);

    conglomerate = new Conglomerate(connections);

    if(conglomerate->polymers.size()!=2){
        cout << 1 << endl;
        delete conglomerate;
        delete poly;
        delete polymer;
        return false;
    }

    if(!(*conglomerate->polymers[0] == *polymer) || !(*conglomerate->polymers[1] == *poly)){
        cout << 2 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    if(conglomerate->polymer_connections[0].size()!=6 || conglomerate->polymer_connections[1].size()!=1){
        cout << 3 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    for(int i=0; i<5; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            cout << 4 << endl;
            delete conglomerate;
            delete polymer;
            delete poly;
            return false;
        }
    }
    if(conglomerate->polymer_connections[0][5].size() != 1 || conglomerate->polymer_connections[1][0].size() != 1){
        cout << 5 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        cout << 6 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->available_free_sites_list[1].empty()){
        cout << 8 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    return true;
}


bool Tests::testConglomerateUpdate() {

    Polymer * polymer = new Polymer(1, 6, 0);
    Conglomerate * conglomerate = new Conglomerate(polymer);
    conglomerate->updateConglomerate();

    if(conglomerate->polymers.size()!=1){
        delete conglomerate;
        delete polymer;
        return false;
    }

    if(!(*conglomerate->polymers[0] == *polymer)){
        delete conglomerate;
        delete polymer;
        return false;
    }

    if(conglomerate->polymer_connections[0].size()!=6){
        delete conglomerate;
        delete polymer;
        return false;
    }

    for(int i=0; i<6; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            delete conglomerate;
            delete polymer;
            return false;
        }
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(conglomerate->available_free_sites_list[0].size()!=6){
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=0){
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(!conglomerate->head_unbinding_list.empty() || !conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        delete conglomerate;
        delete polymer;
        return false;
    }

    delete conglomerate;

    Polymer * poly = new Polymer(2, 1, 1);
    Connection * con = new Connection(polymer, 5, poly, 0);
    vector<Connection *> connections;
    connections.push_back(con);

    conglomerate = new Conglomerate(connections);

    conglomerate->updateConglomerate();

    if(conglomerate->polymers.size()!=2){
        cout << 1 << endl;
        delete conglomerate;
        delete poly;
        delete polymer;
        return false;
    }

    if(!(*conglomerate->polymers[0] == *polymer) || !(*conglomerate->polymers[1] == *poly)){
        cout << 2 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    if(conglomerate->polymer_connections[0].size()!=6 || conglomerate->polymer_connections[1].size()!=1){
        cout << 3 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    for(int i=0; i<5; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            cout << 4 << endl;
            delete conglomerate;
            delete polymer;
            delete poly;
            return false;
        }
    }
    if(conglomerate->polymer_connections[0][5].size() != 1 || conglomerate->polymer_connections[1][0].size() != 1){
        cout << 5 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        cout << 6 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->available_free_sites_list[1].empty()){
        cout << 8 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    return true;
}

bool Tests::testConglomerateAddConnection() {
    Polymer * polymer = new Polymer(1, 6,0);
    Conglomerate * conglomerate = new Conglomerate(polymer);
    Polymer * poly = new Polymer(2, 1, 1);
    Conglomerate * cong = new Conglomerate(poly);
    Connection * con = new Connection(polymer, 5, poly, 0);
    conglomerate->addConnection(cong, con);

    if(conglomerate->polymers.size()!=2){
        cout << 1 << endl;
        delete conglomerate;
        delete poly;
        delete polymer;
        return false;
    }

    if(!(*conglomerate->polymers[0] == *polymer) || !(*conglomerate->polymers[1] == *poly)){
        cout << 2 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    if(conglomerate->polymer_connections[0].size()!=6 || conglomerate->polymer_connections[1].size()!=1){
        cout << 3 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    for(int i=0; i<5; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            cout << 4 << endl;
            delete conglomerate;
            delete polymer;
            delete poly;
            return false;
        }
    }
    if(conglomerate->polymer_connections[0][5].size() != 1 || conglomerate->polymer_connections[1][0].size() != 1){
        cout << 5 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        cout << 6 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->available_free_sites_list[1].empty()){
        cout << 8 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    delete conglomerate;
    delete polymer;
    delete poly;
    return true;

}

bool Tests::testConglomerate() {
    Polymer* p_one = new Polymer(1, 6, 0);
    Polymer* p_two = new Polymer(2, 1, 1);
    Polymer* p_three = new Polymer(3, 1, 1);

    Connection * con_one = new Connection(p_one, 4, p_two, 0);
    Connection * con_two = new Connection(p_one, 5, p_three, 0);

    vector<Connection *> connections;
    connections.push_back(con_one);
    connections.push_back(con_two);

    Conglomerate * conglomerate = new Conglomerate(connections);

    bool failed = false;

    if(conglomerate->polymers.size()!=3){
        cout << 1 << endl;
        failed = true;
    }

    if(!(*conglomerate->polymers[0] == *p_one) || !(*conglomerate->polymers[1] == *p_two) ||!(*conglomerate->polymers[2] == *p_three)){
        cout << 2 << endl;
        failed = true;
    }

    if(conglomerate->polymer_connections[0].size()!=6 || conglomerate->polymer_connections[1].size()!=1 || conglomerate->polymer_connections[2].size()!=1){
        cout << 3 << endl;
        failed = true;
    }

    for(int i=0; i<4; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            cout << 4 << endl;
            failed = true;
        }
    }
    if(conglomerate->polymer_connections[0][4].size() != 1 ||conglomerate->polymer_connections[0][5].size() != 1 || conglomerate->polymer_connections[1][0].size() != 1 || conglomerate->polymer_connections[2][0].size() != 1){
        cout << 5 << endl;
        failed = true;
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        cout << 6 << endl;
        failed = true;
    }
    if(conglomerate->available_free_sites_list[0].size()!=4){
        cout << 7 << endl;
        failed = true;
    }
    if(!conglomerate->available_free_sites_list[1].empty()){
        cout << 8 << endl;
        failed = true;
    }
    if(conglomerate->head_unbinding_list.size()!=2){
        cout << 9 << endl;
        failed = true;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        failed = true;
    }
    if(!conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        failed = true;
    }
    if(conglomerate->unconnected_neighbours_list.size()!=1){
        cout << 12 << endl;
        failed = true;
    }

    delete conglomerate;
    delete p_one;
    delete p_two;
    delete p_three;
    delete con_one;
    delete con_two;
    return !failed;
}

bool Tests::testMiddleTailUnbinding(){
    Polymer* p_one = new Polymer(1, 6, 0);
    Polymer* p_two = new Polymer(2, 6, 1);

    Connection * con_one = new Connection(p_one, 0, p_two, 5);
    Connection * con_two = new Connection(p_one, 1, p_two, 4);
    Connection * con_three = new Connection(p_one, 2, p_two, 3);
    Connection * con_four = new Connection(p_one, 3, p_two, 2);
    Connection * con_five = new Connection(p_one, 4, p_two, 1);
    Connection * con_six = new Connection(p_one, 5, p_two, 0);

    vector<Connection *> connections;
    connections.push_back(con_one);
    connections.push_back(con_two);
    connections.push_back(con_three);
    connections.push_back(con_four);
    connections.push_back(con_five);
    connections.push_back(con_six);

    Conglomerate * conglomerate = new Conglomerate(connections);

    bool failed = false;

    if(conglomerate->polymers.size()!=2){
        cout << 1 << endl;
        failed = true;
    }

    if(!(*conglomerate->polymers[0] == *p_one) || !(*conglomerate->polymers[1] == *p_two) ){
        cout << 2 << endl;
        failed = true;
    }

    if(conglomerate->polymer_connections[0].size()!=6 || conglomerate->polymer_connections[1].size()!=6){
        cout << 3 << endl;
        failed = true;
    }

    for(int i=0; i<6; i++){
        if(conglomerate->polymer_connections[0][i].empty() || conglomerate->polymer_connections[1][i].empty()){
            cout << 4 << endl;
            failed = true;
        }
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        cout << 6 << endl;
        failed = true;
    }
    if(!conglomerate->available_free_sites_list[0].empty()){
        cout << 7 << endl;
        failed = true;
    }
    if(!conglomerate->available_free_sites_list[1].empty()){
        cout << 8 << endl;
        failed = true;
    }
    if(conglomerate->head_unbinding_list.size()!=2){
        cout << 9 << endl;
        failed = true;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        cout << "Head binding size " << conglomerate->head_binding_list.size() << endl;
        cout << "Tail binding size " << conglomerate->tail_binding_list.size() << endl;
        cout << "Tail unbinding size " << conglomerate->tail_unbinding_list.size() << endl;
        failed = true;
    }
    if(conglomerate->connected_neighbours_list.size()!=10){
        cout << 11 << endl;
        failed = true;
    }
    if(!conglomerate->unconnected_neighbours_list.empty()){
        cout << 12 << endl;
        failed = true;
    }

    delete conglomerate;
    delete p_one;
    delete p_two;
    delete con_one;
    delete con_two;
    delete con_three;
    delete con_four;
    delete con_five;
    delete con_six;
    return !failed;
}


bool Tests::testEndTailUnbinding(){
    Polymer* p_one = new Polymer(1, 6, 0);
    Polymer* p_two = new Polymer(2, 5, 1);

    Connection * con_two = new Connection(p_one, 1, p_two, 4);
    Connection * con_three = new Connection(p_one, 2, p_two, 3);
    Connection * con_four = new Connection(p_one, 3, p_two, 2);
    Connection * con_five = new Connection(p_one, 4, p_two, 1);
    Connection * con_six = new Connection(p_one, 5, p_two, 0);

    vector<Connection *> connections;
    connections.push_back(con_two);
    connections.push_back(con_three);
    connections.push_back(con_four);
    connections.push_back(con_five);
    connections.push_back(con_six);

    Conglomerate * conglomerate = new Conglomerate(connections);

    bool failed = false;

    if(conglomerate->polymers.size()!=2){
        cout << 1 << endl;
        failed = true;
    }

    if(!(*conglomerate->polymers[0] == *p_one) || !(*conglomerate->polymers[1] == *p_two) ){
        cout << 2 << endl;
        failed = true;
    }

    if(conglomerate->polymer_connections[0].size()!=6 || conglomerate->polymer_connections[1].size()!=5){
        cout << 3 << endl;
        failed = true;
    }

    for(int i=1; i<6; i++){
        if(conglomerate->polymer_connections[0][i].empty()){
            cout << 4 << endl;
            failed = true;
        }
    }
    for(int i=0; i<5; i++){
        if(conglomerate->polymer_connections[1][i].empty()){
            cout << 4.2 << endl;
            failed = true;
        }
    }
    if(!conglomerate->polymer_connections[0][0].empty()){
        cout << 4 << endl;
        failed = true;
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        cout << 6 << endl;
        failed = true;
    }
    if(conglomerate->available_free_sites_list[0].size() != 1){
        cout << 7 << endl;
        failed = true;
    }
    if(!conglomerate->available_free_sites_list[1].empty()){
        cout << 8 << endl;
        failed = true;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        failed = true;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty()){
        cout << 10 << endl;
        failed = true;
    }
    if(conglomerate->tail_unbinding_list.size()!=1){
        cout << 10.5 << endl;
        failed = true;
    }
    if(conglomerate->connected_neighbours_list.size()!=8){
        cout << 11 << endl;
        failed = true;
    }
    if(!conglomerate->unconnected_neighbours_list.empty()){
        cout << 12 << endl;
        failed = true;
    }

    delete conglomerate;
    delete p_one;
    delete p_two;
    delete con_two;
    delete con_three;
    delete con_four;
    delete con_five;
    delete con_six;
    return !failed;
}