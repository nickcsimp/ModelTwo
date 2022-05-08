//
// Created by Nicholas Simpson on 16/02/2022.
//

#include "Tests.h"

void Tests::run(){

    bool semfz = set_end_monomers_family_zero;
    bool semfo = set_end_monomers_family_one;
    bool sem = set_end_monomer;
    bool sti = set_template_indestructible;
    bool smcic = set_monomer_count_is_constant;
    bool snr = set_no_rebinding;
    bool swte = set_weakened_template_end;
    bool spg = set_parallel_growth;
    double ss = set_seed;
    double sk = set_k;
    double sk0 = set_k0;
    double sGbb = set_G_bb;
    double sGspec = set_G_spec;
    double sGgen = set_G_gen;
    double sMeff = set_M_eff;
    double sGend = set_G_end;
    int smfz = set_monomers_family_zero;
    int smfo = set_monomers_family_one;
    int stl = set_template_length;
    int strl = set_transition_limit;
    int spl = set_polymer_limit;
    double stil = set_time_limit;
    double sv = set_volume;

    bool srt = set_run_tests;
    bool smah = set_make_animated_histogram;
    bool smfh = set_make_final_histogram;
    bool smalg = set_make_average_length_graph;
    bool smldp = set_make_length_distribution_plots;
    bool smi = set_make_images;

    set_end_monomer = false;

    cout << "Testing Conglomerate Initialisation...." << endl;
    if(!testConglomerateInitialissation()){
        cout << "Error: Test conglomerate initialisation" << endl;
    }
    cout << "Testing Conglomerate Update...." << endl;
    if(!testConglomerateUpdate()){
        cout << "Error: Test conglomerate update" << endl;
    }
    cout << "Testing Conglomerate Add Connection...." << endl;
    if(!testConglomerateAddConnection()){
        cout << "Error: Test conglomerate add connection" << endl;
    }
    cout << "Testing Conglomerate...." << endl;
    if(!testConglomerate()){
        cout << "Error: Test conglomerate" << endl;
    }
    cout << "Testing Middle Tail Unbinding...." << endl;
    if(!testMiddleTailUnbinding()){
        cout << "Error: Test middle tail unbinding" << endl;
    }
    cout << "Testing End Tail Unbinding...." << endl;
    if(!testEndTailUnbinding()) {
        cout << "Error: Test End Tail Unbinding()" << endl;
    }
    cout << "Testing Weakened End Bond...." << endl;
    if(!testWeakenedEndBond()) {
        cout << "Error: Test Weakened End Bond" << endl;
    }
    cout << "Testing Two and Three Polymer Conglomerate...." << endl;
    if(!testTwoAndThreePolymerCong()) {
        cout << "Error: Test Two and Three Polymer Conglomerate()" << endl;
    }
    cout << "Testing Polymerisation Equilibrium...." << endl;
    if(!testPolymerisationEquilibrium()) {
        cout << "Error: testPolymerisationEquilibrium()" << endl;
    }
    cout << "Testing Dimerisation Equilibrium...." << endl;
    if(!testDimerisationEquilibrium()) {
        cout << "Error: testDimerisationEquilibrium()" << endl;
    }
    cout << "Testing Tail Binding Equilibrium...." << endl;
    if(!testTailBindEquilibrium()) {
        cout << "Error: testTailBindEquilibrium()" << endl;
    }

    cout << "Testing Growth Directions...." << endl;
    if(!testGrowthDirections()) {
        cout << "Error: testGrowthDirections()" << endl;
    }

    cout << "Testing Weak End With Growth Directions...." << endl;
    if(!testWeakEndWithGrowthDirections()) {
        cout << "Error: testWeakEndWithGrowthDirections()" << endl;
    }
    cout << "Testing End Monomer...." << endl;
    if(!testEndMonomer()) {
        cout << "Error: testEndMonomer()" << endl;
    }
    cout << "Testing Correct Monomer Binding...." << endl;
    if(!testCorrectMonomerBinding()) {
        cout << "Error: testCorrectMonomerBinding()" << endl;
    }


    cout << "Tests Complete" << endl;

    set_end_monomers_family_zero = semfz;
    set_end_monomers_family_one = semfo;
    set_end_monomer = sem;
    set_parallel_growth = spg;
    set_template_indestructible = sti;
    set_monomer_count_is_constant = smcic;
    set_no_rebinding = snr;
    set_seed = ss;
    set_weakened_template_end = swte;
    set_k = sk;
    set_k0 = sk0;
    set_G_bb = sGbb;
    set_G_spec = sGspec;
    set_G_gen = sGgen;
    set_G_end = sGend;
    set_M_eff = sMeff;
    set_monomers_family_zero = smfz;
    set_monomers_family_one = smfo;
    set_template_length = stl;
    set_transition_limit = strl;
    set_polymer_limit = spl;
    set_time_limit = stil;
    set_volume = sv;

    set_run_tests = srt;
    set_make_animated_histogram = smah;
    set_make_final_histogram = smfh;
    set_make_average_length_graph = smalg;
    set_make_length_distribution_plots = smldp;
    set_make_images = smi;
}

bool Tests::testConglomerateInitialissation() {

    set_weakened_template_end = false;

    Polymer * polymer = new Polymer(1, 6, 0, false);
    Conglomerate * conglomerate = new Conglomerate(polymer);

    if(conglomerate->polymers.size()!=1){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }

    if(!(*conglomerate->polymers[0] == *polymer)){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }

    if(conglomerate->polymer_connections[0].size()!=6){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }

    for(int i=0; i<6; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            conglomerate->deleteConglomerate();
            delete conglomerate;
            delete polymer;
            return false;
        }
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(conglomerate->available_free_sites_list[0].size()!=6){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=0){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(!conglomerate->head_unbinding_list.empty() || !conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }

    conglomerate->deleteConglomerate();
    delete conglomerate;

    Polymer * poly = new Polymer(2, 1, 1, false);
    Connection * con = new Connection(polymer, 5, poly, 0);
    vector<Connection *> connections;
    connections.push_back(con);

    conglomerate = new Conglomerate(connections);

    if(conglomerate->polymers.size()!=2){
        cout << 1 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        return false;
    }

    if(!(*conglomerate->polymers[0] == *polymer) || !(*conglomerate->polymers[1] == *poly)){
        cout << 2 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    if(conglomerate->polymer_connections[0].size()!=6 || conglomerate->polymer_connections[1].size()!=1){
        cout << 3 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    for(int i=0; i<5; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            cout << 4 << endl;
            conglomerate->deleteConglomerate();
            delete conglomerate;
            delete polymer;
            delete poly;
            return false;
        }
    }
    if(conglomerate->polymer_connections[0][5].size() != 1 || conglomerate->polymer_connections[1][0].size() != 1){
        cout << 5 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        cout << 6 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->available_free_sites_list[1].empty()){
        cout << 8 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    return true;
}


bool Tests::testConglomerateUpdate() {

    set_weakened_template_end = false;

    Polymer * polymer = new Polymer(1, 6, 0, false);
    Conglomerate * conglomerate = new Conglomerate(polymer);
    conglomerate->updateConglomerate();

    if(conglomerate->polymers.size()!=1){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }

    if(!(*conglomerate->polymers[0] == *polymer)){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }

    if(conglomerate->polymer_connections[0].size()!=6){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }

    for(int i=0; i<6; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            conglomerate->deleteConglomerate();
            delete conglomerate;
            delete polymer;
            return false;
        }
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(conglomerate->available_free_sites_list[0].size()!=6){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=0){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(!conglomerate->head_unbinding_list.empty() || !conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        return false;
    }

    conglomerate->deleteConglomerate();
    delete conglomerate;

    Polymer * poly = new Polymer(2, 1, 1, false);
    Connection * con = new Connection(polymer, 5, poly, 0);
    vector<Connection *> connections;
    connections.push_back(con);

    conglomerate = new Conglomerate(connections);

    conglomerate->updateConglomerate();

    if(conglomerate->polymers.size()!=2){
        cout << 1 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        return false;
    }

    if(!(*conglomerate->polymers[0] == *polymer) || !(*conglomerate->polymers[1] == *poly)){
        cout << 2 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        delete con;
        return false;
    }

    if(conglomerate->polymer_connections[0].size()!=6 || conglomerate->polymer_connections[1].size()!=1){
        cout << 3 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        delete con;
        return false;
    }

    for(int i=0; i<5; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            cout << 4 << endl;
            conglomerate->deleteConglomerate();
            delete conglomerate;
            delete polymer;
            delete poly;
            delete con;
            return false;
        }
    }
    if(conglomerate->polymer_connections[0][5].size() != 1 || conglomerate->polymer_connections[1][0].size() != 1){
        cout << 5 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        delete con;
        return false;
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        cout << 6 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        delete con;
        return false;
    }
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        delete con;
        return false;
    }
    if(!conglomerate->available_free_sites_list[1].empty()){
        cout << 8 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        delete con;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        delete con;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        delete con;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        delete con;
        return false;
    }
    conglomerate->deleteConglomerate();
    delete conglomerate;
    delete polymer;
    delete poly;
    delete con;
    return true;
}

bool Tests::testConglomerateAddConnection() {
    set_weakened_template_end = false;
    Polymer * polymer = new Polymer(1, 6,0, false);
    Conglomerate * conglomerate = new Conglomerate(polymer);
    Polymer * poly = new Polymer(2, 1, 1, false);
    Conglomerate * cong = new Conglomerate(poly);
    Connection * con = new Connection(polymer, 5, poly, 0);
    conglomerate->addConnection(cong, con);

    if(conglomerate->polymers.size()!=2){
        cout << 1 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        cong->deleteConglomerate();
        delete cong;
        return false;
    }

    if(!(*conglomerate->polymers[0] == *polymer) || !(*conglomerate->polymers[1] == *poly)){
        cout << 2 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        cong->deleteConglomerate();
        delete cong;
        return false;
    }

    if(conglomerate->polymer_connections[0].size()!=6 || conglomerate->polymer_connections[1].size()!=1){
        cout << 3 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        cong->deleteConglomerate();
        delete cong;
        return false;
    }

    for(int i=0; i<5; i++){
        if(!conglomerate->polymer_connections[0][i].empty()){
            cout << 4 << endl;
            conglomerate->deleteConglomerate();
            delete conglomerate;
            delete poly;
            delete polymer;
            delete con;
            cong->deleteConglomerate();
            delete cong;
            return false;
        }
    }
    if(conglomerate->polymer_connections[0][5].size() != 1 || conglomerate->polymer_connections[1][0].size() != 1){
        cout << 5 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        cong->deleteConglomerate();
        delete cong;
        return false;
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        cout << 6 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        cong->deleteConglomerate();
        delete cong;
        return false;
    }
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        cong->deleteConglomerate();
        delete cong;
        return false;
    }
    if(!conglomerate->available_free_sites_list[1].empty()){
        cout << 8 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        cong->deleteConglomerate();
        delete cong;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        cong->deleteConglomerate();
        delete cong;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        cong->deleteConglomerate();
        delete cong;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete poly;
        delete polymer;
        delete con;
        cong->deleteConglomerate();
        delete cong;
        return false;
    }
    conglomerate->deleteConglomerate();
    delete conglomerate;
    delete poly;
    delete polymer;
    delete con;
    cong->deleteConglomerate();
    delete cong;
    return true;

}

bool Tests::testConglomerate() {
    set_weakened_template_end = false;
    Polymer* p_one = new Polymer(1, 6, 0, false);
    Polymer* p_two = new Polymer(2, 1, 1, false);
    Polymer* p_three = new Polymer(3, 1, 1, false);

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

    conglomerate->deleteConglomerate();
    delete conglomerate;
    delete p_one;
    delete p_two;
    delete p_three;
    delete con_one;
    delete con_two;
    return !failed;
}

bool Tests::testMiddleTailUnbinding(){
    set_weakened_template_end = false;
    set_parallel_growth = true;
    Polymer* p_one = new Polymer(1, 6, 0, false);
    p_one->is_template = true;
    Polymer* p_two = new Polymer(2, 6, 1, false);

    Connection * con_one = new Connection(p_one, 0, p_two, 0);
    Connection * con_two = new Connection(p_one, 1, p_two, 1);
    Connection * con_three = new Connection(p_one, 2, p_two, 2);
    Connection * con_four = new Connection(p_one, 3, p_two, 3);
    Connection * con_five = new Connection(p_one, 4, p_two, 4);
    Connection * con_six = new Connection(p_one, 5, p_two, 5);

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
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        failed = true;
    }
    if(conglomerate->tail_unbinding_list.size()!=1){
        cout << 9.5 << endl;
        failed = true;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_binding_list.empty()){
        cout << 10 << endl;
        cout << "Head binding size " << conglomerate->head_binding_list.size() << endl;
        cout << "Tail binding size " << conglomerate->tail_binding_list.size() << endl;
        cout << "Tail unbinding size " << conglomerate->tail_unbinding_list.size() << endl;
        failed = true;
    }
    if(set_template_indestructible){
        if (conglomerate->connected_neighbours_list.size() != 5) {
            cout << 11 << endl;
            failed = true;
        }
    } else {
        if (conglomerate->connected_neighbours_list.size() != 10) {
            cout << 11 << endl;
            failed = true;
        }
    }
    if(!conglomerate->unconnected_neighbours_list.empty()){
        cout << 12 << endl;
        failed = true;
    }

    conglomerate->deleteConglomerate();
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
    set_weakened_template_end = false;
    set_parallel_growth = true;
    Polymer* p_one = new Polymer(1, 6, 0, false);
    p_one->is_template = true;
    Polymer* p_two = new Polymer(2, 5, 1, false);

    Connection * con_two = new Connection(p_one, 1, p_two, 0);
    Connection * con_three = new Connection(p_one, 2, p_two, 1);
    Connection * con_four = new Connection(p_one, 3, p_two, 2);
    Connection * con_five = new Connection(p_one, 4, p_two, 3);
    Connection * con_six = new Connection(p_one, 5, p_two, 4);

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
    if(set_template_indestructible){
        if (conglomerate->connected_neighbours_list.size() != 4) {
            cout << 11 << endl;
            failed = true;
        }
    } else {
        if (conglomerate->connected_neighbours_list.size() != 8) {
            cout << 11 << endl;
            failed = true;
        }
    }
    if(!conglomerate->unconnected_neighbours_list.empty()){
        cout << 12 << endl;
        failed = true;
    }

    conglomerate->deleteConglomerate();
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

bool Tests::testWeakenedEndBond(){
    //Repeat the head and tail experiments below but with the weakened bond active

    // 1 free monomer and a template length 1
    // Every interaction is with the `end' monomer so Gbb and Gspec and Ggen should make no difference
    // If G_end = ln(2) then we will be unconnected twice as much as connected
    vector<double> seeds = {101, 201, 301, 401, 501};
    vector<double> Z;
    vector<double> Z_connected;
    vector<double> Z_unconnected;
    int transition_limit = 1000;
    for(auto & seed : seeds) {
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = true;
        set_seed = seed;
        set_volume = 1;
        set_k = 1;
        set_k0 = 1;
        set_G_bb = -1;
        set_G_spec = 1;
        set_G_gen = -(1+log(0.5));
        set_M_eff = 100;
        set_G_end = log(2);
        set_monomers_family_zero = 0;
        set_monomers_family_one = 1;
        set_template_length = 1;
        set_transition_limit = 1000;


        mt19937 gen(seed);

        System *system = new System();

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates[0]->polymers.size() == 2) {
                //Now connected so we add time to unconnected list
                hist[1] = hist[1] + system->simulation_time - previous_time;
            } else if (system->conglomerates[0]->polymers.size() == 1) {
                //Now unconnected so we add time to connected list
                hist[0] = hist[0] + system->simulation_time - previous_time;
            }
            previous_time = system->simulation_time;
            count++;
        }

        // Total transitions will increase with the binding rate and the transition limit but rates are 1/sec so no factor
        Z.push_back((count-transition_limit*set_k0*4/3)/sqrt(transition_limit*set_k0*4/3));
        //We have the two rates the same so the time spent should be equal
        Z_connected.push_back((hist[0]-transition_limit*0.33)/sqrt(transition_limit*0.33));
        Z_unconnected.push_back((hist[1]-transition_limit*0.66)/sqrt(transition_limit*0.66));
    }

    double sum = accumulate(Z.begin(), Z.end(), 0.0);
    double mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 0.1."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 0.1." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 0.1." << endl;
    }


    //Template length 3
    //Copy length 2 attached by the head on index 1
    //Tail should flip between the end monomer
    // If G_end is ln(2*Meff) we are twice as likely to be unconnected
    seeds = {103, 203, 303, 403, 503};
    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    for(auto & seed : seeds) {
        //Initialisations
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = true;
        set_seed = seed;
        set_k = 0.1;
        set_k0 = 0.1;
        set_G_bb = -100000000;
        set_G_spec = -1000000;
        set_G_gen = -1000000;
        set_M_eff = 10;
        set_G_end = log(20);
        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_template_length = 0;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        Polymer *template_polymer = new Polymer(-1, 3, 0, false);
        template_polymer->is_template = true;
        Polymer *copy_polymer = new Polymer(-1, 3, 1, false);

        Connection * con = new Connection(template_polymer, 1, copy_polymer, 2);
        Conglomerate * cong = new Conglomerate({con});
        System *system = new System(cong);

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates.size() == 1) {
                if(system->conglomerates[0]->connections.size()==2){
                    //Now connected so we add time to unconnected list
                    hist[1] = hist[1] + system->simulation_time - previous_time;
                } else if (system->conglomerates[0]->connections.size() == 1) {
                    //Now unconnected so we add time to connected list
                    hist[0] = hist[0] + system->simulation_time - previous_time;
                }
            }
            previous_time = system->simulation_time;
            count++;

        }

        //The effective concentration makes a difference to the rate so must be included
        Z.push_back((count-transition_limit*set_k0*set_M_eff*4/3)/sqrt(transition_limit*set_k0*set_M_eff*4/3));
        //Equal times are expected for each state
        Z_connected.push_back((hist[0]-transition_limit*0.33)/sqrt(transition_limit*033));
        Z_unconnected.push_back((hist[1]-transition_limit*0.66)/sqrt(transition_limit*0.66));
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 0.2."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 0.2." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 0.2." << endl;
    }


    //Template length 3
    //Copy length 2 attached by the head on index 1
    //Tail should flip between the end monomer
    // If G_end is ln(2*Meff) we are twice as likely to be unconnected
    // halve the rate, halve the transitions
    seeds = {103, 203, 303, 403, 503};
    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    for(auto & seed : seeds) {
        //Initialisations
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = true;
        set_seed = seed;
        set_k = 0.1;
        set_k0 = 0.05;
        set_G_bb = -100000000;
        set_G_spec = -1000000;
        set_G_gen = -1000000;
        set_M_eff = 10;
        set_G_end = log(20);
        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_template_length = 0;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        Polymer *template_polymer = new Polymer(-1, 3, 0, false);
        template_polymer->is_template = true;
        Polymer *copy_polymer = new Polymer(-1, 3, 1, false);

        Connection * con = new Connection(template_polymer, 1, copy_polymer, 2);
        Conglomerate * cong = new Conglomerate({con});
        System *system = new System(cong);

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates.size() == 1) {
                if(system->conglomerates[0]->connections.size()==2){
                    //Now connected so we add time to unconnected list
                    hist[1] = hist[1] + system->simulation_time - previous_time;
                } else if (system->conglomerates[0]->connections.size() == 1) {
                    //Now unconnected so we add time to connected list
                    hist[0] = hist[0] + system->simulation_time - previous_time;
                }
            }
            previous_time = system->simulation_time;
            count++;

        }

        //The effective concentration makes a difference to the rate so must be included
        Z.push_back((count-transition_limit*set_k0*set_M_eff*4/3)/sqrt(transition_limit*set_k0*set_M_eff*4/3));
        //Equal times are expected for each state
        Z_connected.push_back((hist[0]-transition_limit*0.33)/sqrt(transition_limit*033));
        Z_unconnected.push_back((hist[1]-transition_limit*0.66)/sqrt(transition_limit*0.66));
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 0.3."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 0.3." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 0.3." << endl;
    }

    // 1 free monomer and a template length 1
    // Every interaction is with the `end' monomer so Gbb and Gspec and Ggen should make no difference
    // If G_end = ln(2) then we will be unconnected twice as much as connected
    // Now change G_end = ln(4) so we spend more time unconnected but increase Ggen loads so we want to be connected
    seeds = {101, 201, 301, 401, 501};
    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    for(auto & seed : seeds) {
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = true;
        set_seed = seed;
        set_k = 1;
        set_k0 = 1;
        set_G_bb = -1;
        set_G_spec = 1;
        set_G_gen = -100;
        set_M_eff = 100;
        set_G_end = log(4);
        set_monomers_family_zero = 0;
        set_monomers_family_one = 1;
        set_template_length = 1;
        set_transition_limit = 1000;


        mt19937 gen(seed);

        System *system = new System();

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates[0]->polymers.size() == 2) {
                //Now connected so we add time to unconnected list
                hist[1] = hist[1] + system->simulation_time - previous_time;
            } else if (system->conglomerates[0]->polymers.size() == 1) {
                //Now unconnected so we add time to connected list
                hist[0] = hist[0] + system->simulation_time - previous_time;
            }
            previous_time = system->simulation_time;
            count++;
        }

        // Total transitions will increase with the binding rate and the transition limit but rates are 1/sec so no factor
        Z.push_back((count-transition_limit*set_k0/0.625)/sqrt(transition_limit*set_k0/0.625));
        //We have the two rates the same so the time spent should be equal
        Z_connected.push_back((hist[0]-transition_limit*0.2)/sqrt(transition_limit*0.2));
        Z_unconnected.push_back((hist[1]-transition_limit*0.8)/sqrt(transition_limit*0.8));
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 0.1."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 0.1." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 0.1." << endl;
    }

    return true;
}

bool Tests::testTwoAndThreePolymerCong(){
    Polymer *template_polymer = new Polymer(1, 3, 0, false);
    Polymer *copy_polymer = new Polymer(2, 2, 1, false);

    Connection * con = new Connection(template_polymer, 1, copy_polymer, 0);
    Conglomerate * conglomerate = new Conglomerate({con});

    bool failed = false;

    if(conglomerate->polymers.size()!=2){
        cout << 1 << endl;
        failed = true;
    }

    if(!(*conglomerate->polymers[0] == *template_polymer) || !(*conglomerate->polymers[1] == *copy_polymer) ){
        cout << 2 << endl;
        failed = true;
    }

    if(conglomerate->polymer_connections[0].size()!=3 || conglomerate->polymer_connections[1].size()!=2){
        cout << 3 << endl;
        failed = true;
    }

    for(int i=0; i<2; i++) {
        if (!conglomerate->polymer_connections[0][2 * i].empty()) {
            cout << 4 << endl;
            failed = true;
        }
    }
    if(!conglomerate->polymer_connections[1][1].empty()){
        cout << 4.2 << endl;
        failed = true;
    }
    if(conglomerate->polymer_connections[0][1].empty()){
        cout << 5 << endl;
        failed = true;
    }
    if(conglomerate->polymer_connections[1][0].empty()){
        cout << 5.5 << endl;
        failed = true;
    }

    if(conglomerate->available_free_sites_list.size()!=2){
        cout << 6 << endl;
        failed = true;
    }
    if(conglomerate->available_free_sites_list[0].size() != 2){
        cout << 7 << endl;
        failed = true;
    }
    if(conglomerate->available_free_sites_list[1].size() != 1){
        cout << 8 << endl;
        failed = true;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        failed = true;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        failed = true;
    }
    if(conglomerate->tail_binding_list.size()!=1){
        cout << 10.5 << endl;
        failed = true;
    }
    if (!conglomerate->connected_neighbours_list.empty()) {
        cout << 11 << endl;
        failed = true;
    }

    if(!conglomerate->unconnected_neighbours_list.empty()){
        cout << 12 << endl;
        failed = true;
    }

    conglomerate->deleteConglomerate();
    delete conglomerate;
    delete template_polymer;
    delete copy_polymer;
    delete con;
    return !failed;
}

bool Tests::testPolymerisationEquilibrium(){
    /* Test the steady state probability of two monomers polymerising
     * Use 5 iterations
     * Use a long time limit
     * Use a template length 2 and 2 free monomers
     * Use a very high G_spec so that the monomers bind to the template and don't unbind
     */

    // If G_bb == G_gen, the monomers should spend approximately equal time connected and unconnected
    // Rate == 1 so average time between events is 1
    vector<double> seeds = {100, 200, 300, 400, 500};
    vector<double> Z;
    vector<double> Z_connected;
    vector<double> Z_unconnected;
    int transition_limit = 1000;
    for(auto & seed : seeds) {

        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = false;
        set_seed = seed;
        set_k = 1;
        set_k0 = 1;
        set_G_bb = -10;
        set_G_spec = -1000;
        set_G_gen = -10;
        set_M_eff = 100;
        set_monomers_family_zero = 0;
        set_monomers_family_one = 2;
        set_template_length = 2;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        System system = System();
        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system.chooseTransition(gen());
            if (system.conglomerates[0]->polymers.size() == 2) {
                //Now connected so we add time to unconnected list
                hist[1] = hist[1] + system.simulation_time - previous_time;
            } else if (system.conglomerates[0]->polymers.size() == 3) {
                //Now unconnected so we add time to connected list
                hist[0] = hist[0] + system.simulation_time - previous_time;
            }
            previous_time = system.simulation_time;
            count++;
        }

        //Transition_limit*k will be the expected number or transitions because the rates are all 1/sec
        Z.push_back((count-transition_limit*set_k)/sqrt(transition_limit*set_k));

        //Both states are equally likely and so we expect that half the time will be spent at each
        //Therefore we use transition_limit*0.5 for both as the expected time
        Z_connected.push_back((hist[0]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        Z_unconnected.push_back((hist[1]-transition_limit*0.5)/sqrt(transition_limit*0.5));

        system.deleteSystem();
    }

    //We take the mean over the 5 tests
    double sum = accumulate(Z.begin(), Z.end(), 0.0);
    double mean = sum/Z.size();

    //If the mean is greater than a specified limit then maybe there is a problem
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 1."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 1." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 1." << endl;
    }

    // If G_bb ==-1 and G_gen == -(1+ln(0.5)), twice as likely to be connected than unconnected
    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();
    for(auto & seed : seeds) {
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = false;
        set_seed = seed;
        set_k = 1;
        set_k0 = 1;
        set_G_bb = -1; //HERE IS THE CHANGE
        set_G_spec = -1000;
        set_G_gen = -(1+log(0.5));//HERE IS THE CHANGE
        set_M_eff = 100;
        set_monomers_family_zero = 0;
        set_monomers_family_one = 2;
        set_template_length = 2;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        System *system = new System();

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates.size() == 1) {
                if (system->conglomerates[0]->polymers.size() == 2) {
                    //Now connected so we add time to unconnected list
                    hist[1] = hist[1] + system->simulation_time - previous_time;
                } else if (system->conglomerates[0]->polymers.size() == 3) {
                    //Now unconnected so we add time to connected list
                    hist[0] = hist[0] + system->simulation_time - previous_time;
                }
            }
            previous_time = system->simulation_time;
            count++;
        }

        // Average time for backbone to form is 1sec
        // Average time for backbone to break is 2sec
        // Average time for any transition to occur is 1.5sec
        // Expected number of transitions = 100000/1.5=66,667 transitions
        Z.push_back((count-transition_limit/1.5)/sqrt(transition_limit/1.5));
        // Time spent connected = 2* time spent unconnected
        // Expected time connected = 66,667
        // Expected time unconnected = 33,333
        Z_connected.push_back((hist[0]-transition_limit/1.5)/sqrt(transition_limit/1.5));
        Z_unconnected.push_back((hist[1]-transition_limit/3)/sqrt(transition_limit/3));
        system->deleteSystem();
        delete system;
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 2."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 2." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 2." << endl;
    }

    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    //If we increase the polymerisation rate, the state probabilities will remain the same but there will be more transitions

    for(auto & seed : seeds) {

        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = false;
        set_seed = seed;
        set_k = 100;//HERE IS THE CHANGE
        set_k0 = 1;
        set_G_bb = -10;
        set_G_spec = -1000;
        set_G_gen = -10;
        set_M_eff = 100;
        set_monomers_family_zero = 0;
        set_monomers_family_one = 2;
        set_template_length = 2;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        System *system = new System();

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates[0]->polymers.size() == 2) {
                //Now connected so we add time to unconnected list
                hist[1] = hist[1] + system->simulation_time - previous_time;
            } else if (system->conglomerates[0]->polymers.size() == 3) {
                //Now unconnected so we add time to connected list
                hist[0] = hist[0] + system->simulation_time - previous_time;
            }
            previous_time = system->simulation_time;
            count++;
        }

        //Transition count will increase as k increases
        Z.push_back((count-transition_limit*set_k)/sqrt(transition_limit*set_k));

        //Still expecting half the time in each state
        Z_connected.push_back((hist[0]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        Z_unconnected.push_back((hist[1]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        system->deleteSystem();
        delete system;
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 1."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 1." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 1." << endl;
    }

    return true;
}

bool Tests::testDimerisationEquilibrium(){
    /* Test the steady state probability of two monomers polymerising
     * Use a template length 1 and 1 free monomer
     */

    // If G_spec == -G_gen, the monomers should spend approximately equal time connected and unconnected
    // Rate == 1 so average time between events is 1
    vector<double> seeds = {101, 201, 301, 401, 501};
    vector<double> Z;
    vector<double> Z_connected;
    vector<double> Z_unconnected;
    int transition_limit = 1000;
    for(auto & seed : seeds) {
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = false;
        set_seed = seed;
        set_k = 1;
        set_k0 = 1;
        set_G_bb = -1;
        set_G_spec = 1;
        set_G_gen = -1;
        set_M_eff = 100;
        set_monomers_family_zero = 0;
        set_monomers_family_one = 1;
        set_template_length = 1;
        set_transition_limit = 1000;


        mt19937 gen(seed);

        System *system = new System();

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates[0]->connections.size() == 1) {
                //Now connected so we add time to unconnected list
                hist[1] = hist[1] + system->simulation_time - previous_time;
            } else if (system->conglomerates[0]->connections.size() == 0) {
                //Now unconnected so we add time to connected list
                hist[0] = hist[0] + system->simulation_time - previous_time;
            }
            previous_time = system->simulation_time;
            count++;
        }
        // Total transitions will increase with the binding rate and the transition limit but rates are 1/sec so no factor
        Z.push_back((count-transition_limit*set_k0)/sqrt(transition_limit*set_k0));
        //We have the two rates the same so the time spent should be equal
        Z_connected.push_back((hist[0]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        Z_unconnected.push_back((hist[1]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        system->deleteSystem();
        delete system;
    }

    double sum = accumulate(Z.begin(), Z.end(), 0.0);
    double mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 3."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 3." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 3." << endl;
    }

    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    //If Gspec == Ggen == -0.5*ln(2) then the connected state is more favourable by a factor 2
    //We expect the monomers are connected for twice as much time as they are separate

    for(auto & seed : seeds) {

        //Initialisations
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = false;
        set_seed = seed;
        set_k = 1;
        set_k0 = 1;
        set_G_bb = -1;
        set_G_spec = -0.5*log(2);//HERE IS THE CHANGE
        set_G_gen = -0.5*log(2);//HERE IS THE CHANGE
        set_M_eff = 100;
        set_monomers_family_zero = 0;
        set_monomers_family_one = 1;
        set_template_length = 1;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        System *system = new System();

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates[0]->connections.size() == 1) {
                //Now connected so we add time to unconnected list
                hist[1] = hist[1] + system->simulation_time - previous_time;
            } else if (system->conglomerates[0]->connections.size() == 0) {
                //Now unconnected so we add time to connected list
                hist[0] = hist[0] + system->simulation_time - previous_time;
            }
            previous_time = system->simulation_time;
            count++;
        }

        // Average time for backbone to form is 1sec
        // Average time for backbone to break is 2sec
        // Average time for any transition to occur is 1.5sec
        // Expected number of transitions = 100000/1.5=66,667 transitions
        Z.push_back((count-transition_limit/1.5)/sqrt(transition_limit/1.5));
        // Time spent connected = 2* time spent unconnected
        // Expected time connected = 66,667
        // Expected time unconnected = 33,333
        Z_connected.push_back((hist[0]-transition_limit/1.5)/sqrt(transition_limit/1.5));
        Z_unconnected.push_back((hist[1]-transition_limit/3)/sqrt(transition_limit/3));
        system->deleteSystem();
        delete system;
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 4."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 4." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 4." << endl;
    }

    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();


    // If we increase k0 then the state probabilities remain the same but the number of transitions will increase
    for(auto & seed : seeds) {

        //Initialisations
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = false;
        set_seed = seed;
        set_k = 1;
        set_k0 = 100;//HERE IS THE CHANGE
        set_G_bb = -1;
        set_G_spec = 1;//HERE IS THE CHANGE
        set_G_gen = -1;
        set_M_eff = 100;
        set_monomers_family_zero = 0;
        set_monomers_family_one = 1;
        set_template_length = 1;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        System *system = new System();

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates[0]->connections.size() == 1) {
                //Now connected so we add time to unconnected list
                hist[1] = hist[1] + system->simulation_time - previous_time;
            } else if (system->conglomerates[0]->connections.size() == 0) {
                //Now unconnected so we add time to connected list
                hist[0] = hist[0] + system->simulation_time - previous_time;
            }
            previous_time = system->simulation_time;
            count++;
        }

        //Same calculations as only k0 makes any change
        Z.push_back((count-transition_limit*set_k0)/sqrt(transition_limit*set_k0));

        Z_connected.push_back((hist[0]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        Z_unconnected.push_back((hist[1]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        system->deleteSystem();
        delete system;
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 4.5."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 4.5." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 4.5." << endl;
    }

    return true;
}

bool Tests::testTailBindEquilibrium(){
    /* Test the steady state probability of a tail binding and unbinding
     * Start with a template of length 3 and a copy polymer of length 2
     * The head of the copy is attached to the tail of the template so the tail of the copy can create and break a 'tail' bond
     * Use a very high G_bb so that the backbone doesn't break
     * Use a very high G_gen so that the head does not detach
     */

    // If G_spec == log(Meff) the tail should spend equal time connected and unconnected
    vector<double> seeds = {103};
    vector<double> Z;
    vector<double> Z_connected;
    vector<double> Z_unconnected;
    int transition_limit = 1000;

    ofstream f_test;
    f_test.open("/Users/nicksimpson/PycharmProjects/MyProject/TestInput.txt");

    for(auto & seed : seeds) {
        vector<double> times;
        vector<double> time_connected;
        //Initialisations
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = false;
        set_seed = seed;
        set_k = 0.1;
        set_k0 = 0.1;
        set_G_bb = -100000;
        set_G_spec = log(100);
        set_G_gen = -1000;
        set_M_eff = 100;
        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_template_length = 0;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        Polymer *template_polymer = new Polymer(-1, 3, 0, false);
        Polymer *copy_polymer = new Polymer(-1, 2, 1, false);

        Connection * con = new Connection(template_polymer, 1, copy_polymer, 0);
        Conglomerate * cong = new Conglomerate({con});
        System *system = new System(cong);

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates.size() == 1) {
                if(system->conglomerates[0]->connections.size()==2){
                    //Now connected so we add time to unconnected list
                    hist[1] = hist[1] + system->simulation_time - previous_time;
                } else if (system->conglomerates[0]->connections.size() == 1) {
                    //Now unconnected so we add time to connected list
                    hist[0] = hist[0] + system->simulation_time - previous_time;
                }
            }
            previous_time = system->simulation_time;
            count++;
            times.push_back(previous_time);
            time_connected.push_back(hist[0]);
        }
        f_test << '[' << times[0];
        for(int i=1; i<times.size(); i++){
            f_test << ',' << times[i] ;
        }
        f_test << ']' << "\n";
        f_test << '[' << time_connected[0];
        for(int i=1; i<time_connected.size(); i++){
            f_test << ',' << time_connected[i] ;
        }
        f_test << ']' << "\n";

        //The effective concentration makes a difference to the rate so must be included
        Z.push_back((count-transition_limit*set_k0*set_M_eff)/sqrt(transition_limit*set_k0*set_M_eff));
        //Equal times are expected for each state
        Z_connected.push_back((hist[0]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        Z_unconnected.push_back((hist[1]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        system->deleteSystem();
        delete system;
    }

    double sum = accumulate(Z.begin(), Z.end(), 0.0);
    double mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 5."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 5." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 5." << endl;
    }

    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    //If we make Gspec = log(Meff/2) then we are more likely to be connected than unconnected
    //Twice as likely connected
    for(auto & seed : seeds) {
        vector<double> times;
        vector<double> time_connected;
        //Initialisations
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = false;
        set_seed = seed;
        set_k = 0.1;
        set_k0 = 0.1;
        set_G_bb = -100000;
        set_G_spec = log(50);//HERE IS THE CHANGE
        set_G_gen = -1000;
        set_M_eff = 100;
        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_template_length = 0;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        Polymer *template_polymer = new Polymer(-1, 3, 0, false);
        Polymer *copy_polymer = new Polymer(-1, 2, 1, false);

        Connection * con = new Connection(template_polymer, 1, copy_polymer, 0);
        Conglomerate * cong = new Conglomerate({con});
        System *system = new System(cong);

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates.size() == 1) {
                if(system->conglomerates[0]->connections.size()==2){
                    //Now connected so we add time to unconnected list
                    hist[1] = hist[1] + system->simulation_time - previous_time;
                } else if (system->conglomerates[0]->connections.size() == 1) {
                    //Now unconnected so we add time to connected list
                    hist[0] = hist[0] + system->simulation_time - previous_time;
                }
            }
            previous_time = system->simulation_time;
            count++;
            times.push_back(previous_time);
            time_connected.push_back(hist[0]);
        }
        f_test << '[' << times[0];
        for(int i=1; i<times.size(); i++){
            f_test << ',' << times[i] ;
        }
        f_test << ']' << "\n";
        f_test << '[' << time_connected[0];
        for(int i=1; i<time_connected.size(); i++){
            f_test << ',' << time_connected[i] ;
        }
        f_test << ']' << "\n";

        //Reduced number of transitions as a more favourable state is available
        Z.push_back((count-transition_limit*set_k0*set_M_eff/1.5)/sqrt(transition_limit*set_k0*set_M_eff/1.5));

        //Twice as likely to be in the connected state as the unconnected state
        Z_connected.push_back((hist[0]-transition_limit/1.5)/sqrt(transition_limit/1.5));
        Z_unconnected.push_back((hist[1]-transition_limit/3)/sqrt(transition_limit/3));
        system->deleteSystem();
        delete system;
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 6."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 6." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 6." << endl;
    }

    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    //If we drop the effective conc (so Gspec = log(2*Meff)) then the unconnected state is more likely
    for(auto & seed : seeds) {
        vector<double> times;
        vector<double> time_connected;
        //Initialisations
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = false;
        set_seed = seed;
        set_k = 0.1;
        set_k0 = 0.1;
        set_G_bb = -100000;
        set_G_spec = log(100);
        set_G_gen = -1000;
        set_M_eff = 50;//HERE IS THE CHANGE
        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_template_length = 0;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        Polymer *template_polymer = new Polymer(-1, 3, 0, false);
        Polymer *copy_polymer = new Polymer(-1, 2, 1, false);

        Connection * con = new Connection(template_polymer, 1, copy_polymer, 0);
        Conglomerate * cong = new Conglomerate({con});
        System *system = new System(cong);

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates.size() == 1) {
                if(system->conglomerates[0]->connections.size()==2){
                    //Now connected so we add time to unconnected list
                    hist[1] = hist[1] + system->simulation_time - previous_time;
                } else if (system->conglomerates[0]->connections.size() == 1) {
                    //Now unconnected so we add time to connected list
                    hist[0] = hist[0] + system->simulation_time - previous_time;
                }
            }
            previous_time = system->simulation_time;
            count++;
            times.push_back(previous_time);
            time_connected.push_back(hist[0]);
        }
        f_test << '[' << times[0];
        for(int i=1; i<times.size(); i++){
            f_test << ',' << times[i] ;
        }
        f_test << ']' << "\n";
        f_test << '[' << time_connected[0];
        for(int i=1; i<time_connected.size(); i++){
            f_test << ',' << time_connected[i] ;
        }
        f_test << ']' << "\n";

        //Transition count stays the same as previous
        Z.push_back((count-transition_limit*set_k0*100/1.5)/sqrt(transition_limit*100*set_k0/1.5));
        //Now the connected state is less likely
        Z_connected.push_back((hist[0]-transition_limit/3)/sqrt(transition_limit/3));
        Z_unconnected.push_back((hist[1]-transition_limit/1.5)/sqrt(transition_limit/1.5));
        system->deleteSystem();
        delete system;
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 7."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 7." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 7." << endl;
    }

    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    //If we increase the binding rate, the number of transitions should increase but the probabilities remain the same
    for(auto & seed : seeds) {
        vector<double> times;
        vector<double> time_connected;
        //Initialisations
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = false;
        set_seed = seed;
        set_k = 0.1;
        set_k0 = 10;//HERE IS THE CHANGE
        set_G_bb = -100000;
        set_G_spec = log(100);
        set_G_gen = -1000;
        set_M_eff = 100;
        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_template_length = 0;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        Polymer *template_polymer = new Polymer(-1, 3, 0, false);
        Polymer *copy_polymer = new Polymer(-1, 2, 1, false);

        Connection * con = new Connection(template_polymer, 1, copy_polymer, 0);
        Conglomerate * cong = new Conglomerate({con});
        System *system = new System(cong);

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates.size() == 1) {
                if(system->conglomerates[0]->connections.size()==2){
                    //Now connected so we add time to unconnected list
                    hist[1] = hist[1] + system->simulation_time - previous_time;
                } else if (system->conglomerates[0]->connections.size() == 1) {
                    //Now unconnected so we add time to connected list
                    hist[0] = hist[0] + system->simulation_time - previous_time;
                }
            }
            previous_time = system->simulation_time;
            count++;
            times.push_back(previous_time);
            time_connected.push_back(hist[0]);
        }
        f_test << '[' << times[0];
        for(int i=1; i<times.size(); i++){
            f_test << ',' << times[i] ;
        }
        f_test << ']' << "\n";
        f_test << '[' << time_connected[0];
        for(int i=1; i<time_connected.size(); i++){
            f_test << ',' << time_connected[i] ;
        }
        f_test << ']' << "\n";

        //Transition count will increase with k0
        Z.push_back((count-transition_limit*set_k0*set_M_eff)/sqrt(transition_limit*set_k0*set_M_eff));

        //Equal chance still for each state
        Z_connected.push_back((hist[0]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        Z_unconnected.push_back((hist[1]-transition_limit*0.5)/sqrt(transition_limit*0.5));
        system->deleteSystem();
        delete system;
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 8."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 8." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 8." << endl;
    }
    f_test.close();
    return true;
}


bool Tests::testGrowthDirections(){
    set_weakened_template_end = false;
    set_parallel_growth = false;

    Polymer * polymer = new Polymer(1, 6, 0, false);
    Polymer * poly = new Polymer(2, 3, 1, false);
    Connection * con = new Connection(polymer, 5, poly, 0);
    vector<Connection *> connections;
    connections.push_back(con);
    Conglomerate * conglomerate = new Conglomerate(connections);

    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_binding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    set_parallel_growth = true;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->tail_binding_list.empty()){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    con->indexes[0] = 0;
    set_parallel_growth = false;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->tail_binding_list.empty()){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }


    set_parallel_growth = true;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_binding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    con->indexes[1] = 1;
    set_parallel_growth = true;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_binding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    con->indexes[1] = 1;
    set_parallel_growth = false;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_binding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    set_parallel_growth = true;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_binding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    con->indexes[0] = 5;
    set_parallel_growth = false;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_unbinding_list.empty()){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty()){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_binding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    set_parallel_growth = true;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_binding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->tail_binding_list.empty() || !conglomerate->head_unbinding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    return true;
}



bool Tests::testWeakEndWithGrowthDirections(){
    //Repeat all of testGrowthDirections but with a weakened end monomer
    set_weakened_template_end = true;
    set_parallel_growth = false;

    Polymer * polymer = new Polymer(1, 6, 0, false);
    Polymer * poly = new Polymer(2, 3, 1, false);
    Connection * con = new Connection(polymer, 5, poly, 0);
    polymer->is_template = true;
    vector<Connection *> connections;
    connections.push_back(con);
    Conglomerate * conglomerate = new Conglomerate(connections);

    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 1 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 2 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->end_unbinding_list.size()!=1){
        cout << 3 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_binding_list.size()!=1){
        cout << 4 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty() || !conglomerate->head_unbinding_list.empty()){
        cout << 5 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 6 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }


    set_parallel_growth = true;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 7 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 8 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 9 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->tail_binding_list.empty()){
        cout << 10 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty() || !conglomerate->end_unbinding_list.empty()){
        cout << 11 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 12 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }


    con->indexes[0] = 0;
    set_parallel_growth = false;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 13 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 14 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 15 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->tail_binding_list.empty()){
        cout << 16 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty() || !conglomerate->end_unbinding_list.empty()){
        cout << 17 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 18 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    set_parallel_growth = true;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 19 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 20 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->end_unbinding_list.size()!=1){
        cout << 21 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_binding_list.size()!=1){
        cout << 22 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty() || !conglomerate->head_unbinding_list.empty()){
        cout << 23 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 24 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    con->indexes[1] = 1;
    set_parallel_growth = false;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 25 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 26 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_unbinding_list.size()!=1){
        cout << 27 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_binding_list.size()!=1){
        cout << 28 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->tail_binding_list.empty() || !conglomerate->tail_unbinding_list.empty() || !conglomerate->end_unbinding_list.empty()){
        cout << 29 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 30 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    set_parallel_growth = true;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 31 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 32 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->end_unbinding_list.size()!=1){
        cout << 33 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_binding_list.size()!=1){
        cout << 34 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty() || !conglomerate->tail_unbinding_list.empty() || !conglomerate->head_unbinding_list.empty()){
        cout << 35 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 36 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    con->indexes[0] = 5;
    set_parallel_growth = false;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 37 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 38 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_unbinding_list.empty()){
        cout << 39 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->head_binding_list.empty()){
        cout << 40 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_binding_list.size()!=1){
        cout << 41 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->end_unbinding_list.size()!=1){
        cout << 42 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 43 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }

    set_parallel_growth = true;
    conglomerate->updateConglomerate();
    if(conglomerate->available_free_sites_list[0].size()!=5){
        cout << 44 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->available_free_sites_list[1].size()!=2){
        cout << 45 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->tail_unbinding_list.size()!=1){
        cout << 46 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(conglomerate->head_binding_list.size()!=1){
        cout << 47 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->tail_binding_list.empty() || !conglomerate->head_unbinding_list.empty()){
        cout << 48 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    if(!conglomerate->unconnected_neighbours_list.empty() || !conglomerate->connected_neighbours_list.empty()){
        cout << 49 << endl;
        conglomerate->deleteConglomerate();
        delete conglomerate;
        delete polymer;
        delete poly;
        return false;
    }
    return true;
}


bool Tests::testEndMonomer(){
    set_parallel_growth = true;
    set_end_monomer = true;
    set_weakened_template_end = false;
    //Repeat the head and tail experiments below but with the weakened bond active

    // 1 free monomer and a template length 1
    // Every interaction is with the `end' monomer so Gbb and Gspec and Ggen should make no difference
    // If G_end = ln(2) then we will be unconnected twice as much as connected
    vector<double> seeds = {101, 201, 301, 401, 501};
    vector<double> Z;
    vector<double> Z_connected;
    vector<double> Z_unconnected;
    int transition_limit = 1000;
    for(auto & seed : seeds) {
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = true;
        set_seed = seed;
        set_volume = 1;
        set_k = 1;
        set_k0 = 1;
        set_G_bb = -1;
        set_G_spec = 1;
        set_G_gen = -(1+log(0.5));
        set_M_eff = 100;
        set_G_end = log(2);
        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_end_monomers_family_zero = 0;
        set_end_monomers_family_one = 1;
        set_template_length = 1;
        set_transition_limit = 1000;


        mt19937 gen(seed);

        System *system = new System();

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates[0]->polymers.size() == 2) {
                //Now connected so we add time to unconnected list
                hist[1] = hist[1] + system->simulation_time - previous_time;
            } else if (system->conglomerates[0]->polymers.size() == 1) {
                //Now unconnected so we add time to connected list
                hist[0] = hist[0] + system->simulation_time - previous_time;
            }
            previous_time = system->simulation_time;
            count++;
        }

        // Total transitions will increase with the binding rate and the transition limit but rates are 1/sec so no factor
        Z.push_back((count-transition_limit*set_k0*4/3)/sqrt(transition_limit*set_k0*4/3));
        //We have the two rates the same so the time spent should be equal
        Z_connected.push_back((hist[0]-transition_limit*0.33)/sqrt(transition_limit*0.33));
        Z_unconnected.push_back((hist[1]-transition_limit*0.66)/sqrt(transition_limit*0.66));
    }

    double sum = accumulate(Z.begin(), Z.end(), 0.0);
    double mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 0.1."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 0.1." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 0.1." << endl;
    }


    //Template length 2
    //Copy length 2 attached by the tail on index 1
    //Tail should flip between the end monomer
    // If G_end is ln(2*Meff) we are twice as likely to be unconnected
    seeds = {103, 203, 303, 403, 503};
    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    for(auto & seed : seeds) {
        //Initialisations
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = true;
        set_seed = seed;
        set_k = 0.1;
        set_k0 = 0.1;
        set_G_bb = -100000000;
        set_G_spec = -1000000;
        set_G_gen = -1000000;
        set_M_eff = 10;
        set_G_end = log(20);
        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_end_monomers_family_zero = 0;
        set_end_monomers_family_one = 0;
        set_template_length = 0;
        set_transition_limit = 1000;

        mt19937 gen(seed);

        Polymer *template_polymer = new Polymer(-1, 2, 0, true);
        template_polymer->is_template = true;
        Polymer *copy_polymer = new Polymer(-1, 2, 1, true);

        Connection * con = new Connection(template_polymer, 1, copy_polymer, 1);
        Conglomerate * cong = new Conglomerate({con});
        System *system = new System(cong);

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates.size() == 1) {
                if(system->conglomerates[0]->connections.size()==2){
                    //Now connected so we add time to unconnected list
                    hist[1] = hist[1] + system->simulation_time - previous_time;
                } else if (system->conglomerates[0]->connections.size() == 1) {
                    //Now unconnected so we add time to connected list
                    hist[0] = hist[0] + system->simulation_time - previous_time;
                }
            }
            previous_time = system->simulation_time;
            count++;

        }

        //The effective concentration makes a difference to the rate so must be included
        Z.push_back((count-transition_limit*set_k0*set_M_eff*4/3)/sqrt(transition_limit*set_k0*set_M_eff*4/3));
        //Equal times are expected for each state
        Z_connected.push_back((hist[0]-transition_limit*0.33)/sqrt(transition_limit*033));
        Z_unconnected.push_back((hist[1]-transition_limit*0.66)/sqrt(transition_limit*0.66));
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 0.2."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 0.2." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 0.2." << endl;
    }


    //Template length 2
    //Copy length 2 attached by the tail on index 1
    //Tail should flip between the end monomer
    // If G_end is ln(2*Meff) we are twice as likely to be unconnected
    // halve the rate, halve the transitions
    seeds = {103, 203, 303, 403, 503};
    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    for(auto & seed : seeds) {
        //Initialisations
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = true;
        set_seed = seed;
        set_k = 0.1;
        set_k0 = 0.05;
        set_G_bb = -100000000;
        set_G_spec = -1000000;
        set_G_gen = -1000000;
        set_M_eff = 10;
        set_G_end = log(20);
        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_end_monomers_family_zero = 0;
        set_end_monomers_family_one = 0;
        set_template_length = 0;
        set_transition_limit = 10;

        mt19937 gen(seed);

        Polymer *template_polymer = new Polymer(-1, 2, 0, true);
        template_polymer->is_template = true;
        Polymer *copy_polymer = new Polymer(-1, 2, 1, true);

        Connection * con = new Connection(template_polymer, 1, copy_polymer, 1);
        Conglomerate * cong = new Conglomerate({con});
        System *system = new System(cong);

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates.size() == 1) {
                if(system->conglomerates[0]->connections.size()==2){
                    //Now connected so we add time to unconnected list
                    hist[1] = hist[1] + system->simulation_time - previous_time;
                } else if (system->conglomerates[0]->connections.size() == 1) {
                    //Now unconnected so we add time to connected list
                    hist[0] = hist[0] + system->simulation_time - previous_time;
                }
            }
            previous_time = system->simulation_time;
            count++;

        }

        //The effective concentration makes a difference to the rate so must be included
        Z.push_back((count-transition_limit*set_k0*set_M_eff*4/3)/sqrt(transition_limit*set_k0*set_M_eff*4/3));
        //Equal times are expected for each state
        Z_connected.push_back((hist[0]-transition_limit*0.33)/sqrt(transition_limit*033));
        Z_unconnected.push_back((hist[1]-transition_limit*0.66)/sqrt(transition_limit*0.66));
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 0.3."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 0.3." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 0.3." << endl;
    }

    // 1 free monomer and a template length 1
    // Every interaction is with the `end' monomer so Gbb and Gspec and Ggen should make no difference
    // If G_end = ln(2) then we will be unconnected twice as much as connected
    // Now change G_end = ln(4) so we spend more time unconnected but increase Ggen loads so we want to be connected
    seeds = {101, 201, 301, 401, 501};
    Z.clear();
    Z_connected.clear();
    Z_unconnected.clear();

    for(auto & seed : seeds) {
        set_template_indestructible = true;
        set_monomer_count_is_constant = false;
        set_no_rebinding = false;
        set_weakened_template_end = true;
        set_seed = seed;
        set_k = 1;
        set_k0 = 1;
        set_G_bb = -1;
        set_G_spec = 1;
        set_G_gen = -100;
        set_M_eff = 100;
        set_G_end = log(4);
        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_end_monomers_family_zero = 0;
        set_end_monomers_family_one = 1;
        set_template_length = 1;
        set_transition_limit = 1000;


        mt19937 gen(seed);

        System *system = new System();

        int count = 0;

        vector<double> hist = {0, 0};
        double previous_time = 0;

        while (previous_time < transition_limit) {
            system->chooseTransition(gen());
            if (system->conglomerates[0]->polymers.size() == 2) {
                //Now connected so we add time to unconnected list
                hist[1] = hist[1] + system->simulation_time - previous_time;
            } else if (system->conglomerates[0]->polymers.size() == 1) {
                //Now unconnected so we add time to connected list
                hist[0] = hist[0] + system->simulation_time - previous_time;
            }
            previous_time = system->simulation_time;
            count++;
        }

        // Total transitions will increase with the binding rate and the transition limit but rates are 1/sec so no factor
        Z.push_back((count-transition_limit*set_k0/0.625)/sqrt(transition_limit*set_k0/0.625));
        //We have the two rates the same so the time spent should be equal
        Z_connected.push_back((hist[0]-transition_limit*0.2)/sqrt(transition_limit*0.2));
        Z_unconnected.push_back((hist[1]-transition_limit*0.8)/sqrt(transition_limit*0.8));
    }

    sum = accumulate(Z.begin(), Z.end(), 0.0);
    mean = sum/Z.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the number of transitions in the simulation 0.1."<< endl;
    }

    sum = accumulate(Z_connected.begin(), Z_connected.end(), 0.0);
    mean = sum/Z_connected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent connected in simulation 0.1." << endl;
    }

    sum = accumulate(Z_unconnected.begin(), Z_unconnected.end(), 0.0);
    mean = sum/Z_unconnected.size();
    if(abs(mean)>=1.65){
        cout << "Confidence level lower than expected." << endl;
        cout << "Z=" << mean << " for the time spent unconnected in simulation 0.1." << endl;
    }

    return true;
}

bool Tests::testCorrectMonomerBinding(){

    set_volume = 1;
    set_k0 = 1;
    set_end_monomer = true;
    set_parallel_growth = true;

    for(int j=1; j<3;j++){
        set_template_length = j;

        set_monomers_family_zero = 0;
        set_monomers_family_one = 0;
        set_end_monomers_family_zero = 0;
        set_end_monomers_family_one = 1;

        System system = System();

        if(system.end_external_rate != 1){
            cout << "1" << endl;
            cout << system.end_external_rate << endl;
            return false;
        }
        if(system.external_connection_rate != 0){
            cout << "2" << endl;
            return false;
        }
        for(int i=0; i<system.transition_rates.size(); i++){
            if(system.transition_rates[i]!=0){
                cout << "3." << i << endl;
                return false;
            }
        }

        set_monomers_family_zero = 1;
        set_monomers_family_one = 1;
        set_end_monomers_family_zero = 1;
        set_end_monomers_family_one = 0;

        system.updateRates(0);

        if(system.end_external_rate != 1){
            cout << "4" << endl;
            return false;
        }

        if(system.external_connection_rate != 0){
            cout << "5.1" << endl;
            return false;
        }


        for(int i=0; i<system.transition_rates.size(); i++){
            if(system.transition_rates[i]!=0){
                cout << "6." << i << endl;
                return false;
            }
        }
        system.deleteSystem();
    }

    return true;
}
