#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include "Polymer.h"
#include "System.h"
#include "Tests.h"
#include "settings.h"

using namespace std;

bool set_template_indestructible;
bool set_monomer_count_is_constant; //Only useful if no_rebinding = true
bool set_no_rebinding;
double set_seed;
double set_k;//Polymerisation rate
double set_k0;//Binding rate
double set_G_bb;//Backbone forming free energy
double set_G_spec;//Specific bond forming free energy
double set_G_gen;//Generic bond forming free energy
double set_M_eff;//Effective concentration of monomers in zipping
double set_G_end;
int set_monomers_family_zero;
int set_monomers_family_one;
int set_template_length;
int set_transition_limit;
int set_time_limit;
int set_polymer_limit;
bool set_run_tests;
bool set_make_animated_histogram;
bool set_make_final_histogram;
bool set_make_average_length_graph;
bool set_make_length_distribution_plots;
bool set_make_images;
double set_volume;
bool set_weakened_template_end;
bool set_parallel_growth;

void read_input(string filename){
    string input(filename);
    string line;
    string delimiter = ",";
    ifstream fin;
    fin.open(input);
    if (fin.fail()) {
        // file could not be opened
        cerr << "Input file could not be opened. \n Exiting... \n\n";
        exit(EXIT_FAILURE);
    }
    int line_counter = 0;
    while (!fin.eof()) {
        fin >> line;
        line += ", ";
        if(line_counter == 1){
            size_t pos = 0;
            string token;
            int variable_count = 0;
            while ((pos = line.find(delimiter)) != std::string::npos) {
                pos = line.find(delimiter);
                token = line.substr(0, pos);
                if(variable_count == 0){
                    set_seed = stod(token);
                } else if(variable_count == 1){
                    set_k = stod(token);
                } else if(variable_count == 2){
                    set_k0 = stod(token);
                } else if(variable_count == 3){
                    set_G_bb = stod(token);
                } else if(variable_count == 4){
                    set_G_spec = stod(token);
                } else if(variable_count == 5){
                    set_G_gen = stod(token);
                } else if(variable_count == 6){
                    set_G_end = stod(token);
                } else if(variable_count == 7){
                    set_M_eff = stod(token);
                } else if(variable_count == 8){
                    set_monomers_family_zero = stoi(token);
                } else if(variable_count == 9){
                    set_monomers_family_one = stoi(token);
                } else if(variable_count == 10){
                    set_template_length = stoi(token);
                } else if(variable_count == 11){
                    set_volume = stod(token);
                }else if(variable_count == 12){
                    set_time_limit = stod(token);
                } else if(variable_count == 13){
                    set_transition_limit = stoi(token);
                } else if(variable_count == 14){
                    set_polymer_limit = stoi(token);
                } else if(variable_count == 15){
                    set_template_indestructible = (token == "TRUE");
                } else if(variable_count == 16){
                    set_monomer_count_is_constant= (token == "TRUE");
                } else if(variable_count == 17){
                    set_no_rebinding = (token == "TRUE");
                } else if(variable_count == 18){
                    set_parallel_growth = (token == "TRUE");
                } else if(variable_count == 19){
                    set_weakened_template_end = (token == "TRUE");
                } else if(variable_count == 20){
                    set_run_tests = (token == "TRUE");
                } else if(variable_count == 21){
                    set_make_animated_histogram = (token == "TRUE");
                } else if(variable_count == 22){
                    set_make_final_histogram = (token == "TRUE");
                } else if(variable_count == 23){
                    set_make_average_length_graph = (token == "TRUE");
                } else if(variable_count == 24){
                    set_make_length_distribution_plots = (token == "TRUE");
                } else if(variable_count == 25){
                    set_make_images = (token == "TRUE");
                }
                variable_count++;
                line.erase(0, pos + delimiter.length());
            }
        }
        line_counter++;
    }

}

int main(int argc, char *argv[]) {

    if(argc!=2){
        cout << "No input file." << endl;
        exit(EXIT_FAILURE);
    }

    string input_file_name = argv[1];
    read_input("/Users/nicksimpson/CLionProjects/ModelTwo/inputlist.csv");

    if(set_run_tests){
        Tests tests;
        tests.run();
    }

    double seed = set_seed;
    mt19937 gen(seed);

    System * system = new System();

    int count = 0;
    int transition_limit = set_transition_limit;
    bool transitions_possible = true;
    ofstream f_hist;

    if(set_make_animated_histogram || set_make_final_histogram || set_make_average_length_graph) {
        f_hist.open("../histogram.txt", std::ios_base::app);
    }

    int polymer_count = 0;

    while(transitions_possible && polymer_count<set_polymer_limit){
        transitions_possible = system->chooseTransition(gen());
        count ++;
        polymer_count = 0;//Reset count
        for(int i=0; i<system->lengths.size(); i++){ //Loop all polymer lengths
            polymer_count = polymer_count + system->lengths[i]; //Count how many polymers there are
        }
        polymer_count = polymer_count - set_monomers_family_zero-set_monomers_family_one-1; // Remove initial monomers and template polymer
    }
    if(set_make_final_histogram) {
        if(f_hist.is_open()){
            f_hist << "\"Gbb=" << set_G_bb << ", Ggen=" << set_G_gen << ", Gspec=" << set_G_spec << ", Gend=" << set_G_end << ", Rebinding=" << !set_no_rebinding << "\"\n";
            if(set_monomer_count_is_constant){
                lengths[0] = lengths[0]-set_monomers_family_one-set_monomers_family_zero;
            }
            f_hist << '[';
            for(int i = 0; i< system->lengths.size()-1; i++){
                f_hist << system->lengths[i] << ", ";
            }
            f_hist << system->lengths[system->lengths.size()-1] << ']' << "\n";
        }
    }

    f_hist.close();

    if(set_make_length_distribution_plots) {
        int full_length_count = system->lengths[system->lengths.size()-1]-1;
        int length_count = 0; //Initialise total length count
        polymer_count = 0;//Reset polymer count
        for(int i=0; i<system->lengths.size(); i++){ //Loop all polymer lengths
            polymer_count = polymer_count + system->lengths[i]; //Count how many polymers there are
            length_count = length_count + system->lengths[i]*(i+1); //Count the lengths of all the polymers
        }
        polymer_count = polymer_count - set_monomers_family_zero-set_monomers_family_one-1; // Remove initial monomers and template polymer
        length_count = length_count - set_monomers_family_zero-set_monomers_family_one-set_template_length; //Remove initial monomers and template polymer

        double average_length = double(length_count) / double(polymer_count);

        ofstream myfile;
        myfile.open("/Users/nicksimpson/CLionProjects/ModelTwo/LengthDist.csv", std::ios_base::app);
        if (myfile.is_open()) {
            myfile << set_G_gen << ',' << set_G_bb << ',' << average_length << "\n";
        }
        myfile.close();
    }

    if(set_make_images) {
        ofstream fw("/Users/nicksimpson/PycharmProjects/MyProject/input.txt", ofstream::out);

        if (fw.is_open()) {
            //Creating all the nodes and joining polymers
            for (int cong = 0; cong < system->conglomerates.size(); cong++) {
                //We only want to view the things that are interesting
                //If there is a single monomer or a dimer then we ignore it
                if (system->conglomerates[cong]->polymers.size() == 1 &&
                    system->conglomerates[cong]->polymers[0]->length == 1) {
                    //TODO make a monomer image?
                    //Monomer
                } else if (system->conglomerates[cong]->polymers.size() == 2 &&
                           system->conglomerates[cong]->polymers[0]->length == 1 &&
                           system->conglomerates[cong]->polymers[1]->length == 1) {
                    //Dimer
                    //TODO make a dimer image?
                } else {

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
                    for (int poly = 0; poly < system->conglomerates[cong]->polymer_connections.size(); poly++) {
                        for (int mon = 0; mon < system->conglomerates[cong]->polymer_connections[poly].size(); mon++) {
                            if (!system->conglomerates[cong]->polymer_connections[poly][mon].empty()) {
                                if (is_first) {
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
                    for (int poly = 0; poly < system->conglomerates[cong]->polymer_connections.size(); poly++) {
                        for (int mon = 0; mon < system->conglomerates[cong]->polymer_connections[poly].size(); mon++) {
                            if (!system->conglomerates[cong]->polymer_connections[poly][mon].empty()) {
                                if (is_first) {
                                    //We need to not have a comma before the first connection
                                    fw << '('
                                       << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[0];
                                    fw << ','
                                       << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[1]
                                       << ")";
                                    is_first = false;
                                } else {
                                    fw << ", ("
                                       << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[0];
                                    fw << ','
                                       << system->conglomerates[cong]->polymer_connections[poly][mon][0]->indexes[1]
                                       << ")";
                                }
                            }
                        }
                    }
                    fw << ']' << endl;
                }
            }

            fw.close();
        }
    }


    if(set_make_images || set_make_final_histogram || set_make_animated_histogram || set_make_average_length_graph) {
        ofstream f_python_main("main.py", ofstream::out);
        f_python_main << "if __name__ == '__main__':" << "\n";
        f_python_main << "    import sys" << "\n";
        f_python_main << "    sys.path.append('../dataAnalysis')" << "\n";
        if (set_make_images) {
            f_python_main << "    import Images" << "\n";
        }
        if (set_make_final_histogram || set_make_average_length_graph) {
            f_python_main << "    import Histogram" << "\n";
        }
        if (set_make_animated_histogram) {
            f_python_main << "    import AnimatedHistogram" << "\n";
        }

        if (set_make_images) {
            f_python_main << "\n" << "    Images.create_plots('input.txt', 'figures/Images')";
        }
        if (set_make_final_histogram) {
            f_python_main << "\n" << "    Histogram.create_histogram('histogram.txt')";
        }
        if (set_make_average_length_graph) {
            f_python_main << "\n" << "    Histogram.create_average_length_graph('histogram.txt')";
        }
        if (set_make_animated_histogram) {
            f_python_main << "\n" << "    AnimatedHistogram.animate_histogram('histogram.txt')";
        }
    }

    system->deleteSystem();
    delete system;
    return 0;
}
