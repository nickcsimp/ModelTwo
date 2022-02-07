#include <iostream>
#include <random>

using namespace std;

int main() {
    cout << "Hello, World!" << std::endl;

    mt19937 gen(200);
    vector<double> v_one;

    for(int i=0; i<100; i++){
        v_one.push_back(gen());
    }

    mt19937 gener(200);
    vector<double> v_two;

    for(int i=0; i<100; i++){
        v_two.push_back(gener());
    }

    for(int i=0; i<100; i++){
        cout << v_one[i] << "   " << v_two[i] << endl;
    }

    //Initialisations
    int count = 0;
    int transition_limit = 1000;
    while(count<transition_limit){
        //chooseTransition
    }
    return 0;
}
