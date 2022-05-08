//
// Created by Nicholas Simpson on 03/02/2022.
//

#ifndef MODELTWO_POLYMER_H
#define MODELTWO_POLYMER_H
#include <vector>

using namespace std;

class Polymer {
public:
    int index; //Identifies the polymer
    int length; //Length used to know how many monomers are in the polymer
    int family; //Family needed to identify which family of monomers are in the polymer
    bool is_template;
    bool has_end_monomer;

    Polymer(int ind, int len, int fam, bool has_end){
        index=ind;
        length=len;
        family=fam;
        is_template = false;
        has_end_monomer = has_end;
    }

    bool operator==(Polymer p);
};


#endif //MODELTWO_POLYMER_H
