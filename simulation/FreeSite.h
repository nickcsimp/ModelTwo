//
// Created by Nicholas Simpson on 03/02/2022.
//

#ifndef MODELTWO_FREESITE_H
#define MODELTWO_FREESITE_H


#include "Polymer.h"

class FreeSite {
public:
    Polymer * polymer; //Identifies the polymer
    int index; //Identifies the location on the polymer

    FreeSite(Polymer * p, int ind){
        polymer = p;
        index = ind;
    }
};


#endif //MODELTWO_FREESITE_H
