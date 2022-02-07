//
// Created by Nicholas Simpson on 03/02/2022.
//

#ifndef MODELTWO_UNCONNECTEDNEIGHBOURS_H
#define MODELTWO_UNCONNECTEDNEIGHBOURS_H


#include "Polymer.h"

class UnconnectedNeighbours {
public:
    Polymer * polymer_one; //This is the tail end of one polymer (index will be getLength()-1)
    Polymer * polymer_two; //This is the head end of one polymer (index will be 0)

    UnconnectedNeighbours(Polymer * p_one, Polymer * p_two){
        polymer_one = p_one;
        polymer_two = p_two;
    }
};


#endif //MODELTWO_UNCONNECTEDNEIGHBOURS_H
