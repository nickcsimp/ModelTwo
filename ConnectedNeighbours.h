//
// Created by Nicholas Simpson on 03/02/2022.
//

#ifndef MODELTWO_CONNECTEDNEIGHBOURS_H
#define MODELTWO_CONNECTEDNEIGHBOURS_H


#include "Polymer.h"

class ConnectedNeighbours {
public:
    Polymer * polymer; //This identifies which polymer the neighbours are part of
    int index; //This identifies the lower index of the connected neighbours (The other index is index+1)

    ConnectedNeighbours(Polymer * p, int ind){
        polymer = p;
        index = ind;
    };

    bool operator==(ConnectedNeighbours c);
};


#endif //MODELTWO_CONNECTEDNEIGHBOURS_H
