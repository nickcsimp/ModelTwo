//
// Created by Nicholas Simpson on 03/02/2022.
//

#include "ConnectedNeighbours.h"



bool ConnectedNeighbours::operator==(ConnectedNeighbours c){
    if(*polymer == *c.polymer && index == c.index){
        return true;
    }
    return false;
}
