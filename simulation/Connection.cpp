//
// Created by Nicholas Simpson on 03/02/2022.
//

#include "Connection.h"

bool Connection::operator==(Connection c){

    if(*polymers_in_connection[0] == *c.polymers_in_connection[0] && *polymers_in_connection[1] == *c.polymers_in_connection[1]){
        if(indexes[0] == c.indexes[0] && indexes[1] == c.indexes[1]){
            return true;
        }
    }

    if(*polymers_in_connection[1] == *c.polymers_in_connection[0] && *polymers_in_connection[0] == *c.polymers_in_connection[1]){
        if(indexes[1] == c.indexes[0] && indexes[0] == c.indexes[1]){
            return true;
        }
    }



    return false;
}