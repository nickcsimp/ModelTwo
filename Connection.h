//
// Created by Nicholas Simpson on 03/02/2022.
//

#ifndef MODELTWO_CONNECTION_H
#define MODELTWO_CONNECTION_H

#include <vector>

#include "Polymer.h"


class Connection {
public:

    vector<Polymer *> polymers_in_connection; //(no particular order)
    vector<int> indexes;

    Connection(Polymer * p_one, int i_one, Polymer * p_two, int i_two){
        polymers_in_connection.push_back(p_one);
        polymers_in_connection.push_back(p_two);
        indexes.push_back(i_one);
        indexes.push_back(i_two);
    }

    bool operator==(Connection c);
};


#endif //MODELTWO_CONNECTION_H
