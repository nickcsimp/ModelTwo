//
// Created by Nicholas Simpson on 03/02/2022.
//

#include "Polymer.h"

bool Polymer::operator==(Polymer p) {
    if (p.index == index && p.length == length && p.family == family) {
        return true;
    }
    return false;
}