//
// Created by Nicholas Simpson on 16/02/2022.
//

#ifndef MODELTWO_TESTS_H
#define MODELTWO_TESTS_H

#include <iostream>
#include <fstream>
#include "Polymer.h"
#include "Conglomerate.h"
#include "System.h"
#include "settings.h"
#include "unistd.h"

class Tests {
public:
    void run();
    bool testConglomerateInitialissation();
    bool testConglomerateUpdate();
    bool testConglomerateAddConnection();
    bool testConglomerate();
    bool testMiddleTailUnbinding();
    bool testEndTailUnbinding();
    bool testWeakenedEndBond();
    bool testPolymerisationEquilibrium();
    bool testDimerisationEquilibrium();
    bool testTailBindEquilibrium();
    bool testTwoAndThreePolymerCong();
    bool testGrowthDirections();
    bool testWeakEndWithGrowthDirections();
    bool testEndMonomer();
    bool testCorrectMonomerBinding();
};


#endif //MODELTWO_TESTS_H
