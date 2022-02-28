//
// Created by Nicholas Simpson on 28/02/2022.
//

#ifndef MODELTWO_SETTINGS_H
#define MODELTWO_SETTINGS_H


class settings {
public:
    bool template_indestructible = true;
    bool monomer_count_is_constant = false; //Only useful if no_rebinding = true
    bool no_rebinding = false;
    double seed = 200;
    double k = 1;//Polymerisation rate
    double k0 = 1;//Binding rate
    double G_bb = -1;//Backbone forming free energy
    double G_spec = 2;//Specific bond forming free energy
    double G_gen = -2;//Generic bond forming free energy
    double M_eff = 100;//Effective concentration of monomers in zipping
    int monomers_family_zero = 0;
    int monomers_family_one = 100;
    int template_length = 10;
    int transition_limit = 100000;
};


#endif //MODELTWO_SETTINGS_H
