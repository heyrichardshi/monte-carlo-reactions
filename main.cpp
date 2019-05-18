#include <iostream>
#include <fstream>

#include "simulation.cpp"

int main() {
    Simulation s(4);

    vector<int> species0;
    vector<int> species1;
    vector<int> species2;
    vector<int> species3;
    vector<int> coeff1;

    species0.push_back(0);
    species1.push_back(1);
    species2.push_back(2);
    species3.push_back(3);
    coeff1.push_back(1);

    Reaction rxn0(5, species0, coeff1, species1, coeff1);
    Reaction rxn1(4, species1, coeff1, species2, coeff1);
    Reaction rxn2(3, species2, coeff1, species1, coeff1);
    Reaction rxn3(3, species1, coeff1, species3, coeff1);

    s.addReaction(rxn0);
    s.addReaction(rxn1);
    s.addReaction(rxn2);
    s.addReaction(rxn3);

    s.run();

    return 0;
}