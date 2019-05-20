#ifndef MONTE_CARLO_REACTIONS_REACTION_HPP
#define MONTE_CARLO_REACTIONS_REACTION_HPP

#include <vector>
#include <cmath>

using namespace std;

class Reaction {
    double _k; // inverse seconds
    vector<int> _reactants;
    vector<int> _products;
    vector<int> _a; // reactant coefficients;
    vector<int> _b; // reactant coefficients;
    bool _dead;
public:
    Reaction(double k, vector<int> reactants, vector<int> a, vector<int> products, vector<int> b);
    int exists(int species);
    int coefficient(int species);
    int reactants();
    int products();
    double k();
    int reactant(int index);
    int product(int index);
    bool dead();
    void kill();
    void unkill();
};


#endif //MONTE_CARLO_REACTIONS_REACTION_HPP
