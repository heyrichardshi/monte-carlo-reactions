#include "reaction.hpp"

using namespace std;

Reaction::Reaction(double k, vector<int> reactants, vector<int> a, vector<int> products, vector<int> b) :
_k(k), _reactants(reactants), _a(a), _products(products), _b(b), _dead(false) { }

int Reaction::exists(int species) {
    int i = -1;
    for (auto j : _reactants) {
        if (j == species) {
            return i;
        }
        i--;
    }
    i = 1;
    for (auto j : _products) {
        if (j == species) {
            return i;
        }
        i++;
    }
    return 0;
}

int Reaction::coefficient(int species) {
    int index = exists(species);
    if (!index) {
        return 0;
    }
    if (index < 0) {
        return _a[abs(index) - 1];
    }
    return _b[index - 1];
}

int Reaction::reactants() {
    return _reactants.size();
}

int Reaction::products() {
    return _products.size();
}

double Reaction::k() {
    return _k;
}

int Reaction::reactant(int index) {
    return _reactants[index];
}

int Reaction::product(int index) {
    return _products[index];
}

bool Reaction::dead() {
    return !_dead;
}

void Reaction::kill() {
    _dead = true;
}

void Reaction::unkill() {
    _dead = false;
}