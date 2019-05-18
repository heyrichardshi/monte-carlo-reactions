#ifndef MONTE_CARLO_REACTIONS_SIMULATION_HPP
#define MONTE_CARLO_REACTIONS_SIMULATION_HPP

#include <vector>
#include <cmath>

#include "reaction.cpp"

using namespace std;

class Simulation {
    int _species;
    vector<Reaction> _reactions;
    double _A;
    double _M1;
    double _M2;
    vector<int> _N; // number of each species
    vector<int> _Ni; // initial
    vector<double> _dNdt; // for each species
    vector<double> _dNidt; // initial
    vector<vector<int>> _Ntrace;
    vector<double> _time;
public:
    Simulation(int n);
    void addReaction(Reaction r);
    int reactions();
    void steadyState();
    void run();
    void cycle();
    double calcR(int alpha);
    double getN(int species);
    double calcdNdt(int species);
    void printN();
};


#endif //MONTE_CARLO_REACTIONS_SIMULATION_HPP
