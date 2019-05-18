#include <cmath>

#include "simulation.hpp"
#include "random.hpp"

using namespace std;

Simulation::Simulation(int n) : _species(n), _A(0.5), _M1(20), _M2(1000) {
    // STEP 1
    for (int i = 0; i < n; i++) {
        _N.push_back(RandomInt(1000, 2000));
        _Ni.push_back(_N[i]); // TODO: should reset each cycle?
        _dNdt.push_back(0);
        _dNidt.push_back(0);
    }
}

void Simulation::addReaction(Reaction r) {
    _reactions.push_back(r);
}

int Simulation::reactions() {
    return _reactions.size();
}

void Simulation::steadyState() {
    //derp;
    cout << "steady state" << endl;
    for (int i = 0; i < _time.size(); i++) {
        cout << _time[i] << " ";
        for (auto j : _Ntrace[i]) {
            cout << j << " ";
        }
        cout << endl;
    }
}

void Simulation::run() {
    _time.push_back(0);
    _Ntrace.push_back(_N);

    int changeless = 0;

    for (int i = 0; i < _M2; i++) {

        for (int j = 0; j < _species; j++) {
            _dNidt[j] = calcdNdt(j); //TODO: reset each cycle (or stay init)
            _Ni[j] = _N[j];
        }

        printN();

        for (int j = 0; j <= changeless; j++) {
            cycle();
        }

        // STEP 6
        vector<bool> restricted;
        vector<double> change;
        for (int j = 0; j < _species; j++) {
            restricted.push_back(abs(_N[j] - _Ni[j]) < 3);
            change.push_back(abs(2 * (_dNdt[j] - _dNidt[j]) / (_dNdt[j] + _dNidt[j])));
//            change.push_back(abs((_dNdt[j] - _dNidt[j]) / _dNidt[j]));
            cout << "species " << j << " | restricted: " << restricted[j] << " | change = " << change[j] << endl;
        }

        int smallest = -1;
        for (int j = 0; j < _species; j++) {
            if (restricted[j]) {
                continue;
            }
            if (smallest == -1) {
                smallest = j;
            }
            if (change[j] < change[smallest]) {
                smallest = j;
            }
        }
        cout << "smallest " << smallest << endl;

        if (smallest == -1) {
            if (changeless == 3) {
                steadyState();
                return;
            }
            changeless++;
            i--;
            continue;
        }
        else {
            changeless = 0;
        }
        double deltat = 2 * (_N[smallest] - _Ni[smallest]) / (_dNdt[smallest] + _dNidt[smallest]);
        cout << "delta_t = " << "2 * (" << _N[smallest] << "-" << _Ni[smallest] << ") / (" << _dNdt[smallest] << "+" << _dNidt[smallest] << ") = " << deltat << endl;
        _time.push_back(_time[_time.size() - 1] + deltat);
        _Ntrace.push_back(_N);

        cout << "t ==> " << _time[_time.size() - 1] << endl;
    }
    cout << "M2 reached" << endl;
}

void Simulation::cycle() {
    for (int i = 0; i < _M1; i++) {

        // STEP 2, 3
        vector<double> R;
        double Rmax = 0;

        cout << "R";
        for (int j = 0; j < reactions(); j++) {
            R.push_back(calcR(j));
            cout << " | " << R[j];
            Rmax = R[j] > Rmax ? R[j] : Rmax;
        }
        cout << endl;

        // STEP 4
        double random = RandomUniform();
        for (int j = 0; j < reactions(); j++) {
            if ((_A * R[j] / Rmax) >= random) {
                for (int k = 0; k < _reactions[j].reactants(); k++) {
                    int reactant = _reactions[j].reactant(k);
                    _N[reactant] -= _reactions[j].coefficient(reactant);
                }
                for (int k = 0; k < _reactions[j].products(); k++) {
                    int product = _reactions[j].product(k);
                    _N[product] += _reactions[j].coefficient(product);
                }
            }
        }

        // STEP 5 (repeat)
    }

    for (int i = 0; i < _species; i++) {
        _dNdt[i] = calcdNdt(i);
    }
}

double Simulation::calcR(int alpha) {
    double R = _reactions[alpha].k();
    for (int i = 0; i < _reactions[alpha].reactants(); i++) {
        int species = _reactions[alpha].reactant(i);
        R *= pow(_N[species], _reactions[alpha].coefficient(species));
    }
    return R;
}

double Simulation::getN(int species) {
    return _N[species];
}

void Simulation::printN() {
    cout << "N";
    for (auto i : _N) {
        cout << " | " << i;
    }
    cout << endl;
}

double Simulation::calcdNdt(int species) {
    double dNdt = 0;
    for (auto i : _reactions) {
        double sumTerm = 0;
        if (int id = i.exists(species)) {
            sumTerm = id > 0 ? 1 : -1;
            sumTerm *= i.coefficient(species) * i.k();
            for (int j = 0; j < i.reactants(); j++) {
                int reactant = i.reactant(j);
                if (reactant != species) {
                    sumTerm *= pow(_N[reactant], i.coefficient(reactant));
                }
            }
            dNdt += sumTerm;
        }
    }
//    cout << "calc dN_" << species << "/dt = " << dNdt << endl;
    return dNdt;
}