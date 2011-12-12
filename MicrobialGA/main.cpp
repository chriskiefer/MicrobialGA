//
//  main.cpp
//  MicrobialGA
//
//  Created by Chris Kiefer on 07/12/2011.
//  Copyright (c) 2011 Goldsmiths Creative Computing. All rights reserved.
//

#include <iostream>
#include "microbialGA.h"

using namespace std;

class gaTest : public objectiveFunctionEvaluator {
public:
    void go() {
        cout << "Microbial GA... testing\n";
        microbialGA ga(30, 15, 10, 0.5, 0.01, this, microbialGA::LOWSCOREISBEST);
//        ga.evolve(10000);
        ga.evolveUntil(0.001);
        genotype g = ga.getFittestIndividual();
        cout << "Winner: ";
        for(int i=0; i < g.size(); i++) cout << g[i] << ",";
        cout << endl;
    }
    float evaluate(genotype &a) {
        //return mean value
        float sum=0;
        vector<float> floatdata;
        microbialGA::genotypeToFloat(a, floatdata);
        for(int i=0; i < a.size(); i++) {
            sum += floatdata[i];
        }
        return sum / a.size();
    }
};

int main (int argc, const char * argv[])
{
    gaTest test;
    test.go();
    return 0;
}

