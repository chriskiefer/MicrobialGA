//
//  microbialGA.h
//  MicrobialGA
//
//  Created by Chris Kiefer on 07/12/2011.
//  Copyright (c) 2011 Goldsmiths, University of London, EAVI. All rights reserved.
//

//based on "The Microbial Genetic Algorithm" by Inman Harvey  
//http://www.sussex.ac.uk/Users/inmanh/MicrobialGA_ECAL2009.pdf

#ifndef MicrobialGA_microbialGA_h
#define MicrobialGA_microbialGA_h
#include <vector>

using namespace std;

typedef vector<unsigned int> genotype;

class objectiveFunctionEvaluator {
public:
    virtual float evaluate(genotype &a) = 0;
};

class microbialGA {
public:
    microbialGA(unsigned int populationSize, unsigned int demeSize, unsigned int geneSize, float recombinationRate, float mutationRate,objectiveFunctionEvaluator *evaluator);
    void evolve(unsigned int numIterations);
    void evolveUntil(float threshold);
    genotype& getFittestIndividual();
private:
    unsigned int populationSize, demeSize, geneSize;
    float recombinationRate, mutationRate;
    vector<genotype> population;
    objectiveFunctionEvaluator *evaluator;
    unsigned int geneCharSize;
    
    void microbialTournament();
    inline unsigned char bitwiseInfectAndMutate(const unsigned char a, const unsigned char b);
    
    float smoothedFitness, bestFitness;
    unsigned int bestFitnessIndex;
};

#endif
