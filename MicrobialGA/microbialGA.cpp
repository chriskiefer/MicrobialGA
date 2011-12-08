//
//  microbialGA.cpp
//  MicrobialGA
//
//  Created by Chris Kiefer on 07/12/2011.
//  Copyright (c) 2011 Goldsmiths Creative Computing. All rights reserved.
//

#include "microbialGA.h"
#include <iostream>
#include <cmath>
#include <bitset>

inline float randUF() {
    return (float)rand() / (float) RAND_MAX;
}

microbialGA::microbialGA(unsigned int _populationSize, unsigned int _demeSize, unsigned int _geneSize, float _recombinationRate, float _mutationRate, objectiveFunctionEvaluator *_eval) : 
    populationSize(_populationSize), demeSize(_demeSize), geneSize(_geneSize), 
    recombinationRate(_recombinationRate), mutationRate(_mutationRate), evaluator(_eval)
{
    srand((int)time(NULL));
    
    //initialise the population randomly
    population.resize(populationSize);
    for(int i=0; i < populationSize; i++) {
        population[i].resize(geneSize);
        for(int j=0; j < geneSize; j++) {
            population[i][j] = randUF() * numeric_limits<unsigned int>::max();
        }
    }
    
    geneCharSize = geneSize * sizeof(unsigned int);
}

void microbialGA::evolve(unsigned int numIterations) {
    smoothedFitness = 0;
    for(int i=0; i < numIterations; i++) {
        microbialTournament();
        cout << "Iteration " << i << ", fitness: " << smoothedFitness << endl;
    }
}

void microbialGA::evolveUntil(float threshold) {
    smoothedFitness = 0;
    bestFitness = 0;
    int i=0;
    while(bestFitness < threshold) {
        microbialTournament();
        cout << "Iteration " << i << ", average fitness: " << smoothedFitness << ", best fitness: " << bestFitness << endl;
        i++;
    }    
}

genotype& microbialGA::getFittestIndividual() {
//    float maxFitness=0;
//    int maxIdx=0;
//    for(int i=0; i < populationSize; i++) {
//        float fitness = evaluator->evaluate(population[i]);
//        if (fitness > maxFitness) {
//            maxFitness = fitness;
//            maxIdx =i;
//        }
//    }
    return population[bestFitnessIndex];
}



unsigned char microbialGA::bitwiseInfectAndMutate(const unsigned char a, const unsigned char b) {
    bitset<8> bitsA(a);
    bitset<8> bitsB(b);
    for(int i=0; i < 8; i++) {
        if (randUF() < recombinationRate)
            bitsB[i] = bitsA[i];
        if (randUF() < mutationRate) 
            bitsB[i].flip();
    }
    return (char)bitsB.to_ulong();
}


void microbialGA::microbialTournament() {
    unsigned int contestant1, contestant2, winner, loser;
    contestant1 = floor(populationSize * randUF() * 0.99999);
    contestant2 = ((unsigned int)(contestant1 + 1 + floor(demeSize * randUF()))) % populationSize;
    float contestant1Fitness = evaluator->evaluate(population[contestant1]);
    float contestant2Fitness = evaluator->evaluate(population[contestant2]);
    float winningFitness;
    if (contestant1Fitness >  contestant2Fitness) {
        winner = contestant1;
        loser = contestant2;
        winningFitness = contestant1Fitness;
    }else{
        winner = contestant2;
        loser = contestant1;        
        winningFitness = contestant2Fitness;
    }
    
    smoothedFitness = (0.9 * smoothedFitness) + (0.1 * winningFitness); //exp moving average
    if (winningFitness > bestFitness) {
        bestFitness = winningFitness;
        bestFitnessIndex = winner;
    }
    
    unsigned char *winnerGene = (unsigned char*) &population[winner][0];
    unsigned char *loserGene = (unsigned char*) &population[loser][0];
    
    for(int i=0; i < geneCharSize; i++) {
        loserGene[i] = bitwiseInfectAndMutate(winnerGene[i], loserGene[i]);
    }
}



