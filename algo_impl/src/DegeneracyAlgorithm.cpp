#include <limits.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <algorithm>

#include "Tools.h"
#include <list>
#include <vector>
#include "MemoryManager.h"
#include "DegeneracyTools.h"

#include "DegeneracyAlgorithm.h"
#include <sys/time.h>

static long long recursiveCallCount(0);

using namespace std;
DegeneracyAlgorithm::DegeneracyAlgorithm(vector<list<int>> const &adjacencyList)
 : Algorithm("degeneracy")
 , m_AdjacencyList(adjacencyList)
{
}

DegeneracyAlgorithm::~DegeneracyAlgorithm()
{
}

long DegeneracyAlgorithm::Run(list<list<int>> &cliques)
{
    return listAllMaximalCliquesDegeneracy(m_AdjacencyList, m_AdjacencyList.size());
}




static unsigned long largestDifference(0);
static unsigned long numLargeJumps;
static unsigned long stepsSinceLastReportedClique(0);
static unsigned long largest_mc_size(0);

long long DegeneracyAlgorithm::listAllMaximalCliquesDegeneracy(vector<list<int>> const &adjList, int size)
{
    // vertex sets are stored in an array like this:
    // |--X--|--P--|
    int* vertexSets = (int*)Calloc(size, sizeof(int));

    // vertex i is stored in vertexSets[vertexLookup[i]]
    int* vertexLookup = (int*)Calloc(size, sizeof(int));

    int** neighborsInP = (int**)Calloc(size, sizeof(int*));

    int* numNeighbors = (int*)Calloc(size, sizeof(int));
    
    NeighborListArray** orderingArray = computeDegeneracyOrderArray(adjList, size);

    int i = 0;

    while(i<size)
    {
        vertexLookup[i] = i;
        vertexSets[i] = i;
        neighborsInP[i] = (int*)Calloc(1, sizeof(int));
        numNeighbors[i] = 1;
        i++;
    }

    int beginX = 0;
    int beginP = 0;
    int beginR = size;

    long long cliqueCount = 0;

    list<int> partialClique;
    clock_t start = clock();

    for(i=0;i<size;i++)
    {
        
        int vertex = (int)orderingArray[i]->vertex;
        largest_mc_size = 0;


        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
        #endif

        // add vertex to partial clique R
        partialClique.push_back(vertex);

        int newBeginX, newBeginP, newBeginR;


        // set P to be later neighbors and X to be be earlier neighbors
        // of vertex
        fillInPandXForRecursiveCallDegeneracy( i, vertex, 
                                               vertexSets, vertexLookup, 
                                               orderingArray,
                                               neighborsInP, numNeighbors,
                                               &beginX, &beginP, &beginR, 
                                               &newBeginX, &newBeginP, &newBeginR);

        long long edges = 0;
        for (int idx = newBeginP; idx < newBeginR; ++idx) {
            int ver = vertexSets[idx];
            edges += numNeighbors[ver];
        }
        assert(edges%2==0);
        edges/=2;

        // recursively compute maximal cliques containing vertex, some of its
        // later neighbors, and avoiding earlier neighbors
        listAllMaximalCliquesDegeneracyRecursive( &cliqueCount,
                                                  partialClique, 
                                                  vertexSets, vertexLookup,
                                                  neighborsInP, numNeighbors,
                                                  newBeginX, newBeginP, newBeginR); 

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        beginR = beginR + 1;
        int sizeofP = orderingArray[i]->laterDegree;
        double dN_density = double(edges) / ((sizeofP-1)*sizeofP/2);
        double portion = double(largest_mc_size) / sizeofP;
        partialClique.pop_back();

    }

    cout << "Recursive Call Count: " << recursiveCallCount << endl;
    clock_t end = clock();
    cout << "real computation time: " << (double)(end-start)/(double)(CLOCKS_PER_SEC) << endl;

    partialClique.clear();

    Free(vertexSets);
    Free(vertexLookup);


    return cliqueCount;
}



void DegeneracyAlgorithm::listAllMaximalCliquesDegeneracyRecursive(long long* cliqueCount,
                                               list<int> &partialClique, 
                                               int* vertexSets, int* vertexLookup,
                                               int** neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR)
{
    recursiveCallCount++;

    stepsSinceLastReportedClique++;

    // if X is empty and P is empty, process partial clique as maximal
    if(beginX >= beginP && beginP >= beginR)
    {
        (*cliqueCount)++;

        if (stepsSinceLastReportedClique > partialClique.size()) {
            numLargeJumps++;
            if (largestDifference < (stepsSinceLastReportedClique - partialClique.size())) {
                largestDifference = stepsSinceLastReportedClique - partialClique.size();
            }
        }

        stepsSinceLastReportedClique = 0;

        //FIXME:
        largest_mc_size = max(largest_mc_size, partialClique.size()-1);

        ExecuteCallBacks(partialClique);
        processClique(partialClique);

        return;
    }

    // avoid work if P is empty.
    if(beginP >= beginR)
        return;

    int* myCandidatesToIterateThrough;
    int numCandidatesToIterateThrough;

    // get the candidates to add to R to make a maximal clique
    findBestPivotNonNeighborsDegeneracy( &myCandidatesToIterateThrough,
                                         &numCandidatesToIterateThrough,
                                         vertexSets, vertexLookup,
                                         neighborsInP, numNeighbors,
                                         beginX, beginP, beginR);

    // add candiate vertices to the partial clique one at a time and 
    // search for maximal cliques
    if(numCandidatesToIterateThrough != 0)
    {
    int iterator = 0;
    while(iterator < numCandidatesToIterateThrough)
    {
        // vertex to be added to the partial clique
        int vertex = myCandidatesToIterateThrough[iterator];

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("%d ", vertex);
        #endif

        int newBeginX, newBeginP, newBeginR;

        // add vertex into partialClique, representing R.
        partialClique.push_back(vertex);

        // swap vertex into R and update all data structures 
        moveToRDegeneracy( vertex, 
                           vertexSets, vertexLookup, 
                           neighborsInP, numNeighbors,
                           &beginX, &beginP, &beginR, 
                           &newBeginX, &newBeginP, &newBeginR);

        // recursively compute maximal cliques with new sets R, P and X
        listAllMaximalCliquesDegeneracyRecursive(cliqueCount,
                                                 partialClique, 
                                                 vertexSets, vertexLookup,
                                                 neighborsInP, numNeighbors,
                                                 newBeginX, newBeginP, newBeginR);

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        // remove vertex from partialClique
        partialClique.pop_back();

        moveFromRToXDegeneracy( vertex, 
                                vertexSets, vertexLookup,
                                &beginX, &beginP, &beginR );

        iterator++;
    }

    // swap vertices that were moved to X back into P, for higher recursive calls.
    iterator = 0;

    while(iterator < numCandidatesToIterateThrough)
    {
        int vertex = myCandidatesToIterateThrough[iterator];
        int vertexLocation = vertexLookup[vertex];

        beginP--;
        vertexSets[vertexLocation] = vertexSets[beginP];
        vertexSets[beginP] = vertex;
        vertexLookup[vertex] = beginP;
        vertexLookup[vertexSets[vertexLocation]] = vertexLocation;

        iterator++;
    }
    }

    // don't need to check for emptiness before freeing, since
    // something will always be there (we allocated enough memory
    // for all of P, which is nonempty)
    Free(myCandidatesToIterateThrough);
}
