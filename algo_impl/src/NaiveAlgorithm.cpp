
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

#include "NaiveAlgorithm.h"
#include <sys/time.h>

static long long recursiveCallCount(0);

using namespace std;
NaiveAlgorithm::NaiveAlgorithm(vector<list<int>> const &adjacencyList)
 : Algorithm("Naive")
 , m_AdjacencyList(adjacencyList)
{
}

NaiveAlgorithm::~NaiveAlgorithm()
{
}

long NaiveAlgorithm::Run(list<list<int>> &cliques)
{
    return listAllMaximalCliquesNaive(m_AdjacencyList, m_AdjacencyList.size());
}

inline int findBestPivotNonNeighborsNaive(  int** pivotNonNeighbors, int* numNonNeighbors,
                                            int* vertexSets, int* vertexLookup,
                                            int** neighborsInP, int* numNeighbors,
                                            int beginX, int beginP, int beginR,
                                            int *isClique)
{
    int pivot = -1;
    int maxIntersectionSize = -1;

    bool PisClique = true;
    int maxdegree = beginR - beginP - 1;

    // iterate over each vertex in P union X 
    // to find the vertex with the most neighbors in P.
    int j = beginX;
    while(j<beginR)
    {
        int vertex = vertexSets[j];
        int numPotentialNeighbors = min(beginR - beginP, numNeighbors[vertex]);

        int numNeighborsInP = 0;

        int k = 0;
        while(k<numPotentialNeighbors)
        {
            int neighbor = neighborsInP[vertex][k];
            int neighborLocation = vertexLookup[neighbor];

            if(neighborLocation >= beginP && neighborLocation < beginR)
            {
                numNeighborsInP++;
            }
            else
            {
                break;
            }

            k++;
        }

        if(numNeighborsInP > maxIntersectionSize)
        {
            pivot = vertex;
            maxIntersectionSize = numNeighborsInP;
        }

        if(j >= beginP && numNeighborsInP != maxdegree)
            PisClique = false;
        /*
        if(j < beginP && numNeighborsInP == maxdegree+1)
            PisClique = false;
        */

        j++;
    }

    // compute non neighbors of pivot by marking its neighbors
    // and moving non-marked vertices into pivotNonNeighbors.
    // we must do this because this is an efficient way
    // to compute non-neighbors of a vertex in 
    // an adjacency list.

    // we initialize enough space for all of P; this is
    // slightly space inefficient, but it results in faster
    // computation of non-neighbors.
    *pivotNonNeighbors = (int*)Calloc(beginR-beginP, sizeof(int));
    memcpy(*pivotNonNeighbors, &vertexSets[beginP], (beginR-beginP)*sizeof(int));

    // we will decrement numNonNeighbors as we find neighbors
    *numNonNeighbors = beginR-beginP;

    int numPivotNeighbors = min(beginR - beginP, numNeighbors[pivot]);

    // mark the neighbors of pivot that are in P.
    j = 0;
    while(j<numPivotNeighbors)
    {
        int neighbor = neighborsInP[pivot][j];
        int neighborLocation = vertexLookup[neighbor];

        if(neighborLocation >= beginP && neighborLocation < beginR)
        {
            (*pivotNonNeighbors)[neighborLocation-beginP] = -1;
        }
        else
        {
            break;
        }

        j++;
    }

    // move non-neighbors of pivot in P to the beginning of
    // pivotNonNeighbors and set numNonNeighbors appriopriately.

    // if a vertex is marked as a neighbor, the we move it
    // to the end of pivotNonNeighbors and decrement numNonNeighbors.
    j = 0;
    while(j<*numNonNeighbors)
    {
        int vertex = (*pivotNonNeighbors)[j];

        if(vertex == -1)
        {
            (*numNonNeighbors)--;
            (*pivotNonNeighbors)[j] = (*pivotNonNeighbors)[*numNonNeighbors];
            continue;
        }

        j++;
    }
    *isClique = PisClique ? 1 : 0;

    return pivot; 
}

inline long long fillInPandXForRecursiveCallNaive( int vertex, int orderNumber,
                              int* vertexSets, int* vertexLookup, 
                              NeighborListArray** orderingArray,
                              int** neighborsInP, int* numNeighbors,
                              int* pBeginX, int *pBeginP, int *pBeginR, 
                              int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{
        long long edges = 0;
        int vertexLocation = vertexLookup[vertex];

        (*pBeginR)--;
        vertexSets[vertexLocation] = vertexSets[*pBeginR];
        vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
        vertexSets[*pBeginR] = vertex;
        vertexLookup[vertex] = *pBeginR;

        *pNewBeginR = *pBeginR;
        *pNewBeginP = *pBeginR;

        // swap later neighbors of vertex into P section of vertexSets
        int j = 0;
        while(j<orderingArray[orderNumber]->laterDegree)
        {
            int neighbor = orderingArray[orderNumber]->later[j];
            int neighborLocation = vertexLookup[neighbor];

            (*pNewBeginP)--;

            vertexSets[neighborLocation] = vertexSets[*pNewBeginP];
            vertexLookup[vertexSets[*pNewBeginP]] = neighborLocation;
            vertexSets[*pNewBeginP] = neighbor;
            vertexLookup[neighbor] = *pNewBeginP;

            j++; 
        }

        *pNewBeginX = *pNewBeginP;

        // swap earlier neighbors of vertex into X section of vertexSets
        j = 0;
        while(j<orderingArray[orderNumber]->earlierDegree)
        {
            int neighbor = orderingArray[orderNumber]->earlier[j];
            int neighborLocation = vertexLookup[neighbor];

            (*pNewBeginX)--;
            vertexSets[neighborLocation] = vertexSets[*pNewBeginX];
            vertexLookup[vertexSets[*pNewBeginX]] = neighborLocation;
            vertexSets[*pNewBeginX] = neighbor;
            vertexLookup[neighbor] = *pNewBeginX;

            Free(neighborsInP[neighbor]);
            neighborsInP[neighbor] = (int*)Calloc(min(*pNewBeginR-*pNewBeginP,orderingArray[neighbor]->laterDegree), sizeof(int));
            numNeighbors[neighbor] = 0;

            // fill in NeighborsInP
            int k = 0;
            while(k<orderingArray[neighbor]->laterDegree)
            {
                int laterNeighbor = orderingArray[neighbor]->later[k];
                int laterNeighborLocation = vertexLookup[laterNeighbor];
                if(laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR)
                {
                    neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                    numNeighbors[neighbor]++;
                }

                k++;
            }

            j++; 

        }

        // reset numNeighbors and neighborsInP for this vertex
        j = *pNewBeginP;
        while(j<*pNewBeginR)
        {
            int vertexInP = vertexSets[j];
            numNeighbors[vertexInP] = 0;
            Free(neighborsInP[vertexInP]);
            neighborsInP[vertexInP]=(int*)Calloc( min( *pNewBeginR-*pNewBeginP, 
                                                 orderingArray[vertexInP]->laterDegree 
                                              + orderingArray[vertexInP]->earlierDegree), sizeof(int));

            j++;
        }

        // count neighbors in P, and fill in array of neighbors
        // in P
        j = *pNewBeginP;
        while(j<*pNewBeginR)
        {
            int vertexInP = vertexSets[j];

            int k = 0;
            while(k<orderingArray[vertexInP]->laterDegree)
            {
                int laterNeighbor = orderingArray[vertexInP]->later[k];
                int laterNeighborLocation = vertexLookup[laterNeighbor];

                if(laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR)
                {
                    neighborsInP[vertexInP][numNeighbors[vertexInP]] = laterNeighbor;
                    numNeighbors[vertexInP]++;
                    neighborsInP[laterNeighbor][numNeighbors[laterNeighbor]] = vertexInP;
                    numNeighbors[laterNeighbor]++;
                    edges++;
                }

                k++;
            }

            j++;
        }
        return edges;
}

static unsigned long largestDifference(0);
static unsigned long numLargeJumps;
static unsigned long stepsSinceLastReportedClique(0);
static unsigned long largest_mc_size(0);


long long NaiveAlgorithm::listAllMaximalCliquesNaive(vector<list<int>> const &adjList, int size)
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
    clock_t start =clock();

    // for each vertex
    // Added by yinuo for counting maximum time

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
        long long edges = fillInPandXForRecursiveCallNaive( i, vertex, 
                                               vertexSets, vertexLookup, 
                                               orderingArray,
                                               neighborsInP, numNeighbors,
                                               &beginX, &beginP, &beginR, 
                                               &newBeginX, &newBeginP, &newBeginR);

        for (int j = newBeginR; j < size; ++j)
            partialClique.push_back(vertexSets[j]);

        long sizeofP = newBeginR - newBeginP;

        double density = 2.0*double(edges) / ((sizeofP-1)*sizeofP);
        if (sizeofP == 0 || sizeofP == 1)
            density = 1.0;

        if (density == 1.0) {
            bool commonNeighborInX = false;
            for (int j = newBeginX; j < newBeginP; ++j) {
                 int vertexInX = vertexSets[j];
                 int numNeighborsInP = 0;
                 for (int k = 0; k < numNeighbors[vertexInX]; ++k) {
                     int neighborLocation = vertexLookup[neighborsInP[vertexInX][k]];
                     if (neighborLocation >= newBeginP  && neighborLocation < newBeginR)
                         numNeighborsInP++;
                 }
                 if (numNeighborsInP == sizeofP) {
                     commonNeighborInX = true;
                     break;
                 }
            }
            if (!commonNeighborInX) {
                for (int k = newBeginP; k < newBeginR; ++k)
                    partialClique.push_back(vertexSets[k]);
                ExecuteCallBacks(partialClique);
                processClique(partialClique);
                cliqueCount++;
                for (int k = newBeginP; k < newBeginR; ++k)
                    partialClique.pop_back();
            }
        } else {
            // recursively compute maximal cliques containing vertex, some of its
            // later neighbors, and avoiding earlier neighbors
            listAllMaximalCliquesNaiveRecursive( &cliqueCount,
                                                 partialClique, 
                                                 vertexSets, vertexLookup,
                                                 neighborsInP, numNeighbors,
                                                 newBeginX, newBeginP, newBeginR); 
        }

        beginR = beginR + 1;
        for (int j = newBeginR; j < size; ++j)
            partialClique.pop_back();
        partialClique.pop_back();

    }

    cout << "Recursive Call Count: " << recursiveCallCount << endl;
    clock_t end = clock();
    cout << "real computation time: " << (double)(end-start)/(double)(CLOCKS_PER_SEC) << endl;

    partialClique.clear();

    Free(vertexSets);
    Free(vertexLookup);

    for(i = 0; i<size; i++)
    {
        delete orderingArray[i];
    }

    Free(orderingArray);
    Free(numNeighbors);

    return cliqueCount;
}

inline void moveToRNaive( int vertex, 
                          int* vertexSets, int* vertexLookup, 
                          int** neighborsInP, int* numNeighbors,
                          int* pBeginX, int *pBeginP, int *pBeginR, 
                          int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{

        int vertexLocation = vertexLookup[vertex];

        (*pBeginR)--;
        vertexSets[vertexLocation] = vertexSets[*pBeginR];
        vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
        vertexSets[*pBeginR] = vertex;
        vertexLookup[vertex] = *pBeginR;

        // this is not a typo, initially newX is empty
        *pNewBeginX = *pBeginP;
        *pNewBeginP = *pBeginP;
        *pNewBeginR = *pBeginP;

        int sizeOfP = *pBeginR - *pBeginP;

        int j = *pBeginX;
        while(j<*pNewBeginX)
        {
            int neighbor = vertexSets[j];
            int neighborLocation = j;

            int incrementJ = 1;

            //int numPotentialNeighbors = min(sizeOfP, numNeighbors[neighbor]);
            int numPotentialNeighbors = numNeighbors[neighbor];

            int k = 0;
            while(k<numPotentialNeighbors)
            {
                if(neighborsInP[neighbor][k] == vertex)
                {
                    (*pNewBeginX)--;
                    vertexSets[neighborLocation] = vertexSets[(*pNewBeginX)];
                    vertexLookup[vertexSets[(*pNewBeginX)]] = neighborLocation;
                    vertexSets[(*pNewBeginX)] = neighbor;
                    vertexLookup[neighbor] = (*pNewBeginX);
                    incrementJ=0;
                }
                
                k++;
            }

            if(incrementJ) j++;
        }

        j = (*pBeginP);
        while(j<(*pBeginR))
        {
            int neighbor = vertexSets[j];
            int neighborLocation = j;

            //int numPotentialNeighbors = min(sizeOfP, numNeighbors[neighbor]);
            int numPotentialNeighbors =  numNeighbors[neighbor];

            int k = 0;
            while(k<numPotentialNeighbors)
            {
                if(neighborsInP[neighbor][k] == vertex)
                {
                    vertexSets[neighborLocation] = vertexSets[(*pNewBeginR)];
                    vertexLookup[vertexSets[(*pNewBeginR)]] = neighborLocation;
                    vertexSets[(*pNewBeginR)] = neighbor;
                    vertexLookup[neighbor] = (*pNewBeginR);
                    (*pNewBeginR)++;
                }

                k++;
            }

            j++;
        }

        j = (*pNewBeginX);

        while(j < *pNewBeginR)
        {
            int thisVertex = vertexSets[j];

            //int numPotentialNeighbors = min(sizeOfP, numNeighbors[thisVertex]);
            int numPotentialNeighbors = numNeighbors[thisVertex];

            int numNeighborsInP = 0;

            int k = 0;
            while(k < numPotentialNeighbors)
            {
                int neighbor = neighborsInP[thisVertex][k];
                int neighborLocation = vertexLookup[neighbor];
                if(neighborLocation >= *pNewBeginP && neighborLocation < *pNewBeginR)
                {
                    neighborsInP[thisVertex][k] = neighborsInP[thisVertex][numNeighborsInP];
                    neighborsInP[thisVertex][numNeighborsInP] = neighbor;
                    numNeighborsInP++;
                }
                k++;
            }

            j++;
        }

}


inline void moveFromRToXNaive( int vertex, 
                                    int* vertexSets, int* vertexLookup, 
                                    int* pBeginX, int* pBeginP, int* pBeginR )
{
    int vertexLocation = vertexLookup[vertex];

    //swap vertex into X and increment beginP and beginR
    vertexSets[vertexLocation] = vertexSets[*pBeginP];
    vertexLookup[vertexSets[*pBeginP]] = vertexLocation;
    vertexSets[*pBeginP] = vertex;
    vertexLookup[vertex] = *pBeginP;

    *pBeginP = *pBeginP + 1;
    *pBeginR = *pBeginR + 1;

}


void NaiveAlgorithm::listAllMaximalCliquesNaiveRecursive(long long* cliqueCount,
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

    int pisclique = 0;

    // get the candidates to add to R to make a maximal clique
    // pisclique 会返回当前的 P 集合是否已经是一个 Clique
    findBestPivotNonNeighborsNaive( &myCandidatesToIterateThrough,
                                    &numCandidatesToIterateThrough,
                                    vertexSets, vertexLookup,
                                    neighborsInP, numNeighbors,
                                    beginX, beginP, beginR, &pisclique);

    if (pisclique) {
        for (int j = beginP; j < beginR; ++j)
            partialClique.push_back(vertexSets[j]);
        ExecuteCallBacks(partialClique);
        processClique(partialClique);
        (*cliqueCount)++;
        
        for (int j = beginP; j < beginR; ++j)
            partialClique.pop_back();

        return ;
    }

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
        moveToRNaive( vertex, 
                       vertexSets, vertexLookup, 
                       neighborsInP, numNeighbors,
                       &beginX, &beginP, &beginR, 
                       &newBeginX, &newBeginP, &newBeginR);

        // recursively compute maximal cliques with new sets R, P and X
        listAllMaximalCliquesNaiveRecursive(cliqueCount,
                                                 partialClique, 
                                                 vertexSets, vertexLookup,
                                                 neighborsInP, numNeighbors,
                                                 newBeginX, newBeginP, newBeginR);

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        // remove vertex from partialClique
        partialClique.pop_back();

        moveFromRToXNaive( vertex, 
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
