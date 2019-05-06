
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

#include "newHybrid.h"
#include "BKReverse.h"
#include "DegeneracyAlgorithm.h"
#include <sys/time.h>

//#define RB_DEBUG

static long long recursiveCallCount(0);

using namespace std;
newHybrid::newHybrid(vector<list<int>> const &adjacencyList, double _density_threshold)
 : Algorithm("hybrid")
 , m_AdjacencyList(adjacencyList)
 , density_threshold(_density_threshold)
{}

newHybrid::~newHybrid()
{}

long newHybrid::Run(list<list<int>> &cliques)
{
    return newHybridMain(m_AdjacencyList, m_AdjacencyList.size());
}

inline long long fillInPandX( int vertex, int orderNumber,
                              int* vertexSets, int* vertexLookup, 
                              NeighborListArray** orderingArray,
                              vector<vector<int>>& neighborsInP, int* numNeighbors,
                              int* pBeginX, int *pBeginP, int *pBeginR, 
                              int* pNewBeginX, int* pNewBeginP, int *pNewBeginR,
                              long& scount, long& kcount, long& noise, long& dive)
{
        int vertexLocation = vertexLookup[vertex];

        (*pBeginR)--;
        vertexSets[vertexLocation] = vertexSets[*pBeginR];
        vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
        vertexSets[*pBeginR] = vertex;
        vertexLookup[vertex] = *pBeginR;

        *pNewBeginR = *pBeginR;
        *pNewBeginP = *pBeginR;

        // |-------------------|vertex|
        //                     ^
        //                 pNewBeginP

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

        // reset numNeighbors and neighborsInP for this vertex
        j = *pNewBeginP;
        while(j<*pNewBeginR)
        {
            int vertexInP = vertexSets[j];
            numNeighbors[vertexInP] = 0;

            neighborsInP[vertexInP].resize(min( *pNewBeginR-*pNewBeginP, orderingArray[vertexInP]->laterDegree + orderingArray[vertexInP]->earlierDegree));

            j++;
        }

        // count neighbors in P, and fill in array of neighbors
        // in P
        j = *pNewBeginP;
        long long edges = 0;
        while(j<*pNewBeginR)
        {
            int vertexInP = vertexSets[j];

            int k = 0;
            while(k < orderingArray[vertexInP]->laterDegree)
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

        int oldBeginR = *pNewBeginR;
        int sizeofP = *pNewBeginR - *pNewBeginP;
        int vers_max_degree = 0;
        int max_degree = sizeofP - 1;

        scount = 0;
        kcount = 0;
        noise = 0;
        dive = 0;

        j = *pNewBeginP;
        while(j < *pNewBeginR)
        {
            int vertexInP = vertexSets[j];

            if (numNeighbors[vertexInP] == max_degree) {
                scount++;
                vers_max_degree += 1;
                int vertexLocation = vertexLookup[vertexInP];

                (*pNewBeginR)--;
                
                vertexSets[vertexLocation] = vertexSets[*pNewBeginR];
                vertexLookup[vertexSets[*pNewBeginR]] = vertexLocation;
                vertexSets[*pNewBeginR] = vertexInP;
                vertexLookup[vertexInP] = *pNewBeginR;


                continue;
            }

            ++j;
        }

        sizeofP = *pNewBeginR - *pNewBeginP;
        kcount = sizeofP;

        j = *pNewBeginP;
        while(j < *pNewBeginR) 
        {
            int vertexInP = vertexSets[j];
            int k = 0;
            while (k < numNeighbors[vertexInP]) 
            {
                if ( vertexLookup[neighborsInP[vertexInP][k]] >= *pNewBeginR ) 
                {
                    swap(neighborsInP[vertexInP][k], neighborsInP[vertexInP][--numNeighbors[vertexInP]]);
                    continue;
                }
                ++k;
            }

            noise += numNeighbors[vertexInP];
            dive = max((int)dive, numNeighbors[vertexInP]);

            ++j;
        }

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

            neighborsInP[neighbor].resize(min(*pNewBeginR-*pNewBeginP,orderingArray[neighbor]->laterDegree));
            numNeighbors[neighbor] = 0;

            // fill in NeighborsInP
            int k = 0;
            int degreeInR = 0;
            while(k<orderingArray[neighbor]->laterDegree)
            {
                int laterNeighbor = orderingArray[neighbor]->later[k];
                int laterNeighborLocation = vertexLookup[laterNeighbor];
                if(laterNeighborLocation >= *pNewBeginP && laterNeighborLocation < *pNewBeginR)
                {
                    neighborsInP[neighbor][numNeighbors[neighbor]] = laterNeighbor;
                    numNeighbors[neighbor]++;
                } else if (laterNeighborLocation >= *pNewBeginR && laterNeighborLocation < oldBeginR) {
                    degreeInR += 1;
                }

                k++;
            }

            if ( degreeInR < vers_max_degree ) {
                (*pNewBeginX)++;
            }


            j++; 

        }

        return edges;
}

inline int findBestPivotNonNeighborsDgcy( int** pivotNonNeighbors, int* numNonNeighbors,
                                                int* vertexSets, int* vertexLookup,
                                                vector<vector<int>>& neighborsInP, int* numNeighbors,
                                                int beginX, int beginP, int beginR)
{
    int pivot = -1;
    int maxIntersectionSize = -1;

    // iterate over each vertex in P union X 
    // to find the vertex with the most neighbors in P.
    int j = beginX;
    while(j<beginR)
    {
        int vertex = vertexSets[j];
        //int numPotentialNeighbors = min(beginR - beginP, numNeighbors[vertex]);
        int numPotentialNeighbors = numNeighbors[vertex];

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

    //int numPivotNeighbors = min(beginR - beginP, numNeighbors[pivot]);
    int numPivotNeighbors = numNeighbors[pivot];

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
    return pivot; 
}

inline void moveToRDgcy( int vertex, 
                         int* vertexSets, int* vertexLookup, 
                         vector<vector<int>>& neighborsInP, int* numNeighbors,
                         int* pBeginX, int *pBeginP, int *pBeginR, 
                         int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{

////    clock_t clockStart = clock();
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

/*! \brief Move a vertex from the set R to the set X, and update all necessary pointers
           and arrays of neighbors in P

    \param vertex The vertex to move from R to X.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param pBeginX The index where set X begins in vertexSets.
 
    \param pBeginP The index where set P begins in vertexSets.

    \param pBeginR The index where set R begins in vertexSets.

*/

inline void moveFromRToXDgcy( int vertex, 
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

static int _size;

long long newHybrid::newHybridMain(vector<list<int>> const &adjList, int size)
{
    // vertex sets are stored in an array like this:
    // vertex i is stored in vertexSets[vertexLookup[i]]
    // |--X--|--P--|
    int*  vertexSets   = (int*) Calloc(size, sizeof(int));
    int*  vertexLookup = (int*) Calloc(size, sizeof(int));
    vector<vector<int>> neighborsInP(size, vector<int>(0));
    int*  numNeighbors = (int*) Calloc(size, sizeof(int));

    // compute the degeneracy order
    NeighborListArray** orderingArray = computeDegeneracyOrderArray(adjList, size);

    int i = 0;

    while(i<size)
    {
        vertexLookup[i] = i;
        vertexSets[i]   = i;
        //neighborsInP_[i] = (int*)Calloc(1, sizeof(int));
        numNeighbors[i] = 1;
        i++;
    }

    int beginX = 0;
    int beginP = 0;
    int beginR = size;

    long long cliqueCount = 0, precliqueCount = 0;
    long long BKreverseCount = 0, BKpivotCount = 0;

    // scount 是全联通的点的个数, kcount是其余点的个数, noise是kcount所有点之间的边数
    long scount = 0, kcount = 0, noise = 0, dive = 0;

    list<int> partialClique;
    clock_t start =clock();

    for (int i = 0; i < size; ++i)
    {

        int vertex = (int)orderingArray[i]->vertex;

        // add vertex to partial clique R
        int newBeginX, newBeginP, newBeginR;
        int newBeginX_, newBeginP_, newBeginR_;


        // set P to be later neighbors and X to be be earlier neighbors of vertex
        long long edges = fillInPandX( i, vertex, 
                                       vertexSets, vertexLookup, 
                                       orderingArray,
                                       neighborsInP, numNeighbors,
                                       &beginX, &beginP, &beginR, 
                                       &newBeginX, &newBeginP, &newBeginR,
                                       scount, kcount, noise, dive);

        // 将新加入R集合中的点都加入 partialClique 中(因为其中含有全联通的点)
        for (int j = newBeginR; j < size; ++j)
            partialClique.push_back(vertexSets[j]);

        long long sizeofP = scount+kcount;

        if (sizeofP == 0) {
            if (newBeginX == newBeginP) {
                cliqueCount++;
                ExecuteCallBacks(partialClique);
                processClique(partialClique);
            }
            beginR = beginR + 1;
            for (int j = newBeginR; j < size; ++j)
                partialClique.pop_back();
            continue;
        }

        /*
        if (sizeofP == 1){
            int j;
            for (j = newBeginX; j < newBeginP; ++j)
                if (numNeighbors[vertexSets[j]] == 1)
                    break;
            if (j == newBeginP) {
                cliqueCount++;
                ExecuteCallBacks(partialClique);
                processClique(partialClique);
            }
            beginR = beginR + 1;
            for (int j = newBeginR; j < size; ++j)
                partialClique.pop_back();
            continue;
        }
        */


        // 是否调用 BKrcd 的判断条件
        //double benefit = (scount*scount*scount)/3.0+(1-kcount)*scount*scount+(-kcount*kcount+kcount*2-1)*scount+kcount;
        //double benefit = (scount*scount*scount)/2.0+(2-kcount/2.0)*scount*scount+(-2.5*kcount*kcount+kcount*2+3*kcount-3.5)*scount+2;
        //double benefit = double(scount)/(scount+kcount)-0.76;
        //double benefit = scount+3.0-2.7*kcount;

        double benefit = 0.0;
        //if (scount == 0) {
        //    benefit = kcount;
        //    //benefit = scount+3.0-2.7*kcount;
        //} else {
            if (dive == 0) {
                benefit = scount+4.5-2.8*kcount;
            } else if (dive == 1) {
                benefit = scount+8.0-2.8*kcount;
            } else if (dive >= 2) {
                benefit = scount+11.0-2.8*kcount;
            }
        //}

        if (benefit >= 0.0) {
            BKreverseCount++;
            BKReverseRecursive( &cliqueCount, partialClique,
                                vertexSets, vertexLookup,
                                neighborsInP, numNeighbors,
                                newBeginX, newBeginP, newBeginR); 
        } else {
            listAllMaximalCliquesDgcyRecursive( &cliqueCount,
                                              partialClique, 
                                              vertexSets, vertexLookup,
                                              neighborsInP, numNeighbors,
                                              newBeginX, newBeginP, newBeginR); 
            BKpivotCount++;
        }

        beginR = beginR + 1;
        for (int j = newBeginR; j < size; ++j)
            partialClique.pop_back();

    }

    cout << endl << "BKreverseCount BKpivotCount BKreversePortion BKpivotPortion" << endl;
    cout << BKreverseCount << " " << BKpivotCount << " ";
    cout << double(BKreverseCount)/(BKreverseCount+BKpivotCount) << " " << double(BKpivotCount) / (BKreverseCount+BKpivotCount) << endl;
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

void newHybrid::BKReverseRecursive(long long* cliqueCount,
                                list<int> &partialClique, 
                                int* vertexSets, int* vertexLookup,
                                vector<vector<int>>& neighborsInP, int* numNeighbors,
                                int beginX, int beginP, int beginR)
{
    // if X is empty and P is empty, process partial clique as maximal
    if(beginX >= beginP && beginP >= beginR)
    {
        (*cliqueCount)++;

        ExecuteCallBacks(partialClique);
        processClique(partialClique);

        return;
    }

    // avoid work if P is empty.
    if(beginP >= beginR)
        return;

    // add candiate vertices to the partial clique one at a time and 
    // search for maximal cliques
    list<int> pivots;
    while(true)
    {
        int pivot = 0;
        if (findMostNonNeighbors ( vertexSets, vertexLookup, neighborsInP, numNeighbors, beginX, beginP, beginR, &pivot))
        {
            if (!existCommonNeighborOfPinX(vertexSets, vertexLookup, neighborsInP, numNeighbors, beginX, beginP, beginR)) 
            {
                //R U P is a maximal clique now
                for (int j = beginP; j < beginR; ++j)
                    partialClique.push_back(vertexSets[j]);
                (*cliqueCount)++;

                ExecuteCallBacks(partialClique);
                processClique(partialClique);

                for (int j = beginP; j < beginR; ++j)
                    partialClique.pop_back();
            }

            // swap vertices that were moved to X back into P, for higher recursive calls.
            for (auto p : pivots)
            {
                int vertexLocation = vertexLookup[p];
                swap(vertexSets[vertexLocation], vertexSets[--beginP]);
                swap(vertexLookup[p], vertexLookup[vertexSets[vertexLocation]]);
            }
            return ;
        }
        pivots.push_back(pivot);

        int newBeginX, newBeginP, newBeginR;

        // add vertex into partialClique, representing R.
        partialClique.push_back(pivot);

        // swaps pivot into R and update all data structures 
        moveToR( pivot, 
                 vertexSets, vertexLookup, 
                 neighborsInP, numNeighbors,
                 &beginX, &beginP, &beginR, 
                 &newBeginX, &newBeginP, &newBeginR);

        // recursively compute maximal cliques with new sets R, P and X
        BKReverseRecursive(cliqueCount, 
                           partialClique, 
                           vertexSets, vertexLookup,
                           neighborsInP, numNeighbors,
                           newBeginX, newBeginP, newBeginR);


        // remove vertex from partialClique
        partialClique.pop_back();

        moveFromRToX( pivot, 
                      vertexSets, vertexLookup,
                      neighborsInP, numNeighbors,
                      &beginX, &beginP, &beginR );
        int sizeOfP = beginR-beginP;
    }

}

void newHybrid::listAllMaximalCliquesDgcyRecursive(long long* cliqueCount,
                                               list<int> &partialClique, 
                                               int* vertexSets, int* vertexLookup,
                                               vector<vector<int>>& neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR)
{
    recursiveCallCount++;


    // if X is empty and P is empty, process partial clique as maximal
    if(beginX >= beginP && beginP >= beginR)
    {
        (*cliqueCount)++;
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
    findBestPivotNonNeighborsDgcy( &myCandidatesToIterateThrough,
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
        moveToRDgcy( vertex, 
                     vertexSets, vertexLookup, 
                     neighborsInP, numNeighbors,
                     &beginX, &beginP, &beginR, 
                     &newBeginX, &newBeginP, &newBeginR);

        // recursively compute maximal cliques with new sets R, P and X
        listAllMaximalCliquesDgcyRecursive(cliqueCount,
                                           partialClique, 
                                           vertexSets, vertexLookup,
                                           neighborsInP, numNeighbors,
                                           newBeginX, newBeginP, newBeginR);

        #ifdef PRINT_CLIQUES_TOMITA_STYLE
        printf("b ");
        #endif

        // remove vertex from partialClique
        partialClique.pop_back();

        moveFromRToXDgcy( vertex, 
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

