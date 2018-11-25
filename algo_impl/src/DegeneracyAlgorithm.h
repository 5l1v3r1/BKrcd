#ifndef _DJS_DEGENERACY_ALGORITHM_H_
#define _DJS_DEGENERACY_ALGORITHM_H_

/* 
    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by 
    the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version. 
 
    This program is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
    GNU General Public License for more details. 
 
    You should have received a copy of the GNU General Public License 
    along with this program.  If not, see <http://www.gnu.org/licenses/> 
*/

// local includes
#include "Algorithm.h"
#include "Tools.h"
#include "MemoryManager.h"
#include "DegeneracyTools.h"

// system includes
#include <list>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>

using namespace std;

/*! \file DegeneracyAlgorithm.h

    \brief see DegeneracyAlgorithm.cpp

    \author Darren Strash (first name DOT last name AT gmail DOT com)

    \copyright Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    \image html gplv3-127x51.png

    \htmlonly
    <center>
    <a href="gpl-3.0-standalone.html">See GPL 3.0 here</a>
    </center>
    \endhtmlonly
*/

class DegeneracyAlgorithm : public Algorithm
{
public:
    DegeneracyAlgorithm(std::vector<std::list<int>> const &adjacencyList);
    virtual ~DegeneracyAlgorithm();

    virtual long Run(std::list<std::list<int>> &cliques);

    DegeneracyAlgorithm           (DegeneracyAlgorithm const &) = delete;
    DegeneracyAlgorithm& operator=(DegeneracyAlgorithm const &) = delete;

    void listAllMaximalCliquesDegeneracyRecursive(long long* cliqueCount,
                                               std::list<int> &partialClique, 
                                               int* vertexSets, int* vertexLookup,
                                               int** neighborsInP, int* numNeighbors,
                                               //vector<vector<int>>& neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR);

    long long listAllMaximalCliquesDegeneracy(std::vector<std::list<int>> const &adjList, int size);

private:
    std::vector<std::list<int>> const &m_AdjacencyList;
};

inline int findBestPivotNonNeighborsDegeneracy( int** pivotNonNeighbors, int* numNonNeighbors,
                                                int* vertexSets, int* vertexLookup,
                                                int** neighborsInP, int* numNeighbors,
                                                //vector<vector<int>>& neighborsInP, int* numNeighbors,
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

inline long long  fillInPandXForRecursiveCallDegeneracy( int vertex, int orderNumber,
                                                   int* vertexSets, int* vertexLookup, 
                                                   NeighborListArray** orderingArray,
                                                   int** neighborsInP, int* numNeighbors,
                                                   //vector<vector<int>>& neighborsInP, int* numNeighbors,
                                                   int* pBeginX, int *pBeginP, int *pBeginR, 
                                                   int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{
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
            //neighborsInP[neighbor].resize(min(*pNewBeginR-*pNewBeginP,orderingArray[neighbor]->laterDegree));
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
            //neighborsInP[vertexInP].resize(min( *pNewBeginR-*pNewBeginP, 
            //                                     orderingArray[vertexInP]->laterDegree 
            //                                   + orderingArray[vertexInP]->earlierDegree));

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
                }

                k++;
            }

            j++;
        }
}

inline void moveToRDegeneracy( int vertex, 
                               int* vertexSets, int* vertexLookup, 
                               int** neighborsInP, int* numNeighbors,
                               //vector<vector<int>>& neighborsInP, int* numNeighbors,
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

inline void moveFromRToXDegeneracy( int vertex, 
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

#endif
