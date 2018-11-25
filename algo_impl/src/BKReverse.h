
#ifndef _DJS_BKREVERSE_ALGORITHM_H_
#define _DJS_BKREVERSE_ALGORITHM_H_

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
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <unordered_map>
#include <unordered_set>

using namespace std;

class BKReverse : public Algorithm
{
public:
    BKReverse(std::vector<std::list<int>> const &adjacencyList);
    virtual ~BKReverse();

    virtual long Run(std::list<std::list<int>> &cliques);

    BKReverse           (BKReverse const &) = delete;
    BKReverse& operator=(BKReverse const &) = delete;

    void BKReverseRecursive(long long* cliqueCount,
                            std::list<int> &partialClique, 
                            int* vertexSets, int* vertexLookup,
                            //int** neighborsInP, int* numNeighbors,
                            vector<vector<int>>& neighborsInP, int* numNeighbors,
                            int beginX, int beginP, int beginR);

    long long BKReverseMain(std::vector<std::list<int>> const &adjList, int size);

private:
    std::vector<std::list<int>> const &m_AdjacencyList;
};

inline bool existCommonNeighborOfPinX( int* vertexSets, int *vertexLookup, 
                                       //int** neighborsInP, int* numNeighbors,
                                       vector<vector<int>>& neighborsInP, int *numNeighbors, 
                                       int beginX, int beginP, int beginR)
{
    int j = beginX, sizeofP = beginR-beginP;
    while (j < beginP) {
        int numNeighborsInP = 0;
        int vertex = vertexSets[j];
        //int numPotentialNeighbors = min(sizeofP, numNeighbors[vertex]);
        int numPotentialNeighbors = numNeighbors[vertex];

        int k = 0;
        while(k < numPotentialNeighbors)
        {
            int neighbor = neighborsInP[vertex][k];
            int neighborLocation = vertexLookup[neighbor];
            if(neighborLocation >= beginP && neighborLocation < beginR)
                numNeighborsInP++;
            /*
            else
                break;
                */
            k++;
        }
        if (numNeighborsInP == sizeofP) 
            return true;
        ++j;
    }
    return false;
}

/*
 * return true if the remaining vertices in P is a clique
 */

inline bool findMostNonNeighbors ( int* vertexSets, int *vertexLookup, 
                                   vector<vector<int>>& neighborsInP, int *numNeighbors, 
                                   //int **neighborsInP, int *numNeighbors,
                                   int beginX, int beginP, int beginR,
                                   int *pivot)
{
    int j       = beginP;
    int degree  = beginR-beginP;
    int sizeofP = beginR-beginP;

    while (j < beginR) {
        int numNeighborsInP = 0;
        int vertex = vertexSets[j];
        //int numPotentialNeighbors = min(sizeofP, numNeighbors[vertex]);
        int numPotentialNeighbors = numNeighbors[vertex];

        int k = 0;
        while(k<numPotentialNeighbors)
        {
            int neighbor = neighborsInP[vertex][k];
            int neighborLocation = vertexLookup[neighbor];

            if(neighborLocation >= beginP && neighborLocation < beginR)
                numNeighborsInP++;
            k++;
        }

        if (numNeighborsInP < degree) {
            *pivot = vertexSets[j];
            degree = numNeighborsInP;
        } 

        ++j;
    }

    if (degree == beginR-beginP-1)
        return true;

    return false;
}

inline bool newfindMostNonNeighbors ( int* vertexSets, int *vertexLookup, 
                                   //int **neighborsInP, int *numNeighbors,
                                   vector<vector<int>>& neighborsInP, int *numNeighbors, 
                                   int beginX, int beginP, int beginR,
                                   int *pivot, map<int, unordered_set<int>>& pivots, unordered_map<int,int>& mapping)
{
    int j       = beginP;
    int degree  = beginR-beginP;
    int sizeofP = beginR-beginP;

    for (auto it = pivots.begin(); it != pivots.end(); ++it) {
        if (!it->second.empty()) {
            *pivot = *(it->second.begin());
            it->second.erase(it->second.begin());
            degree = 0;
            for (int j = 0; j < numNeighbors[*pivot]; ++j) {
                int neighbor = neighborsInP[*pivot][j];
                int neighborLocation = vertexLookup[neighbor];
                if (neighborLocation >= beginP && neighborLocation < beginR) {
                    degree++;
                    const int& neighbordegree = mapping[neighbor];
                    pivots[neighbordegree].erase(neighbor);
                    pivots[neighbordegree-1].insert(neighbor);
                    --mapping[neighbor];
                }
            }
            break;
        }
    }

    if (degree == sizeofP-1)
        return true;

    return false;
}

inline void moveToR( int vertex, 
                     int* vertexSets, int* vertexLookup, 
                     //int **neighborsInP, int *numNeighbors,
                     vector<vector<int>>& neighborsInP, int* numNeighbors,
                     int* pBeginX, int *pBeginP, int *pBeginR, 
                     int* pNewBeginX, int* pNewBeginP, int *pNewBeginR)
{
        int vertexLocation = vertexLookup[vertex];

        /* move vertex into R */
        (*pBeginR)--;
        swap(vertexLookup[vertexSets[*pBeginR]], vertexLookup[vertex]);
        swap(vertexSets[vertexLocation], vertexSets[*pBeginR]);

        // this is not a typo, initially newX is empty
        *pNewBeginX = *pBeginP;
        *pNewBeginP = *pBeginP;
        *pNewBeginR = *pBeginP;

        int sizeOfP = *pBeginR - *pBeginP;
        int j = *pBeginX;

        j = (*pBeginX);
        while(j < *pNewBeginX)
        {
            int neighbor = vertexSets[j];
            int neighborLocation = j;

            int incrementJ = 1;
            //int numPotentialNeighbors = min(sizeOfP+1, numNeighbors[neighbor]);
            int numPotentialNeighbors = numNeighbors[neighbor];

            int k = 0;
            while(k<numPotentialNeighbors)
            {
                if(neighborsInP[neighbor][k] == vertex)
                {
                    (*pNewBeginX)--;
                    swap(vertexLookup[neighbor], vertexLookup[vertexSets[*pNewBeginX]]);
                    swap(vertexSets[neighborLocation], vertexSets[*pNewBeginX]);
                    incrementJ=0;
                    break;
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

            //int numPotentialNeighbors = min(sizeOfP+1, numNeighbors[neighbor]);
            int numPotentialNeighbors = numNeighbors[neighbor];

            int k = 0;
            while(k<numPotentialNeighbors)
            {
                if(neighborsInP[neighbor][k] == vertex)
                {
                    swap(vertexLookup[neighbor], vertexLookup[vertexSets[*pNewBeginR]]);
                    swap(vertexSets[neighborLocation], vertexSets[*pNewBeginR]);
                    (*pNewBeginR)++;
                    break;
                }
                k++;
            }
            j++;
        }
}

inline void moveFromRToX( int vertex, 
                          int* vertexSets, int* vertexLookup,               
                          //int**  neighborsInP, int* numNeighbors,
                          vector<vector<int>>& neighborsInP, int* numNeighbors,
                          int* pBeginX, int* pBeginP, int* pBeginR )
{
    int vertexLocation = vertexLookup[vertex];

    //swap vertex into X and increment beginP and beginR
    swap(vertexLookup[vertexSets[*pBeginP]], vertexLookup[vertex]);
    swap(vertexSets[vertexLocation], vertexSets[*pBeginP]);

    ++*pBeginP;
    ++*pBeginR;
}

#endif

