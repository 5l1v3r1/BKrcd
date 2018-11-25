
#ifndef _DJS_NEWHYBRID_ALGORITHM_H_
#define _DJS_NEWHYBRID_ALGORITHM_H_

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

using namespace std;

class newHybrid : public Algorithm
{
public:
    newHybrid(std::vector<std::list<int>> const &adjacencyList, double _density_threshold);
    virtual ~newHybrid();

    virtual long Run(std::list<std::list<int>> &cliques);

    newHybrid           (newHybrid const &) = delete;
    newHybrid& operator=(newHybrid const &) = delete;

    void BKReverseRecursive(long long* cliqueCount,
                            std::list<int> &partialClique, 
                            int* vertexSets, int* vertexLookup,
                            vector<vector<int>>& neighborsInP, int* numNeighbors,
                            int beginX, int beginP, int beginR);

    void listAllMaximalCliquesDgcyRecursive(long long* cliqueCount,
                                               list<int> &partialClique, 
                                               int* vertexSets, int* vertexLookup,
                                               vector<vector<int>>& neighborsInP, int* numNeighbors,
                                               int beginX, int beginP, int beginR);


    long long newHybridMain(std::vector<std::list<int>> const &adjList, int size);

private:
    std::vector<std::list<int>> const &m_AdjacencyList;
    double density_threshold;
};

#endif

