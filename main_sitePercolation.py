""" main_sitePercolation.py

Contains function defintions for union find based percolation algorithm. 
Albeit the implementation is slightly different, it is based on the idea
of the fast percolation algorithm detailed in

        "Fast Monte Carlo algorithm for site and bond percolation"
        Newman, M.E.J. and Ziff, R.M.
        PRE 64 (2001) 016706

Author : Oliver Melchert
Date   : 07/28/2015 
"""
import sys
import os
import gzip
from random import random, seed, shuffle


# MISC 1{{{

class Graph(object):
        pass


def shuffleSites(myList):
        """iterator over shuffled list elements"""
        shuffle(myList)
        ctr = 0
        for x in myList:
            ctr += 1
            yield ctr, x

# 1}}}

# READ CFG 1{{{

def fetchGraph_Zheng(fName):
        """construct adjacengy list graph from edge list in DIMACS formatted file"""
        G     = Graph()
        adjList = dict()
        eCtr = 0
        with gzip.open(fName,'r') as f:
            for line in f:
                c = line.split()
                if len(c)==1:
                        G.v = int(c[0])
                        adjList = {i:[] for i in range(G.v)}
                elif len(c[0])==3:
                        vi = int(c[0])
                        vj = int(c[1])
                        eCtr += 1

                        #if vi not in adjList: adjList[vi]=[]
                        adjList[vi].append(vj)

                        #if vj not in adjList: adjList[vj]=[]
                        adjList[vj].append(vi)

        G.e = eCtr
        G.adjList = adjList
        return G

def fetchGraph(fName):
        """construct adjacengy list graph from edge list in DIMACS formatted file"""
        G     = Graph()
        adjList = dict()
        with gzip.open(fName,'r') as f:
            for line in f:
                c = line.split()
                if c[0]=='p':
                        G.v   = int(c[2])
                        G.e   = int(c[3])
                        eType = int(c[4])

                elif c[0]=='e':
                        vi = int(c[1])
                        vj = int(c[2])

                        if vi not in adjList: adjList[vi]=[]
                        adjList[vi].append(vj)

                        if eType==0:
                            if vj not in adjList: adjList[vj]=[]
                            adjList[vj].append(vi)

        G.adjList = adjList
        return G
# 1}}}

# UNION-FIND 1{{{

class UnionFind_byRank_pathCompression(object):
        """ union find data structure implementing union-by-rank 
        and path compression

        NOTE: 
        -# keeps track of largest cluster size
        -# keeps track of the number of components
        """
        def __init__(self):
                self.par   = dict()
                self.size  = dict()
                self.rank  = dict() 

                self.nComp = 0 
                self.cMax  = 0


        def newElement(self,i):
                """initialize new element"""
                self.par[i]  = i
                self.size[i] = 1
                self.rank[i] = 1
                self.nComp  += 1


        def find(self,i):
                """find root node of element and perform compression"""
                if self.par[i]!=self.par[self.par[i]]: 
                        self.par[i] = self.find(self.par[i])
                return self.par[i]

        
        def union(self,i,j):
                """union two elements by rank"""
                ii, jj = self.find(i), self.find(j)
                if ii != jj:

                        if self.rank[ii] > self.rank[jj]: 
                                ii,jj = jj,ii

                        if self.rank[ii] == self.rank[jj]: 
                                self.rank[jj]+=1

                        self.size[jj] += self.size[ii]
                        self.par[ii]   = jj

                        self.cMax   = max(self.size[jj],self.cMax)
                        self.nComp -=1

                        del self.size[ii]
# 1}}}


def main_sitePercolation_unionFind():

        bondFile = sys.argv[1]
        mySeed   = int(sys.argv[2])

        seed(mySeed)
        print "# bondFile = %s"%(bondFile)
        print "# seed = %d"%(mySeed)
        print "# (nSites) (G.v) (nComp) (cMax)"
        G   = fetchGraph(bondFile)
        uf  = UnionFind_byRank_pathCompression()
        occ = [0 for i in range(G.v)]
        for nSites,i in shuffleSites(range(G.v)):
            uf.newElement(i); occ[i] = 1
            for j in G.adjList[i]:
                if occ[j]: uf.union(i,j)
            print "%4d %4d %4d %4d"%(nSites, G.v, uf.nComp, uf.cMax)


main_sitePercolation_unionFind()
# EOF: main_sitePercolation.py
