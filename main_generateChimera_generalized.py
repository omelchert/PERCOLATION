import sys
import itertools


def KnnBasedChimeraGraph_eList(n,Lx,Ly):

        def V_K_nn(ctr,n):
                """ordered pairs of K_nn vertices
                
                ctr -- id of the first vertex in current K_nn
                n   -- size of K_nn
                """
                return zip(range(ctr,ctr+n),
                        range(ctr+n,ctr+2*n))

        def E_K_nn(ctr,n):
                """edge set of K_nn
                
                ctr -- id of the first vertex in current K_nn
                n   -- size of K_nn
                """
                return list(itertools.product(range(ctr,ctr+n),
                        range(ctr+n,ctr+2*n)))

        # arrange ordered pairs of K_nn subgraphs in Lx x Ly grid
        V = { (x,y): V_K_nn(2*n*(x+y*Ly),n) 
               for x,y in itertools.product(range(Lx),range(Ly)) }

        # accumulate edges
        innerEdges = []
        interEdges = []
        for x,y in itertools.product(range(Lx),range(Ly)):
             # inner edges of current K_nn
             innerEdges += E_K_nn(2*n*(x+y*Ly),n)
             # interconnecting K_nn (K prime)_nn edges
             for ni in range(n):
                if x < Lx-1:
                        interEdges.append((V[(x,y)][ni][1],V[(x+1,y)][ni][1]))
                if y < Ly-1:
                        interEdges.append((V[(x,y)][ni][0],V[(x,y+1)][ni][0]))

        # full edge set of K_nn based lattice graph
        return innerEdges + interEdges


def main():

        Lx = int(sys.argv[1])
        Ly = int(sys.argv[2])
        n  = int(sys.argv[3])


        E = KnnBasedChimeraGraph_eList(n,Lx,Ly) 

        print "c K_{%d,%d} based (Chimera-like) lattice graphs"%(n,n)
        print "c grid layout with side lengths Lx = %d, Ly = %d"%(Lx,Ly)
        print "p edge %d %d 0"%(2*n*Lx*Ly,len(E))
        for e in sorted(E):
                print "e %d %d"%(e[0],e[1])
        

main()
