import sys
import itertools

def KnnBasedChimeraGraph_eList(n,Lx,Ly):
        E = []
        for x,y in itertools.product(range(Lx),range(Ly)):
             i0 = 2*n*(x+y*Lx)
             # inner edges of current K_nn
             E += list(itertools.product(range(i0,i0+n),range(i0+n,i0+2*n)))
             # interconnecting K_nn (K prime)_nn edges
             for ni in range(n):
                if x < Lx-1:
                        E.append((i0+n+ni,i0+2*n+n+ni))
                if y < Ly-1:
                        E.append((i0+ni,i0+2*n*Lx+ni))
        # full edge set of K_nn based lattice graph
        return E 


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
