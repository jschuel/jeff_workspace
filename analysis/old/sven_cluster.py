# file: cluster.py
# author: Sven Vahsen, Nov. 2018
#
# synopsis if main functions:
#
#   nclusters ()
#      clusters 3d hit distributions
#      each hit is specified by a list of three integers: [x,y,z]
#      returns a list of clusters
#              where each cluster is a list of 3d hits
#   n_gaps ()
#      finds gaps in x,y,z
#      return xlist, ylist, zlist
#             where each list contain one entry per gap, with the
#             value equal to the gap size
#
#   the other functions are helper functions.
#   find_gaps_1d might be useful if you only care about gaps / clusters
#   in one direction

from copy import deepcopy

def n_gaps (unclustered_hits):
    x = []
    y = []
    z = []
    for hit in unclustered_hits:
        x.append(hit[0])
        y.append(hit[1])
        z.append(hit[2])
    return find_gaps_1d (x), find_gaps_1d (y), find_gaps_1d (z)

def find_gaps_1d (positions):
    gaps = []
    gapsize = 0
    # move through all positions and search for gap
    for x in range(min(positions), max(positions)+1, 1):
        if not x in positions:
            gapsize +=1
        elif gapsize > 0:
            gaps += [gapsize]
            gapsize = 0
    return gaps
                 

# work properly if some pixels are present more than once?
# need to test that

def nclusters (unclustered_hits):
    return len (find_clusters(unclustered_hits))

def find_clusters (unclustered_hits, debug=False):
    clusters = []

    # make list of cluster numbers to make life easier
    
    # identify nearest neighbours of each hit
    # and create single-pixel seed clusters
    
    for i, hit in enumerate (unclustered_hits):
        clusters +=[[i]]
        neighbrs = find_neighbours (unclustered_hits)

    # combine clusters with common neighbours
    # until no more mergers happen
    
    

    while True:
        just_merged_so_keep_going = False    
        if debug:
            print ("clusters=", clusters)
            print ("neighbrs=", neighbrs)
        for k, (neig, clus) in enumerate (zip(neighbrs,clusters)):
            for j, (neig2,clus2) in enumerate (zip(neighbrs,clusters)):
                if j != k:
                    # is current cluster, clus2, is a neighbour to cluster clus?
                    if k in neig2:
                        # yes --> remove neighbour entries from clus2
                        # and then add merge remaining cluster and neighbour entries
                        neig2.remove(k)
                        neig += deepcopy(neig2)
                        clus += deepcopy(clus2)
                        #if debug:
                        #    print ("clusters=", clusters)
                        #    print ("neighbrs=", neighbrs)
                        neig2.clear()
                        clus2.clear()
                        just_merged_so_keep_going = True
        if just_merged_so_keep_going == False: 
            break

    # now we have the final clusters

    clusters = remove_empty_clusters(clusters)
    return clusters

    # unpack them from pixel number to coordinates
    # and return the results

# got this from stackoverflow.com
def remove_empty_clusters (list1):
    list2 = [x for x in list1 if x != []]
    return list2

# for each pixel, create list of direct neighbours
def find_neighbours (unclustered_hits):
    neighbours = []
 
    for hit in unclustered_hits:
        current_neighbours = []
        for j, hit2 in enumerate (unclustered_hits):
            if hit2 != hit:
                if are_neighbours (hit,hit2):
                    current_neighbours += [j]
        neighbours += [current_neighbours]
    return neighbours

# this function defines being a nearest neigbour
# note: currently allowing diagonal neigbours in 2D or 2D

def are_neighbours (hit,hit2):
    if abs (hit[0]-hit2[0])<=1 and abs (hit[1]-hit2[1]) <=1 and abs (hit[2]-hit2[2]) <= 1:
        return True
    else:
        return False

