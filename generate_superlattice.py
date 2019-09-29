import numpy as np
import generate_HNF
import spglib
from copy import deepcopy

rutile = ([(4, 0, 0),

           (0, 4, 0),

           (0, 0, 3)],

          [(0, 0, 0),

           (0.5, 0.5, 0.5),

           (0.3, 0.3, 0.0),

           (0.7, 0.7, 0.0),

           (0.2, 0.8, 0.5),

           (0.8, 0.2, 0.5)],

          [14, 14, 8, 8, 8, 8])
o_sublattice = ([(4, 0, 0),

           (0, 4, 0),

           (0, 0, 3)],

          [(0.3, 0.3, 0.0),

           (0.7, 0.7, 0.0),

           (0.2, 0.8, 0.5),

           (0.8, 0.2, 0.5)],

          [8, 8, 8, 8])
parent_sym = spglib.get_symmetry(rutile)
HNF_list = generate_HNF.generate_all_superlattices(2)
reduced_HNF = generate_HNF.reduce_HNF_list_by_parent_lattice_symmetry(HNF_list, parent_sym["rotations"])

def get_superlattice_vectors(o_sublattice, HNF):
    lattice = np.asarray(o_sublattice[0])
    lattice = HNF.dot(lattice)
    lattice = lattice.tolist()
    return lattice 

def zurasu(vector, motono):
    t = deepcopy(motono)
    for i in range(len(motono)):
        t[i] += vector
    d = np.concatenate([motono,t])
    return d

def tukuru(o_sublattice, HNF):  #need repair#########################################################################
    a = np.sum(HNF, axis=0)
    def zurasu(vector, motono):
        t = deepcopy(motono)
        for i in range(len(motono)):
            t[i] += vector
        d = np.concatenate([motono,t])
        return d
    new_zahyou = deepcopy(o_sublattice[1])
    for i in range(3):
        for j in range(1, a[i]):# a[i]が整数行列でないとエラーを起こす。　dtype=np.int で解決
            vec = np.zeros(3)
            vec[i] = j
            new_zahyou = zurasu(vec, new_zahyou)
    
    g = deepcopy(new_zahyou)
    renew = np.linalg.inv(HNF) #逆行列 新しいsuperlatticevector で分立座標にあらわせるように　
    for i in range(len(new_zahyou)):
        g[i] = g[i].dot(renew)
        
    #消すべき行を記録
    h = []
    for kesu in range(len(g)):
        if  0 <= g[kesu][0] < 1 and 0 <= g[kesu][1] < 1 and 0 <= g[kesu][2] < 1 :
            continue
        else:
            h.append(kesu)
        
    new_zahyou = np.delete(g, h, 0)
    return new_zahyou


#superlattice wo tukuru
def get_superlattice(parent_lattice, HNF, index):
    def get_parent_gensi(index): #要修正　　
        empty = []
        for i in range(index*2):
            empty.append(14)
        for j in range(index*4):
            empty.append(8)
        return empty
    a = [[], [], []] # a は　新しいsuperlattice
    a[0] = get_superlattice_vectors(parent_lattice, HNF)
    a[1] = tukuru(parent_lattice, HNF)
    a[1] = np.round(a[1], 4)
    a[1] = a[1].tolist()
    a[2] = get_parent_gensi(index)
    a = tuple(a)
    return a
#superlattice no o_sublattice wo tukuru
def get_o_subsuperlattice(o_sublattice, HNF, index):
    a = [[], [], []]
    a[0] = get_superlattice_vectors(o_sublattice, HNF)
    a[1] = tukuru(o_sublattice, HNF)
    a[1] = np.round(a[1], 4)
    a[1] = a[1].tolist()
    def get_o_subsuper_gensi(index):
        empty = []
        for j in range(index*4):
            empty.append(8)
        return empty
    a[2] = get_o_subsuper_gensi(index)
    a = tuple(a)
    return a