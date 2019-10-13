#generate permutation
import spglib
import numpy as np
import pymatgen
from copy import deepcopy

def gene_per(rotation, translation, dic_zahyou):
    """
    対象操作について置換の移り先を示す辞書を作成

    parameters

    retruns

    """
    dic2 = deepcopy(dic_zahyou)
    for l in range(len(dic_zahyou)):
        dic2[l] = (rotation.dot(dic_zahyou[l]) + translation) % 1
    tin = dict()
    for l in range(len(dic_zahyou)):
        for j in range(len(dic_zahyou)):
            if all(np.round(dic_zahyou[j], 4) == np.round(dic2[l], 4)):
                tin[l] = j
    return tin

def generate_per(superlattice, o_sublattice):
    parent_sym = spglib.get_symmetry(superlattice)
    tin = dict()
    arr = np.asarray(o_sublattice[1])   #lattcie no 行列表示
    dic = dict()                        #分率座標の行列表示
    for i in range(len(o_sublattice[1])):
        dic[i] = np.asarray(o_sublattice[1][i])

    jun = dict()
    for i in range(len(parent_sym["translations"])):
        jun[i] = gene_per(parent_sym["rotations"][i], parent_sym["translations"][i], dic)
    return jun

# 置換の積の式を作る [{0, 1, 2, 3}, {4, 5, 6, 7}]
def junkai(tin):
    c = [0, 1, 2, 3, 4, 5, 6, 7]
    d = []
    f = tin[0]
    k = -1
    for i in range(8):
        if i in c:
            d.append([])
            k += 1
            f = tin[i]
            d[k].append(i)
            while f in c:
                d[k].append(f)
                c.remove(f)
                f = tin[f]
    for i in range(len(d)):
        d[i] = set(d[i])
    return d
#置換の積の式の辞書を作成
def tikan():
    jun = dict()
    for i in range(len(o_subsim["translations"])):
        tin = gene_per(i)
        jun[i] = junkai(tin)
    return jun

#置換の積の式からpolyaの数え上げの定理を使って配色のパターン数を出す
def polya(jun):
    sum = 0
    for i in range(len(jun)):
        sum += (number_of_colors)**(len(jun[i]))
    polya = sum / len(jun)
    return polya

#unique
def unique(o_sublattice, parent_sym_jun, omomi4):
    omomi4_neo = deepcopy(omomi4)#copy wo sakusei
    lis = set()

    for i in range(len(omomi4)):#through all candidate
        if omomi4[i] in omomi4_neo:#kouho ga mada 生き残ってるかチェック
            for j in range(1, len(parent_sym_jun)):#置換操作について回す　
                d = dict()
                sta = ""
                for k in range(len(o_sublattice[2])):
                    d[parent_sym_jun[j][k]] = omomi4[i][k]#str辞書の作成
                for r  in range(len(o_sublattice[2])):#staの作成 sta is made from tikan[j]
                    sta += d[r]
                if int(sta) is not int(omomi4[i]):
                    if sta in omomi4_neo:
                        omomi4_neo.remove(sta)
                        lis.add(i)
    manko = []
    can = omomi4
    for i in range(len(lis)):
        lis = list(lis)
        a = lis[i]
        if can[a][0] == can[a][4] and can[a][1] == can[a][5] and can[a][2] == can[a][6] and can[a][3] == can[a][7]:
            continue
        else:
            manko.append(can[a])
    return manko

def make_candidate(o_sublattice, l):
    omomi4 = []
    for i in range(2**len(o_sublattice[2])):
        k = format(i,'b').zfill(len(o_sublattice[2]))
        sum = 0
        for m in range(len(k)):
            sum += int(k[m])
        if sum ==l :
            omomi4.append(k)
    return omomi4






