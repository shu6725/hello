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

def gene_trans(HNF):
    """
    HNFから並進操作の一覧を作成

    parameters:HNF

    retruns

    """
    int_HNF = HNF.astype('int64')
    vectors = np.zeros(3)
    hantei = True
    for dimention in range(3):
        cp_vecs = deepcopy(vectors)
        if HNF[dimention][dimention] > 1:
            for king in range(1, HNF[dimention][dimention]):
                cp_vecs2 = deepcopy(cp_vecs)
                if hantei:
                    cp_vecs2[dimention] += king/HNF[dimention][dimention]
                else:
                    for num_vecs in range(cp_vecs2.shape[0]):
                        cp_vecs2[num_vecs][dimention] += king/HNF[dimention][dimention]
                vectors = np.vstack((vectors, cp_vecs2))
            hantei = False
    return vectors

def get_trans_permuation(dic_zahyou, translation):
    """
    並進操作について置換の移り先を示す辞書を作成

    parameters

    retruns

    """
    tin = dict()
    dic2 = deepcopy(dic_zahyou)
    for l in range(len(dic_zahyou)):
        dic2[l] = (dic_zahyou[l] + translation)
        dic2[l] = np.round(dic2[l], 3)%1
    for l in range(len(dic_zahyou)):
        for j in range(len(dic_zahyou)):
            if all(np.round(dic_zahyou[j], 3) == np.round(dic2[l], 3)):
                tin[l] = j
    return tin

def get_trans_perms(dic_zahyou, translations):
    """
    並進操作について置換の移り先を示す辞書の辞書を作成

    parameters

    retruns

    """
    trans_perms = dict()
    for i in range(len(translations)):
        trans_perms[i] = get_trans_permuation(dic_zahyou, translations[i])
    return trans_perms

def kumiawase(rot_tikan, trans_tikan):
    """
    回転操作と並進操作の置換を組合す

    parameters

    retruns

    """
    hosii = dict()
    for i in range(len(rot_tikan)):
        hosii[i] = trans_tikan[rot_tikan[i]]
    return hosii

def generate_abs_permuatation(parent_lattice_jun, trans_perms):
    """
    回転操作と並進操作の置換を組合せた究極の辞書を作成

    parameters

    retruns

    """
    tikan_list = list()
    for i in parent_lattice_jun:
        for j in trans_perms:
            if i == 1:
                print(j, trans_perms[j])
            tikan = kumiawase(parent_lattice_jun[i], trans_perms[j])
            if i == 1:
                print(tikan)
            tikan_list.append(tikan)
    return tikan_list

def generate_per(superlattice, o_sublattice):##各座標をnumpyに変換して、gene_perを使う
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

def unique(parent_sym_jun, omomi4):
    """
    得られた究極の対称操作の辞書に基づいて、配列をユニークしていく

    parameters

    parent_sym_jun

    omomi4 は各空孔数ごとのcandidate集合

    retruns

    """
    omomi4_neo = deepcopy(omomi4)#copy wo sakusei
    lis = set()

    for i in range(len(omomi4)):#through all candidate
        if omomi4[i] in omomi4_neo:#kouho ga mada 生き残ってるかチェック
            for j in range(1, len(parent_sym_jun[0])):#置換操作について回す　
                d = dict()
                sta = ""
                for k in range(len(parent_sym_jun[0])):
                    d[parent_sym_jun[j][k]] = omomi4[i][k]#str辞書の作成
                for r  in range(len(parent_sym_jun[0])):#staの作成 sta is made from tikan[j]
                    sta += d[r]
                if int(sta) is not int(omomi4[i]):
                    if sta in omomi4_neo:
                        omomi4_neo.remove(sta)
                        lis.add(omomi4[i])

    return lis









