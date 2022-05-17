import pandas as pd
import numpy as np
import random
from genetorch import finder


def creat_list(h, n):
    return random.sample(range(1, h + 1), n)


def fake_data(a, gene_name):
    dist = list(a.result[a.result['gene'] != gene_name]['size'])
    h = len(a.taglist)
    faker = []
    for i in dist:
        faker.append(list(creat_list(h, i)))
    return faker


def inters(list1, list2):
    tmp = [a for a in list1 if a in list2]
    return len(tmp)


def give_p(snp_sup, snp_all, num_gene1, num_gene2, num_geneall, inter):
    if num_gene1 >= num_gene2:
        p1 = num_gene1 / num_geneall
        p2 = num_gene2 / num_geneall
    else:
        p1 = num_gene2 / num_geneall
        p2 = num_gene1 / num_geneall
    c_num = comb(num_geneall, inter)
    pdouble = p1 * p2
    p = c_num * (pdouble ** inter) * ((1 - pdouble) ** (num_geneall - inter)) * (
            (snp_all - snp_sup) / (snp_all - num_gene1))
    return p


def comb(m, n):
    if n != 0 and m - n != 0:
        b = re(m) / (re(n) * re(m - n))
    if n == 0 or m - n == 0:
        b = 1
    return b


def re(n):  # 阶乘
    x = 1
    if n == 0:
        return 0
    elif n == 1:
        return 1
    else:
        for i in range(1, n + 1):
            x *= i
        return x


def test(snp, sup, dist, probe, threshold):
    p2 = (snp - sup) / (snp - len(probe))
    sup_gene = dist
    m = []
    for i in range(len((sup_gene))):
        inter = inters(sup_gene[i], probe)
        p = finder.p_val(sup, len(sup_gene[i]), len(probe), inter, snp) * p2
        if p <= threshold:
            m.append(p)
    return len(m)


def false_positive(a, gene_name, threshold=None):
    if threshold is None:
        threshold = [0.05, 0.04, 0.03, 0.02, 0.01]
    snp = len(a.co_data)
    sup = len(a.taglist)
    probelen = int(a.result[a.result['gene'] == gene_name]['size'])

    m = []
    for i in range(100):
        probe = list(random.sample(range(1, sup), probelen))
        faked = fake_data(a, gene_name)
        line = []
        for j in threshold:
            p = test(snp, sup, faked, probe, j)
            line.append(p)
            m.append(line)
    print(m)
    mat = np.asarray(m)
    mat = mat.transpose()
    pre = pd.DataFrame(index=threshold, data=mat)
    pre['average'] = pre.mean(axis=1)
    print(pre['average'])
    return pre
