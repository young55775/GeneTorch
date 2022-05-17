import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from genetorch import reader


def search(df, col, kw):
    return df[col] == kw


def find(a):
    #  合并df
    reader.get_impact(a)
    a.candidate = []
    a.suppressor_group = []
    co_data = pd.DataFrame()
    for i in range(len(a.taglist)):
        co_data = pd.concat([co_data, a.taglist[i]])
        a.co_data = co_data
    a.result = get(co_data)
    print(a.result)
    return get(co_data)


# get为核心分析函数
def get(co_data):
    co_data.drop_duplicates()
    result = []
    index = list(set(list(co_data['gene'])))
    for i in range(len(index)):
        res = []
        ind = search(co_data, 'gene', index[i])
        co = pd.DataFrame(co_data.loc[ind, :])
        mut_ind = []
        var_ind = []
        var_fin = []
        var = []
        name = []
        for k, v in co.iterrows():
            name = v['gene']
            mut_ind.append(v['tag'])
            var_ind.append(v['protein'])
            var = list(set(var_ind))

        for q in range(len(var)):
            var_count_ind = search(co, 'protein', var[q])
            var_count = pd.DataFrame(co.loc[var_count_ind, :]).shape[0]
            var_fin.append([var[q], var_count])
        res.append([list(set(mut_ind)), len(set(mut_ind)), name, var_fin, len(var_fin)])
        result.extend(res)
    result = pd.DataFrame(result, columns=['sample', 'size', 'gene', 'variation', 'variation_number'])
    result = result[~result['variation'].isin([[]])]
    result = result.sort_values(by=['size'], ascending=False)
    result = result.reset_index(drop=True)
    return result


def filter(a, lengthlimit=0.6, rid=['ttn-1', 'cla-1']):
    reader.get_impact(a)
    a.candidate = []
    a.suppressor_group = []
    conc = pd.DataFrame()
    for i in range(len(a.taglist)):
        conc = pd.concat([conc, a.taglist[i]])
    filter = conc.groupby(['gene', 'protein', 'ID', 'base']).filter(lambda x: len(x) <= lengthlimit * len(a.taglist))
    a.co_data = filter
    result = get(filter)
    result = result[~result['gene'].isin(rid)]
    n = result.shape[0]
    result.index = list(range(n))
    a.result = result
    return result


def plot(result, size=800, c='Reds', title='variation_number', sub_c='RdPu', intersection_factor=0.03):
    list_var = []
    samplesize = 1.4 * max(result['variation_number'])
    m = result
    result2 = result[0:size]
    y = list(result2['size'])
    # 预设图像各种信息
    large = 22
    med = 16
    small = 12
    fig, ax = plt.subplots()
    xlist = list(result2['gene'])
    ylist = list(result2['variation_number'])
    slist = list(500 * result2['variation_number'] / max(result2['variation_number']))
    clist = list(0.7 * result2['variation_number'] / result2['size'])
    vlist = list(result2['variation'])
    samlist = list(result2['sample'])
    for i in range(len(vlist)):
        list_var = str(vlist[i]).split(',')
        if len(list_var) / 2 < 20:
            for j in range(0, len(list_var), 3):
                if j % 2 != 0:
                    list_var[j] = str(list_var[j]) + '\n'
            vlist[i] = ','.join(list_var)
        else:
            vlist[i] = 'variance number =' + str(len(list_var) / 2)

    for i in range(len(samlist)):
        list_sam = str(samlist[i]).split(',')
        if len(list_sam) < 20:
            for j in range(2, len(list_sam), 3):
                list_sam[j] = str(list_sam[j]) + '\n'
            samlist[i] = ','.join(list_sam)
        else:
            samlist[i] = 'sample size =' + str(len(list_sam))
    tlist2 = list(map(str, vlist))
    tlist1 = list(map(str, xlist))
    tlist3 = list(map(str, samlist))
    sc = plt.scatter(x=xlist, y=ylist,
                     s=slist,
                     alpha=0.7,
                     c=plt.get_cmap(c)(clist))
    print(sc.get_offsets()[1][0])
    annot = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="-",
                                        color="0.5",
                                        shrinkB=5,
                                        ),
                        fontsize=7)
    annot.set_visible(False)

    def update_annot(ind):

        pos = sc.get_offsets()[ind['ind'][0]]
        annot.xy = pos
        text = "{},\n{},\n{}".format("".join([tlist1[n] for n in ind["ind"]]),
                                     "".join([tlist2[n] for n in ind["ind"]]),
                                     "".join([tlist3[n] for n in ind["ind"]]))
        annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    # mouse event to get intersection
    def getintersection(event, color=sub_c, num=intersection_factor):
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                genename = str([tlist1[n] for n in ind["ind"]][0])
                print(ind['ind'])
                print([tlist1[n] for n in ind["ind"]])
                a = gsamplename(m, [tlist1[n] for n in ind["ind"]][0])
                print(a)
                b = intersect(m, a, i_num=num)
                print(b)
                plot(b, title=title + '---' + genename, c=color)

    fig.canvas.mpl_connect('button_press_event', getintersection)
    fig.canvas.mpl_connect("motion_notify_event", hover)
    # 装饰图像
    plt.gca().set(xlabel='mutation_frequency', ylabel='variation_number')
    plt.xticks([])
    plt.ylim(0, samplesize)
    plt.yticks(fontsize=12)
    plt.title(title, fontsize=12)
    plt.show()


def gsamplename(result, genename):
    if genename in list(result['gene']):
        ind = result[result['gene'].isin([genename])]
        return list(ind['sample'])


def intersect(result, samplename, i_num):
    result.reset_index()
    print(result)
    sn = list(samplename[0])
    sm = list(result['sample'])
    indlist = []
    for i in range(len(sm)):
        tmp = [a for a in sn if a in list(sm[i])]
        if len(tmp) > int(i_num * len(sn)):
            indlist.append(i)
    result = pd.DataFrame(result).drop(index=indlist, axis=0)
    n = result.shape[0]
    result.index = list(range(n))
    return result


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


def comb(m, n):
    if n != 0 and m - n != 0:
        b = re(m) / (re(n) * re(m - n))
    if n == 0 or m - n == 0:
        b = 1
    return b


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


def inters(list1, list2):
    tmp = [a for a in list1 if a in list2]
    return len(tmp)


def p_val(total, casea, caseb, n, snp):
    p = comb(casea, n) * comb(total - casea, caseb - n) / comb(total, caseb) * (snp - total) / (snp - casea)
    return p


def get_p_all(file, path):
    a = len(file.taglist)
    b = file.co_data.shape[0]
    m = []
    n = []
    sup_gene1 = list(file.result['sample'])
    sup_gene2 = list(file.result['sample'])
    inter = 0
    for i in range(len((sup_gene1))):
        for j in range(len(sup_gene2)):
            inter = inters(sup_gene1[i], sup_gene2[j])
            p = give_p(a, b, len(sup_gene1[i]), len(sup_gene2[j]), a, inter=inter)
            m.append(p)
        n.append(m)
    q = np.asarray(m).reshape(len(sup_gene1), len(sup_gene2))
    m = pd.DataFrame(q, index=list(file.result['gene']), columns=list(file.result['gene']))
    m.to_csv(path)


def get_p_sup(gene, file):
    a = len(file.taglist)
    b = file.co_data.shape[0]
    m = []
    n = []
    sup_gene = list(file.result['sample'])
    namelist = list(file.result['gene'])
    num = namelist.index(gene)
    inter = 0
    for i in range(len((sup_gene))):
        inter = inters(sup_gene[i], sup_gene[num])
        p = give_p(a, b, len(sup_gene[i]), len(sup_gene[num]), a, inter=inter)
        m.append(p)
    n.append(m)
    q = np.asarray(m).reshape(len(sup_gene), 1)
    m = pd.DataFrame(q, index=list(file.result['gene']), columns=['p_val'])
    return m


def merge(a, gene1, gene2):
    ga = a.result[a.result['gene'] == gene1]['sample'].tolist()[0]
    gb = a.result[a.result['gene'] == gene2]['sample'].tolist()[0]
    a.result = a.result[~a.result['gene'].isin([gene1, gene2])]
    lista = ga + gb
    lista = list(set(lista))
    size = len(lista)
    i = a.result.shape[0]
    a.result.loc[i + 2] = [lista, size, gene1 + '+' + gene2, '?', '?']
    print(a.result)
    return a.result


def split(a, candidate_list):
    candi = []
    no_candi = []
    for i in a.taglist:
        if list(i['tag'])[1] in candidate_list:
            candi.append(i)
        else:
            no_candi.append(i)
    return candi, no_candi


def candidate_list(a, gene_name):
    candidate = []
    for i in gene_name:
        candidate.extend(a.result[a.result['gene'] == i]['sample'].tolist()[0])
    candidate = list(set(candidate))
    return candidate


def important_gene(p_frame):
    return p_frame[p_frame['p_val'] < 0.05].index.tolist()


def classification(a, gene):
    p = get_p_sup(gene, a)
    lista = important_gene(p)
    candidate = candidate_list(a, lista)
    m = reader.readfile()
    n = reader.readfile()
    m.taglist, n.taglist = split(a, candidate)
    filter(m, 0.5)
    filter(n, 0.5)
    print(p)
    return m, n, p


def get_p_between(a, genea, geneb):
    if genea != geneb:
        sa = a.result[a.result['gene'] == genea]['sample'].tolist()[0]
        sb = a.result[a.result['gene'] == geneb]['sample'].tolist()[0]
        inter = inters(sa, sb)
        sup = len(a.taglist)
        snp = a.co_data.shape[0]
        p = give_p(sup, snp, len(sa), len(sb), sup, inter=inter)
    else:
        p = 1
    return p


def suppressor_group_p(a, genea):
    sa = a.result[a.result['gene'] == genea]['sample'].tolist()[0]
    sb = a.suppressor_group
    inter = inters(sa, sb)
    sup = len(a.taglist)
    snp = a.co_data.shape[0]
    p = p_val(len(a.taglist), len(sb), len(sa), inter, len(a.co_data))
    return p


def intersect_pro(a, genename):
    c = a
    c.result = c.result.reset_index(drop=True)
    lista = c.result['sample'].tolist()
    listb = c.result[c.result['gene'] == genename]['sample'].tolist()[0]
    c.candidate.append(genename)
    c.candidate = list(set(c.candidate))
    c.suppressor_group.extend(listb)
    c.suppressor_group = list(set(c.suppressor_group))
    print(c.candidate)
    print(len(c.suppressor_group))
    rid = []
    gene = []
    b = reader.readfile()
    for i in range(len(lista)):
        tmp = inters(lista[i], c.suppressor_group)
        if tmp != 0:
            p = get_p_between(c, genename, c.result['gene'].tolist()[i])
            if p > 0.05:
                rid.append(i)
                gene.append(c.result['gene'].tolist()[i])
    b.result = pd.DataFrame(c.result).drop(index=rid, axis=0)
    b.result = b.result[~b.result['gene'].isin(c.candidate)]
    b.taglist = c.taglist
    b.co_data = c.co_data[~c.co_data['gene'].isin(gene)]
    c = a
    return b, c


def plotpro(a, size=800, c='Reds', title='variation_number', sub_c='RdPu', origin=None):
    if origin is None:
        origin = a
    samplesize = 1.4 * max(a.result['variation_number'])
    result2 = a.result[0:size]
    fig, ax = plt.subplots()
    xlist = list(result2['gene'])
    ylist = list(result2['variation_number'])
    slist = list(500 * result2['variation_number'] / max(result2['variation_number']))
    clist = list(0.7 * result2['variation_number'] / result2['size'])
    vlist = list(result2['variation'])
    samlist = list(result2['sample'])
    for i in range(len(vlist)):
        list_var = str(vlist[i]).split(',')
        if len(list_var) / 2 < 20:
            for j in range(0, len(list_var), 3):
                if j % 2 != 0:
                    list_var[j] = str(list_var[j]) + '\n'
            vlist[i] = ','.join(list_var)
        else:
            vlist[i] = 'variation number =' + str(len(list_var) / 2)

    for i in range(len(samlist)):
        list_sam = str(samlist[i]).split(',')
        if len(list_sam) < 20:
            for j in range(2, len(list_sam), 3):
                list_sam[j] = str(list_sam[j]) + '\n'
            samlist[i] = ','.join(list_sam)
        else:
            samlist[i] = 'sample size =' + str(len(list_sam))
    tlist2 = list(map(str, vlist))
    tlist1 = list(map(str, xlist))
    tlist3 = list(map(str, samlist))
    sc = plt.scatter(x=xlist, y=ylist,
                     s=slist,
                     alpha=0.7,
                     c=plt.get_cmap(c)(clist))
    print(sc.get_offsets()[1][0])
    annot = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="-",
                                        color="0.5",
                                        shrinkB=5,
                                        ),
                        fontsize=7)
    annot.set_visible(False)

    def update_annot(ind):

        pos = sc.get_offsets()[ind['ind'][0]]
        annot.xy = pos
        text = "{},\n{},\n{}".format(",".join([tlist1[n] for n in ind["ind"]]),
                                     ",".join([tlist2[n] for n in ind["ind"]]),
                                     ",".join([tlist3[n] for n in ind["ind"]]))
        annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    # mouse event to get intersection
    def getintersection(event, color=sub_c):
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                genename = str([tlist1[n] for n in ind["ind"]][0])
                print(ind['ind'])
                b, c = intersect_pro(origin, genename)
                plotpro(b, title=title + '---' + genename, c=color, origin=c)

    fig.canvas.mpl_connect('button_press_event', getintersection)
    fig.canvas.mpl_connect("motion_notify_event", hover)
    plt.gca().set(xlabel='mutation_frequency', ylabel='variation_number')
    plt.xticks([])
    plt.ylim(0, samplesize)
    plt.yticks(fontsize=12)
    plt.title(title, fontsize=12)
    plt.show()
