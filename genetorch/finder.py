import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def search(df, col, kw):
    return df[col] == kw


def find(taglist):
    #  合并df
    co_data = pd.DataFrame()
    for i in range(len(taglist)):
        co_data = pd.concat([co_data, taglist[i]])
    print(co_data)
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
        res.append([set(mut_ind), len(set(mut_ind)), name, var_fin, len(var_fin)])
        result.extend(res)
    result = pd.DataFrame(result, columns=['sample', 'size', 'gene', 'variation', 'variation_number'])
    result = result[~result['variation'].isin([[]])]
    result = result.sort_values(by=['size'], ascending=False)
    result = result.reset_index(drop=True)
    return result


def filter(taglist, lengthlimit=0.6, threshold=1.05, rid=['ttn-1', 'cla-1']):
    conc = pd.DataFrame()
    for i in range(len(taglist)):
        conc = pd.concat([conc, taglist[i]])
    filter = conc.groupby(['gene', 'protein', 'ID', 'base']).filter(lambda x: len(x) <= lengthlimit * len(taglist))
    n = threshold * len(taglist)
    result = get(filter)
    if int(n) < max(list(result['variation_number'])) + 1:
        for i in range(int(n), max(list(result['variation_number'])) + 1):
            plus = pd.DataFrame(result.loc[result['variation_number'] == i, :])
            result = pd.concat([result, plus])
        result.drop_duplicates(subset=['gene'], keep=False, inplace=True)
    result = result[~result['gene'].isin(rid)]
    n = result.shape[0]
    result.index = list(range(n))
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
    plt.gca().set(xlabel='gene', ylabel='variation_number')
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
        if len(tmp) > int(i_num*len(sn)):
            indlist.append(i)
    result = pd.DataFrame(result).drop(index=indlist, axis=0)
    n = result.shape[0]
    result.index = list(range(n))
    return result
