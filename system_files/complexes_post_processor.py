import copy
from system_files.extract_ssdna_from_data import extract_ssdna


# post-process of the given file
# source_path, file_name, des(destination path) are necessary arguments parsed by browse_and_parse_v03.py
def complexes_post_process(source_path, file_name, des, adv):

    ssdna_and_n = extract_ssdna(source_path, 'run_vis', '')

    data_set = list(ssdna_and_n[0])
    n_complexes = list(ssdna_and_n[1])

    # add new ID for each complex (as appeared on a single line on the species file)
    def id_data():
        id_set = []
        for i, e in zip(range(1, len(data_set) + 1), data_set):
            e.insert(0, [str(i)])
            id_set.append(e)
        return id_set

    id_dat = list(id_data())

    # create a dictionary of data by adding IDs to single list in a sub-list
    def dic_id_data():
        sd = []
        for i in id_dat:
            a = dict(enumerate(i, start=1))
            sd.append(a)
        return sd

    # lexicographical sorting by replacing 'x' to complimented agents' labels
    # return a list of their IDs
    def lex_sorting():
        p = []
        for vvv in dic_id_data():
            s = {}
            p.append(s)
            for k, v in vvv.items():
                    n_list = [x[0] + 'x' if x[0].isalpha()
                              and x[1:].isdigit() else x for x in v]
                    n_dic = {k: n_list}
                    s.update(n_dic)

        lexi_sorte = []
        for i in p:
            l_s = dict(sorted(i.items(), key=lambda x: x[1]))
            lexi_sorte.append(list(l_s.keys()))
        return lexi_sorte

    so_li_keys = list(lex_sorting())

    # lex_sorting() method returns the IDs of the lists, after lexicographically sorting
    # arrange the lists according to sorted IDs and create a new list
    def sorted_id():
        lex_f = []
        for c, cc in zip(so_li_keys, dic_id_data()):
            n_li = [cc[x] for x in c]
            lex_f.append(n_li)
        return lex_f

    lexi_sorted = list(sorted_id())

    # reset the IDs of the sub-lists ascending manner
    def re_id_data():
        id_set_comp_2, id_set_raw_2,  = [], []

        for c in lexi_sorted:
            id_set_comp_1, id_set_raw = [], []
            id_set_comp_1.append(c[0])
            id_set_raw.append(c[0])
            id_set_comp_2.append(id_set_comp_1)
            id_set_raw_2.append(id_set_raw)
            count = []
            count.append(str(0))
            for e in c[1:]:
                raw = [x[0] for x in e]
                n = int(count[-1]) + 1
                new_li_raw = [[str(c[0][0]) + '.' + str(n)], raw]
                new_li_comp_1 = [[str(c[0][0]) + '.' + str(n)], e]
                count.append(n)
                id_set_raw.append(new_li_raw)
                id_set_comp_1.append(new_li_comp_1)

        return id_set_comp_2

    re_ids = list(re_id_data())

    # recognise the compliment agent and complimented
    # create a new list of that
    def bind():
        comp_dic = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        orie_dic = {'5': '3', '3': '5'}
        ll, q, temp = [], [], []

        def look(li, co):
            for k in li:
                for i in k[1]:
                    if i == co:
                        return k[0]

        for z in re_ids:
            n = 0
            temp = [None] * int(len(z) - 1)
            ll_1, q = [], []
            ll_1.append([z[1][0][0], '5'])
            q.append(z[1][0][0])
            q.append('5')
            ll.append(ll_1)
            temp.remove(None)
            temp.insert(n, z[1][0])
            n += 1
            for t in temp:
                if None in temp:
                    for y in z[1:]:
                        if y[0] == t:
                            for u in y[1]:
                                if len(u) > 1:
                                    if u[0] in comp_dic:
                                        comp = comp_dic[u[0]] + u[1:]
                                        comp_li_id = look(z[1:], comp)
                                        if comp_li_id not in temp:
                                            orie_ind = q[q.index(y[0][0]) + 1]
                                            ex_orie = orie_dic[orie_ind]
                                            ll_1.append([comp_li_id[0], ex_orie])
                                            q.append(comp_li_id[0])
                                            q.append(ex_orie)
                                            temp.insert(n, comp_li_id)
                                            temp.remove(None)
                                            n += 1

        return ll

    key_find_lis = list(bind())

    # set orientation identifier to the list (complexes) as ['5', '-', '3']
    def set_orient():
        orie_li = []
        for i in re_ids:
            orie_li_1 = []
            orie_li_1.append(i[0])
            orie_li.append(orie_li_1)
            for ii in i[1:]:
                orie = ['5', '-', '3']
                orie_1 = [ii[0], orie, ii[1]]
                orie_li_1.append(orie_1)
        return orie_li

    set_orie = list(set_orient())

    # create a dictionary of items with it's orientation
    # ID set as dictionary key
    def ori_copy_dic():
        n_dic = []
        for cc in set_orie:
            n_dic_1 = {}
            n_dic.append(n_dic_1)
            for i in cc[1:]:
                orie_li_1 = [i[1:], i[2:]]
                re_dic = dict(zip(i[0], orie_li_1))
                n_dic_1.update(re_dic)
        return n_dic

    ori_list = list(ori_copy_dic())

    # set the list orientation ['5', '-', '3'] or ['3', '-', '5']
    # if ['3', '-', '5'] the list is written in reverse
    def comp_sort():
        ab = []
        for h, e in zip(key_find_lis, ori_list):
            li = [e[x[0]] for x in h]
            li_k = [h, li]
            ab.append(li_k)

        a_li = []
        for h in ab:
            a_s = []
            a_li.append(a_s)
            for e, c in zip(h[0], h[1]):
                num_i = [e, c]
                a_s.append(num_i)

        f_li = []
        for n in a_li:
            n_li = []
            f_li.append(n_li)
            for vv in n:
                if vv[0][1] == '3':
                    tw = vv[1][1][::-1]
                    tw_1 = vv[1][0][::-1]
                    tw_li = [[vv[0][0]], tw_1, tw]
                    n_li.append(tw_li)
                elif vv[0][1] == '5':
                    tw_li = [[vv[0][0]], vv[1][0], vv[1][1]]
                    n_li.append(tw_li)

        return f_li

    f_sorted_list = list(comp_sort())

    # re-write the compliments ID starting from 1, in ascending order
    def re_write():
        li_4, n_3 = [], []
        for a in f_sorted_list:
            li_3 = []
            li_4.append(li_3)
            q, n = [], 0

            for aa in a:
                li_2 = []
                li_3.append(li_2)
                li_2.append(aa[0])
                li_2.append(aa[1])
                for gg in aa[2:]:
                    li_1 = []
                    li_2.append(li_1)
                    for e in gg:
                        if len(e) > 1:
                            c_n = '|' + str(e[1:]) + '|'
                            if c_n not in q:
                                n += 1
                                let = e[0]
                                n_v = let + str(n)
                                li_1.append(n_v)
                                q.append(c_n)
                                q.append(str(n))

                            elif c_n in q:
                                let = e[0]
                                c_n_2 = q[q.index(c_n) + 1]
                                n_v = let + str(c_n_2)
                                li_1.append(n_v)
                        else:
                            li_1.append(e)
            n_3.append(n)

        return [li_4, n_3]

    f_list = list(re_write()[0])
    n_of_comps = list(re_write()[1])

    # check if same list as duplicates exists
    # if so recount them and sum their count
    def cal_same_comp():
        f_list_1 = list(copy.deepcopy(f_list))
        li_1 = []
        for c, noc in zip(f_list_1, n_of_comps):
            li_se_1 = []
            li_1.append([noc, li_se_1])
            for cc in c:
                li_se_1.append(cc[1:])

        n_li_1, n_amo, n_ssDNAs = [], [], []
        for uu in li_1:
            if uu not in n_li_1:
                n_li_1.append(uu)
                n_ssDNAs.append(len(uu[1]))
                i_ind_in = li_1.index(uu)
                n_amo.append(n_complexes[i_ind_in])
                li_1[i_ind_in] = ['x']

            elif uu in n_li_1:
                i_ind_in = n_li_1.index(uu)
                i_ind = li_1.index(uu)
                new_amo = n_complexes[i_ind]
                n_amo[i_ind_in] += new_amo
                li_1[i_ind] = ['x']

        return n_li_1, n_amo, n_ssDNAs

    comp_re_calculated = list(cal_same_comp())

    # this method is for creating and saving a bngl syntax for results obtained
    # will be saved at the given destination directory as a 'txt' file
    def make_bngl():
        n_o_c = []

        def convert(o, a, n_o):
            n_l = copy.deepcopy(n_o)
            n_r = n_l + 1

            if o == '5':
                w_com_1 = str('!' + str(a[0][1:]) if len(a[0]) != 1 else '')
                f_ele = 'N(b~{},5,3{},W{})'.format(a[0][0], '!' + str(n_r), w_com_1)

                bgl = []
                bgl.append(f_ele)
                n_l += 1
                n_r += 1

                for s in a[1:-1]:
                    w_com_2 = str('!' + str(s[1:]) if len(s) != 1 else '')
                    bn = 'N(b~{},5{},3{},W{})'.format(s[0], '!' + str(n_l), '!' + str(n_r), w_com_2)

                    bgl.append(bn)
                    n_l += 1
                    n_r += 1

                w_com_3 = str('!' + str(a[-1][1:]) if len(a[-1]) != 1 else '')
                l_ele = 'N(b~{},5{},3,W{})'.format(a[-1][0], '!' + str(n_l), w_com_3)

                bgl.append(l_ele)
                n_o_c.clear(), n_o_c.append(n_l)

                return '.'.join(bgl)

            elif o == '3':
                w_com_4 = str('!' + str(a[0][1:])if len(a[0]) != 1 else '')
                f_ele = 'N(b~{},5{},3,W{})'.format(a[0][0], '!' + str(n_r), w_com_4)

                bgl = []
                bgl.append(f_ele)
                n_l += 1
                n_r += 1

                for s in a[1:-1]:
                    w_com_5 = str('!' + str(s[1:])if len(s) != 1 else '')
                    bn = 'N(b~{},5{},3{},W{})'.format(s[0], '!' + str(n_r), '!' + str(n_l), w_com_5)

                    bgl.append(bn)
                    n_l += 1
                    n_r += 1

                w_com_6 = str('!' + str(a[-1][1:])if len(a[-1]) != 1 else '')
                l_ele = 'N(b~{},5,3{},W{})'.format(a[-1][0], '!' + str(n_l), w_com_6)

                bgl.append(l_ele)
                n_o_c.clear(), n_o_c.append(n_l)

                return '.'.join(bgl)

        comps_bngl, n = [], 0
        for a, b in zip(comp_re_calculated[0], comp_re_calculated[1]):
            n_o_c.append(a[0])
            no_num = ['.'.join([convert(x[0][0], x[1], n_o_c[-1]) for x in a[1]]), str(b)]
            comps_bngl.append(no_num)

        return comps_bngl

    return make_bngl()
