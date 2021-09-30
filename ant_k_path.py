# -*-coding:utf8-*-
#
# @autor:iniwaper@gmail.com
# ACO algorithm for K shortest path
#
# 研究问题：
# 针对北京市轨道交通的部分网络图，求解出O（起点站）-D（终点站）的K短路问题。本研究采用蚁群算法求解K短路的方案，
# 不仅考虑到蚁群算法的各种优越性，更在于其能很好解决该问题，并为诸多相关的问题研究提供一种可行的解决方案
# 考虑到所处理的问题的特殊性，为增强算法效率，本设计对地图进行了简化，但仍能保证正确的表述问题。
#
# 总体算法概述：
# 蚁群算法：
# 其是一种用来在图中寻找优化路径的机率型算法，一种模拟进化算法，具有一种新的模拟进化优化方法的有效性和应用价值。
# 生物学告诉我们，蚂蚁寻找食物的过程中，会向周围环境中释放一定的信息素，以吸引其他蚂蚁一起找食物，慢慢的越来越多的
# 蚂蚁找到食物，如果某个时刻的某条食物路线短，可以预见的是往返于这条路的蚂蚁就多，那么散发的激素也就多，最终越来越多
# 的蚂蚁会到这条线路上来，也就是我们期望的大部分蚂蚁最终重复着这条路线，也即最短路。那么我们在计算机世界里如何利用这样
# 一种智能呢？其实，大多数看似复杂的事物均由简单的规则组成，只是数量的不同而已。所以不必为此担心，来看看这些简单的规则：
# a)范围：蚂蚁观察到的范围是一个方格世界(square)，假定为3*3，即真实世界中的八个方向，那么蚂蚁移动的范围也就确定了
# b)环境：蚂蚁所处的环境被模拟为一个世界(world),其中有障碍，其他蚂蚁，信息素(食物和窝)，并且蚂蚁仅能感知其移动范围内
# 的信息，显然这符合客观事实，当然环境会以一定的速率是信息素挥发掉。
# c)觅食(窝)规则:find_food_nest_rule，在蚂蚁感知的范围内寻找是否有食物，有则直接找到，否则检查是否有信息素，有则比较信息素量，
# 并朝大的方向走。我们必须让蚂蚁以小概率犯错而选择非最大的点，这符合客观事实，不然蚂蚁就丧失了创新性，探索性，这不利于其发展，当然也
# 不利于算法的完备性。另找窝和食物规则基本一致，不再赘述
# d)移动规则:move_rule,每只蚂蚁都朝信息素最多的方向移动(这与觅食规则有所重复，本算法做了适当优化)，当周围没有信息素指引的时候，蚂蚁则
# 惯性的按原来的方向运动，并且在此方向上会有随机扰动，这可以理解为环境因素的作用，毕竟方向感总有不正确的时候。那么有一个问题，，蚂蚁很可能
# 会在原地打转，如何解决呢？我们采用简单的记忆策略，即蚂蚁会记录下最近走过的地方，不再重复走，而是尽量避开。
# e)避障规则:avoid_obstacle_rule，如果蚂蚁要移动的方向有障碍物挡住，那么它会随机选择一个方向，有信息素指引则按照觅食规则(这里也有重复行为)，如果
# 找不到可运动方向，则认为该蚂蚁走入死胡同(当然仅在本研究问题出现)。
# f)释放信息素规则：spread_pheromone_rule，不难想象，蚂蚁找到食物的时候释放的食物信息素最多，找到窝的时候释放窝信息素最多，并随着其走远，释放会减少。
# 综合以上规则(可以有更多规则，但超过本研究范围，不在讨论)，可以看到，蚂蚁之间并没有必然联系，而是与环境发生信息交互，正式通过这个客观想象，使得蚂蚁关联起来。
# 举个例子，当一只蚂蚁找到食物时，并没有也无法告知其他蚂蚁哪里有食物，而是向环境释放信息素，当其他蚂蚁靠近时，就会跟随过来，从而也找到了食物。
# 关于某些问题的解释(算法部分参数)：在没有任何信息素指引的情况下，蚂蚁是如何相对有效的找到食物的呢？答案是移动规则，这样保证了蚂蚁尽量向前，这就避免的原地打转
# 如果只是一味的向前，不太符合基本常理，所以算法使用PERTURBATION扰动参数来干扰蚂蚁的运动方向，这样就有了一定的目的性，保持原来的方向，又一定的试探，尤其是
# 遇到障碍的时候回立即改变方向。当然一只蚂蚁找到食物时，其他的会跟着找到食物。至于如何找到最短路的，上面已经介绍，归结于信息素，相同的时间内往返的蚂蚁多，释放的信息素
# 多，而且越来越多，趋于稳定。那么这样的规则下自然会产生局部最短和全局最短的问题，其实这个要归功于蚂蚁会选择错误的信息素点，这是一种创新，而一旦当这样创新找到的
# 路径更短时，那么经过充分的时间，这条路线上的蚂蚁一定会多起来，经过这样一系列的寻找，探索，最短路便被找到了。
# 本设计为增加扩展性，算法采用面向对象编程，方便扩展，包括以下几部分
# 设置类：主要用于所要求解问题的输入参数设置，以及蚁群算法的参数设置
# 蚂蚁类：主要用于表示蚂蚁个体，包括蚂蚁所具备的属性，以及特有生理现象
# 食物类：如果单针对本研究问题，食物类不是必须的，但考虑到扩展性可应用性，增加食物类，即世界中可以有不同种类的不同大小的食物，蚂蚁会做适当的选择，决定哪种更应先搬回窝。
# 蚁巢类：同样，蚁巢类不是必须的，但世界中可能有多个蚁巢，这是一个仿生环境，已不局限于所要解决所研究的问题
# 激素类：蚂蚁信息素也不是必须类，但为了更能表现蚂蚁世界的多样性，也考虑到有其他生物加入的可能性，激素类对增加算法的通用性是必须的
# 范围类:即蚂蚁移动的感知的范围，更确切的说与蚂蚁无关，它完全可以是其他生物的范围，所以为其不应携带蚂蚁信息。包括坐标信息，食物，窝，甚至可以包括敌窝等信息，它表征的是世界组成信息
# 世界类：显然这是模拟世界，里面不仅可以有蚂蚁，同样可以由各种生物，当然也可以理解为蚂蚁的世界，部落，世界之外还有世界等信息。其又很多范围组成。
#
# K短路求解算法：
# 对于K短路的求解，本设计采用了比较简单的策略，即在蚂蚁搜索食物的过程中，对于返回的路径进行记录，并保存当前以及以前找到的最短路径，在搜索的过程中实时的更新这个值，那么最终得到
# 基本可以保证是最短路径(经过有效的迭代次数以后)，这在上面已经论述过。在迭代的过程中，还做了以下操作来求解K短路：排除记录表中已经存在的路线，排除路径损耗大于当前最短路线10分钟
# 的路线，将得到的路线按顺序插入路径表，显然这个操作是非常高效的，最重要的一点是这样最终得到的表中的路径正是我们需要的路径顺序，最短，第2,3,4...短路。当然为了最终排除大于最短路径10分钟
# 的路线，还需做最后一步操作，即去除表中大于最短路10分钟的路线以及以后的所有路线，这个非常容易实现。这样经过上述步骤，我们就获得了所有K短路线了。
# 算法函数：加载地图，执行初始化操作，以及最重要的求K短路的算法
#
# 算法详细步骤：
# 1)算法开始
# 2)初始化各个变量参数，初始化蚂蚁，世界等，重要参数为记录路线表k_paths,最短路线best_route，最短损耗shortest_len
# 3)执行算法，对于所有OD对，执行3)
# 4)迭代。在最大迭代次数内，对于每一只蚂蚁执行蚁群算法的各种规则，直到找到食物，返回该路线，并转为回窝操作，规则基本和觅食一致，当找到窝时，返回路线
# (本问题不使用该路线，所有排除，但考虑到算法适用性，仍然给出，而相应的窝激素和食物激素权重设为1，即互不影响，显然当二者影响时，可能会得到更佳的效果，但同样会增加控制的难度)
# 5)对于3)中，返回的每条路线，执行K算法求解策略，记录路线，更新最短路线等。另对于那些死亡的蚂蚁，则置为死亡状态，不再处理其运动操作。
# 6)迭代完毕，获得路线表，输出所有K短路
# 7)算法结束
#
# 算法总结以及分析：
# 从以上的分析可见，由于蚁群算法具有一定的随机性，针对本问题则显示出了一定的局限性，由于站点的不可重复性和不可回头性，蚂蚁容易走入死胡同，对蚂蚁的各项参数要求较高。
# 而且迭代次数无法预知，蚁群算法非一定从最大值向最小值收敛，不过基本可以保证最短路径的有效性，第K短路则不能保证解的完备性，这与蚁群算法本身特性相关。不过总体来看，
# 本解决方案在很好的解决所研究问题的同时，增加了算法的适用性，为更多的应用领域做了准备。当然算法不可能做到完美，更多的性能有待探讨和发现，限于水平，不再对蚁群算法以及相关领域做深入探讨
#

import random
import math


class Settings(object):  # 设置类
    ROUTE_LINE_1 = ["b", "t", "n", "e", "v", "k", "u"]
    ROUTE_LINE_2 = ["d", "i", "j", "k", "l", "m", "f", "n"]  # 2号线站点表
    ROUTE_LINE_4 = ["a", "c", "d", "e", "f", "g", "h"]  # 4号线
    ROUTE_LINE_13 = ["d", "o", "p", "q", "r", "s", "j"]  # 13号线
    STATIONS_ROUTE_TIME = {
        "ac": 13,
        "ca": 13,
        "cd": 12,
        "dc": 12,
        "de": 10,
        "ed": 10,
        "ef": 2,
        "fe": 2,
        "fg": 6,
        "gf": 6,
        "gh": 5,
        "hg": 5,
        "di": 9,
        "id": 9,
        "ij": 4,
        "ji": 4,
        "jk": 7,
        "kj": 7,
        "kl": 2,
        "lk": 2,
        "lm": 6,
        "ml": 6,
        "mf": 4,
        "fm": 4,
        "fn": 4,
        "nf": 4,
        "nd": 7,
        "dn": 7,
        "do": 5,
        "od": 5,
        "op": 12,
        "po": 12,
        "pq": 7,
        "qp": 7,
        "qr": 6,
        "rq": 6,
        "rs": 13,
        "sr": 13,
        "sj": 6,
        "js": 6,
        "bt": 23,
        "tb": 23,
        "tn": 6,
        "nt": 6,
        "ne": 3,
        "en": 3,
        "ev": 6,
        "ve": 6,
        "vk": 5,
        "kv": 5,
        "ku": 13,
        "uk": 13,
    }  # 各站点间耗时

    # #######修OD结构,必须填写完整的OD################################
    O_D = {
        "ac": 1480,
        "ad": 2432,
        "ae": 1704,
        "af": 816,
        "ag": 928,
        "ah": 472,
        "ai": 968,
        "aj": 848,
        "ak": 1328,
        "al": 1456,
        "am": 512,
        "an": 272,
        "ao": 808,
        "ap": 744,
        "aq": 512,
        "ar": 536,
        "as": 424,
        "ab": 576,
        "at": 1696,
        "av": 664,
        "au": 432,
        "ca": 1056,
        "cd": 1320,
        "ce": 1432,
        "cf": 896,
        "cg": 1040,
        "ch": 672,
        "ci": 1488,
        "cj": 1608,
        "ck": 1904,
        "cl": 1040,
        "cm": 760,
        "cn": 848,
        "co": 864,
        "cp": 1072,
        "cq": 752,
        "cr": 960,
        "cs": 680,
        "cb": 184,
        "ct": 1104,
        "cv": 1704,
        "cu": 1080,
        "da": 1392,
        "dc": 1344,
        "de": 832,
        "df": 2520,
        "dg": 2432,
        "dh": 1712,
        "di": 2136,
        "dj": 2352,
        "dk": 2824,
        "dl": 1976,
        "dm": 1278,
        "dn": 960,
        "do": 1128,
        "dp": 1784,
        "dq": 1656,
        "dr": 1888,
        "ds": 1696,
        "db": 1104,
        "dt": 2144,
        "dv": 2432,
        "du": 2208,
        "ea": 2144,
        "ec": 2920,
        "ed": 984,
        "ef": 608,
        "eg": 1136,
        "eh": 712,
        "ei": 808,
        "ej": 1088,
        "ek": 1624,
        "el": 1872,
        "em": 1384,
        "en": 312,
        "eo": 968,
        "ep": 1048,
        "eq": 832,
        "er": 872,
        "es": 1040,
        "eb": 2496,
        "et": 3088,
        "ev": 1232,
        "eu": 3208,
        "fa": 472,
        "fc": 920,
        "fd": 2296,
        "fe": 128,
        "fg": 360,
        "fh": 448,
        "fi": 1144,
        "fj": 1008,
        "fk": 1928,
        "fl": 816,
        "fm": 360,
        "fn": 448,
        "fo": 992,
        "fp": 816,
        "fq": 224,
        "fr": 344,
        "fs": 464,
        "fb": 1336,
        "ft": 1544,
        "fv": 984,
        "fu": 1520,
        "ga": 624,
        "gc": 1008,
        "gd": 1472,
        "ge": 848,
        "gf": 464,
        "gh": 248,
        "gi": 1744,
        "gj": 840,
        "gk": 1184,
        "gl": 448,
        "gm": 360,
        "gn": 688,
        "go": 1264,
        "gp": 984,
        "gq": 848,
        "gr": 192,
        "gs": 824,
        "gb": 864,
        "gt": 1552,
        "gv": 1680,
        "gu": 1008,
        "ha": 752,
        "hc": 1264,
        "hd": 824,
        "he": 1168,
        "hf": 1080,
        "hg": 408,
        "hi": 1632,
        "hj": 992,
        "hk": 1736,
        "hl": 1016,
        "hm": 832,
        "hn": 1312,
        "ho": 992,
        "hp": 544,
        "hq": 672,
        "hr": 560,
        "hs": 864,
        "hb": 672,
        "ht": 1728,
        "hv": 1288,
        "hu": 1472,
        "ia": 432,
        "ic": 552,
        "id": 784,
        "ie": 1912,
        "if": 1920,
        "ig": 1272,
        "ih": 448,
        "ij": 328,
        "ik": 1392,
        "il": 544,
        "im": 1296,
        "in": 832,
        "io": 1136,
        "ip": 256,
        "iq": 464,
        "ir": 248,
        "is": 544,
        "ib": 288,
        "it": 944,
        "iv": 1248,
        "iu": 1968,
        "ja": 624,
        "jc": 1472,
        "jd": 1704,
        "je": 1344,
        "jf": 1448,
        "jg": 1824,
        "jh": 808,
        "ji": 328,
        "jk": 808,
        "jl": 568,
        "jm": 712,
        "jn": 1232,
        "jo": 1128,
        "jp": 168,
        "jq": 280,
        "jr": 216,
        "js": 512,
        "jb": 808,
        "jt": 944,
        "jv": 1248,
        "ju": 1792,
        "ka": 904,
        "kc": 3064,
        "kd": 2848,
        "ke": 944,
        "kf": 2192,
        "kg": 2408,
        "kh": 1248,
        "ki": 1880,
        "kj": 2016,
        "kl": 1264,
        "km": 1312,
        "kn": 2032,
        "ko": 1368,
        "kp": 2928,
        "kq": 1704,
        "kr": 1880,
        "ks": 1840,
        "kb": 912,
        "kt": 1968,
        "kv": 816,
        "ku": 2024,
        "la": 840,
        "lc": 1488,
        "ld": 2688,
        "le": 1536,
        "lf": 1160,
        "lg": 992,
        "lh": 1048,
        "li": 680,
        "lj": 184,
        "lk": 904,
        "lm": 448,
        "ln": 696,
        "lo": 720,
        "lp": 1248,
        "lq": 808,
        "lr": 1008,
        "ls": 848,
        "lb": 992,
        "lt": 2280,
        "lv": 960,
        "lu": 1680,
        "ma": 536,
        "mc": 1104,
        "md": 2776,
        "me": 992,
        "mf": 800,
        "mg": 1136,
        "mh": 696,
        "mi": 808,
        "mj": 984,
        "mk": 1472,
        "ml": 960,
        "mn": 1080,
        "mo": 904,
        "mp": 592,
        "mq": 856,
        "mr": 840,
        "ms": 808,
        "mb": 752,
        "mt": 1496,
        "mv": 1000,
        "mu": 1024,
        "na": 432,
        "nc": 1312,
        "nd": 1704,
        "ne": 328,
        "nf": 272,
        "ng": 1128,
        "nh": 1072,
        "ni": 1424,
        "nj": 1312,
        "nk": 2608,
        "nl": 1424,
        "nm": 1712,
        "no": 1072,
        "np": 808,
        "nq": 992,
        "nr": 504,
        "ns": 784,
        "nb": 1472,
        "nt": 992,
        "nv": 1440,
        "nu": 1728,
        "oa": 816,
        "oc": 904,
        "od": 1432,
        "oe": 1264,
        "of": 1136,
        "og": 1048,
        "oh": 368,
        "oi": 1232,
        "oj": 1088,
        "ok": 1304,
        "ol": 1968,
        "om": 1112,
        "on": 1040,
        "op": 216,
        "oq": 504,
        "or": 408,
        "os": 456,
        "ob": 1336,
        "ot": 2008,
        "ov": 1280,
        "ou": 992,
        "pa": 592,
        "pc": 368,
        "pd": 1128,
        "pe": 968,
        "pf": 1832,
        "pg": 808,
        "ph": 752,
        "pi": 1344,
        "pj": 1704,
        "pk": 1312,
        "pl": 1176,
        "pm": 304,
        "pn": 320,
        "po": 456,
        "pq": 272,
        "pr": 544,
        "ps": 1136,
        "pb": 280,
        "pt": 776,
        "pv": 816,
        "pu": 1360,
        "qa": 840,
        "qc": 904,
        "qd": 856,
        "qe": 512,
        "qf": 1088,
        "qg": 904,
        "qh": 696,
        "qi": 776,
        "qj": 856,
        "qk": 1768,
        "ql": 1264,
        "qm": 2240,
        "qn": 864,
        "qo": 1000,
        "qp": 688,
        "qr": 320,
        "qs": 800,
        "qb": 584,
        "qt": 1616,
        "qv": 912,
        "qu": 1312,
        "ra": 352,
        "rc": 392,
        "rd": 592,
        "re": 1728,
        "rf": 808,
        "rg": 768,
        "rh": 480,
        "ri": 864,
        "rj": 1160,
        "rk": 1528,
        "rl": 1256,
        "rm": 1312,
        "rn": 2048,
        "ro": 1424,
        "rp": 960,
        "rq": 456,
        "rs": 1072,
        "rb": 432,
        "rt": 1808,
        "rv": 832,
        "ru": 1392,
        "sa": 512,
        "sc": 952,
        "sd": 1112,
        "se": 1296,
        "sf": 1144,
        "sg": 864,
        "sh": 904,
        "si": 936,
        "sj": 648,
        "sk": 1608,
        "sl": 1472,
        "sm": 1344,
        "sn": 824,
        "so": 1256,
        "sp": 1192,
        "sq": 1056,
        "sr": 872,
        "sb": 320,
        "st": 1144,
        "sv": 968,
        "su": 1472,
        "ta": 912,
        "tc": 1584,
        "td": 2432,
        "te": 1126,
        "tf": 1312,
        "tg": 1096,
        "th": 992,
        "ti": 1432,
        "tj": 1736,
        "tk": 2576,
        "tl": 1312,
        "tm": 1192,
        "tn": 408,
        "to": 992,
        "tp": 280,
        "tq": 696,
        "tr": 752,
        "ts": 1016,
        "tb": 856,
        "tv": 992,
        "tu": 1192,
        "va": 448,
        "vc": 1216,
        "vd": 2048,
        "ve": 544,
        "vf": 1312,
        "vg": 1008,
        "vh": 344,
        "vi": 1176,
        "vj": 1528,
        "vk": 896,
        "vl": 832,
        "vm": 2800,
        "vn": 1472,
        "vo": 1120,
        "vp": 672,
        "vq": 2304,
        "vr": 936,
        "vs": 536,
        "vb": 608,
        "vt": 1664,
        "vu": 2512,
        "ua": 368,
        "uc": 1008,
        "ud": 3264,
        "ue": 2760,
        "uf": 952,
        "ug": 1560,
        "uh": 1656,
        "ui": 1432,
        "uj": 1496,
        "uk": 3472,
        "ul": 984,
        "um": 1576,
        "un": 1080,
        "uo": 1728,
        "up": 2048,
        "uq": 1960,
        "ur": 2512,
        "us": 2064,
        "ub": 368,
        "ut": 1048,
        "uv": 1048,
        "ba": 1192,
        "bc": 1496,
        "bd": 2136,
        "be": 1336,
        "bf": 1072,
        "bg": 968,
        "bh": 752,
        "bi": 1072,
        "bj": 1256,
        "bk": 3400,
        "bl": 1664,
        "bm": 1272,
        "bn": 1496,
        "bo": 776,
        "bp": 856,
        "bq": 1024,
        "br": 1168,
        "bs": 832,
        "bt": 1744,
        "bu": 1968,
        "bv": 1072,
    }
    # 需求解的OD对以及总人数
    # ###所有相邻的两点,这里必须填上所有的相邻点，初始值都是0############
    ALL_TWO = {
        "ac": 0,
        "cd": 0,
        "de": 0,
        "ef": 0,
        "fg": 0,
        "gh": 0,
        "di": 0,
        "ij": 0,
        "jk": 0,
        "kl": 0,
        "lm": 0,
        "mf": 0,
        "fn": 0,
        "nd": 0,
        "do": 0,
        "op": 0,
        "pq": 0,
        "qr": 0,
        "rs": 0,
        "sj": 0,
        "bt": 0,
        "tn": 0,
        "ne": 0,
        "ev": 0,
        "vk": 0,
        "ku": 0,
        "ca": 0,
        "dc": 0,
        "ed": 0,
        "fe": 0,
        "gf": 0,
        "hg": 0,
        "id": 0,
        "ji": 0,
        "kj": 0,
        "lk": 0,
        "ml": 0,
        "fm": 0,
        "nf": 0,
        "dn": 0,
        "od": 0,
        "po": 0,
        "qp": 0,
        "rq": 0,
        "sr": 0,
        "js": 0,
        "tb": 0,
        "nt": 0,
        "en": 0,
        "ve": 0,
        "kv": 0,
        "uk": 0,
    }

    D_CHANGE_ROUTE = {
        "4->13": 14.5,
        "4->2": 7.5,
        "13->2": 15,
        "13->4": 12.75,
        "2->4": 8.25,
        "2->13": 14.25,
    }  # D站点换乘耗时
    F_CHANGE_ROUTE = {"4->2": 7.5, "2->4": 9.75}  # F
    J_CHANGE_ROUTE = {"13->2": 9, "2->13": 11.25}  # J
    N_CHANGE_ROUTE = {"1->2": 7.5, "2->1": 7.5}  # J
    E_CHANGE_ROUTE = {"1->4": 7.5, "4->1": 6}  # J
    K_CHANGE_ROUTE = {"1->2": 9, "2->1": 7.5}  # J
    MAX_NC = 100  # 最大迭代次数
    ANTS_NUM = 50  # 蚂蚁个数
    PHEROMONE_WEIGHT = 1.0  # 表征信息素重要程度的参数(这里指食物和蚁巢素所占权重)
    MISTAKE_RATE = 1.0 / ANTS_NUM  # 表征启发式因子重要程度的参数
    PERTURBATION = 1.0  # 表征移动时的扰动参数,由于这里的线路限制，扰动为100%出现，可避免大量死蚂蚁
    MEMORY_ABILITY = 8  # 表征蚂蚁记忆力参数
    RHO = 1.0 / MAX_NC  # 信息素蒸发系数
    PHEROMONE = 100  # 信息素总量参数


class Ant(object):  # 蚂蚁类
    def __init__(self, world, square, ant_id):
        self.world = world  # 蚂蚁所处的世界
        self.square = square  # 蚂蚁所在方格
        self.behavior = "FOOD"  # 当前蚂蚁行为，初始化为觅食
        self.prev_direction = -1  # 0-7分别为八个方向
        self.total_food = 0  # 所搬的食物数
        self.ant_id = ant_id  # 蚂蚁的编号
        self.pos_memory = []  # 蚂蚁记忆最近走过的路
        self.total_route_len = 1  # 总路线长度
        self.food_route = []  # 觅食路线
        self.nest_route = []  # 寻巢路线
        self.route = []  # 路线
        self.is_dead = False  # 是否已经进入死角
        self.scope = self.get_scope(square)  # 蚂蚁当前能观察到的范围

    # ########以下为扩展############
    def get_square(self):
        return self.square

    def set_square(self, square):
        self.square = square

    def get_total_food(self):
        return self.total_food

    def set_total_food(self):
        pass

    def get_ant_id(self):
        return self.ant_id

    def set_ant_id(self, ant_id):
        self.ant_id = ant_id

    # #####扩展结束###############
    def get_scope(self, square):  # 获取当前蚂蚁的【范围】
        x, y = square.x, square.y
        self.route.append(square.label)  #
        self.food_route.append(square)  #
        self.pos_memory.append(square)  #
        scope = [
            self.world.SQUARES[y][x - 1],
            self.world.SQUARES[y + 1][x - 1],
            self.world.SQUARES[y + 1][x],
            self.world.SQUARES[y + 1][x + 1],
            self.world.SQUARES[y][x + 1],
            self.world.SQUARES[y - 1][x + 1],
            self.world.SQUARES[y - 1][x],
            self.world.SQUARES[y - 1][x - 1],
        ]
        return scope

    # 蚂蚁所具备的本能【前进->觅食或寻窝】
    def go_next_square(self):
        next_direction = -1
        square = None
        route_line = ""

        def change_route_cost():  # 计算换乘代价
            change_cost = 0
            # 获得当前换乘情况
            if self.prev_direction >= 4:
                tmp = self.prev_direction - 4
            else:
                tmp = self.prev_direction + 4
            prev_station = self.scope[tmp].label
            cur_station = self.square.label
            next_station = square.label

            if cur_station == "d":
                if prev_station in Settings.ROUTE_LINE_4 and next_station in Settings.ROUTE_LINE_13:
                    change_cost = Settings.D_CHANGE_ROUTE["4->13"]
                elif (
                    prev_station in Settings.ROUTE_LINE_4 and next_station in Settings.ROUTE_LINE_2
                ):
                    change_cost = Settings.D_CHANGE_ROUTE["4->2"]
                elif (
                    prev_station in Settings.ROUTE_LINE_13 and next_station in Settings.ROUTE_LINE_2
                ):
                    change_cost = Settings.D_CHANGE_ROUTE["13->2"]
                elif (
                    prev_station in Settings.ROUTE_LINE_13 and next_station in Settings.ROUTE_LINE_4
                ):
                    change_cost = Settings.D_CHANGE_ROUTE["13->4"]
                elif (
                    prev_station in Settings.ROUTE_LINE_2 and next_station in Settings.ROUTE_LINE_4
                ):
                    change_cost = Settings.D_CHANGE_ROUTE["2->4"]
                elif (
                    prev_station in Settings.ROUTE_LINE_2 and next_station in Settings.ROUTE_LINE_13
                ):
                    change_cost = Settings.D_CHANGE_ROUTE["2->13"]
            elif cur_station == "f":
                if prev_station in Settings.ROUTE_LINE_4 and next_station in Settings.ROUTE_LINE_2:
                    change_cost = Settings.F_CHANGE_ROUTE["4->2"]
                elif (
                    prev_station in Settings.ROUTE_LINE_2 and next_station in Settings.ROUTE_LINE_4
                ):
                    change_cost = Settings.F_CHANGE_ROUTE["2->4"]
            elif cur_station == "j":
                if prev_station in Settings.ROUTE_LINE_13 and next_station in Settings.ROUTE_LINE_2:
                    change_cost = Settings.J_CHANGE_ROUTE["13->2"]
                elif (
                    prev_station in Settings.ROUTE_LINE_2 and next_station in Settings.ROUTE_LINE_13
                ):
                    change_cost = Settings.J_CHANGE_ROUTE["2->13"]
            elif cur_station == "n":
                if prev_station in Settings.ROUTE_LINE_1 and next_station in Settings.ROUTE_LINE_2:
                    change_cost = Settings.N_CHANGE_ROUTE["1->2"]
                elif (
                    prev_station in Settings.ROUTE_LINE_2 and next_station in Settings.ROUTE_LINE_1
                ):
                    change_cost = Settings.N_CHANGE_ROUTE["2->1"]
            elif cur_station == "e":
                if prev_station in Settings.ROUTE_LINE_1 and next_station in Settings.ROUTE_LINE_4:
                    change_cost = Settings.E_CHANGE_ROUTE["1->4"]
                elif (
                    prev_station in Settings.ROUTE_LINE_4 and next_station in Settings.ROUTE_LINE_1
                ):
                    change_cost = Settings.E_CHANGE_ROUTE["4->1"]
            elif cur_station == "k":
                if prev_station in Settings.ROUTE_LINE_1 and next_station in Settings.ROUTE_LINE_2:
                    change_cost = Settings.K_CHANGE_ROUTE["1->2"]
                elif (
                    prev_station in Settings.ROUTE_LINE_2 and next_station in Settings.ROUTE_LINE_1
                ):
                    change_cost = Settings.K_CHANGE_ROUTE["2->1"]

            return change_cost

        if self.is_dead:  # 如果蚂蚁状态为死亡，则不再移动
            return (None, 10000)

        (next_direction, square, route_line, find_food_nest) = self.find_food_nest_rule()  # 执行觅食规则
        if find_food_nest:  # 找到食物或者窝
            if self.behavior == "FOOD":
                self.food_route.append(square)
                self.route.append(square.label)
                route = self.route

                self.behavior = "NEST"
                self.nest_route = []
                self.food_route = []
                self.nest_route.append(square)
            else:
                self.nest_route.append(square)
                self.route.append(square.label)
                route = self.route

                self.behavior = "FOOD"
                self.total_food += 1
                self.nest_route = []
                self.food_route = []
                self.food_route.append(square)
            # 计算最终路线耗时
            self.total_route_len += Settings.STATIONS_ROUTE_TIME[route_line]
            if (
                self.square.label == "d"
                or self.square.label == "f"
                or self.square.label == "j"
                or self.square.label == "n"
                or self.square.label == "e"
                or self.square.label == "k"
            ) and self.prev_direction != -1:
                self.total_route_len += change_route_cost()
            total_route_len = self.total_route_len - 1
            # 置位所有参数，进行相反的行为(觅食寻窝转换)
            self.square = square
            self.prev_direction = -1
            self.total_route_len = 1
            self.pos_memory = []
            self.route = []

            self.pos_memory.append(square)  # 蚂蚁记忆最近走过的路
            self.route.append(square.label)  # 路线
            x, y = self.square.x, self.square.y
            self.scope = [
                self.world.SQUARES[y][x - 1],
                self.world.SQUARES[y + 1][x - 1],
                self.world.SQUARES[y + 1][x],
                self.world.SQUARES[y + 1][x + 1],
                self.world.SQUARES[y][x + 1],
                self.world.SQUARES[y - 1][x + 1],
                self.world.SQUARES[y - 1][x],
                self.world.SQUARES[y - 1][x - 1],
            ]
            self.spread_pheromone_rule(self.square.x, self.square.y)  # 执行播撒信息素规则

            return (route, total_route_len)  # 返回获得的路线
        elif square is None:  # 如果觅食规则无效，则执行移动规则(无效即无信息素指引)
            (next_direction, square, route_line) = self.move_rule()  # 执行移动规则
            if square is None:  # 如果移动规则无效(即撞墙)，则执行避障规则
                (next_direction, square, route_line) = self.avoid_obstacle_rule()  # 执行避障规则
        if self.is_dead:  # 检查避障规则后蚂蚁是否死亡，是则直接返回
            return (None, 10000)
        # 计算当前移动后，路径耗时
        self.total_route_len += Settings.STATIONS_ROUTE_TIME[route_line]
        if (
            self.square.label == "d"
            or self.square.label == "f"
            or self.square.label == "j"
            or self.square.label == "n"
            or self.square.label == "e"
            or self.square.label == "k"
        ) and self.prev_direction != -1:
            self.total_route_len += change_route_cost()
        # 执行移动
        self.square = square
        self.prev_direction = next_direction
        x, y = self.square.x, self.square.y
        self.scope = [
            self.world.SQUARES[y][x - 1],
            self.world.SQUARES[y + 1][x - 1],
            self.world.SQUARES[y + 1][x],
            self.world.SQUARES[y + 1][x + 1],
            self.world.SQUARES[y][x + 1],
            self.world.SQUARES[y - 1][x + 1],
            self.world.SQUARES[y - 1][x],
            self.world.SQUARES[y - 1][x - 1],
        ]
        memory_len = len(self.pos_memory) - 1  # 检查当前已记忆值
        if memory_len >= Settings.MEMORY_ABILITY:
            self.pos_memory[Settings.MEMORY_ABILITY - 1] = self.square
        else:
            self.pos_memory.append(self.square)
        # 根据当前行为记录当前路径到对应表
        if self.behavior == "FOOD":
            self.food_route.append(self.square)
        if self.behavior == "NEST":
            self.nest_route.append(self.square)
        # 记录该位置到路径表
        self.route.append(self.square.label)
        self.spread_pheromone_rule(self.square.x, self.square.y)  # 执行播撒信息素规则
        return (None, 10000)

    def move_rule(self):
        """
        移动规则：每只蚂蚁都朝向外激素最多的方向移，并且，当周围没有外激素指引的时候，蚂蚁会按照自己原来
        运动的方向惯性的运动下去，并且，在运动的方向有一个随机的小的扰动。为了防止蚂蚁原地转圈，
        它会记住最近刚走过了哪些点，如果发现要走的下一点已经在最近走过了，它就会尽量避开。
        """
        station_routes = Settings.STATIONS_ROUTE_TIME.keys()
        routes = []
        square = None
        next_direction = -1
        line = ""

        def search_route_line():  # 查找上次所在的线路，获得惯性方向
            next_station = None
            route_lines = [
                Settings.ROUTE_LINE_1,
                Settings.ROUTE_LINE_2,
                Settings.ROUTE_LINE_4,
                Settings.ROUTE_LINE_13,
            ]
            next_line = None

            if self.prev_direction >= 4:
                tmp = self.prev_direction - 4
            else:
                tmp = self.prev_direction + 4
            for i in range(4):
                if self.scope[tmp].label in route_lines[i] and self.square.label in route_lines[i]:
                    next_line = i  # 哪一条线路

                    index_p = route_lines[i].index(self.scope[tmp].label)
                    index_c = route_lines[i].index(self.square.label)
                    if index_p == 0 and index_c == len(route_lines[i]) - 1:
                        # 0-->last-->
                        next_station = index_c - 1
                    elif index_c == 0 and index_p == len(route_lines[i]) - 1:
                        # last-->0-->
                        next_station = index_c + 1
                    elif index_p > index_c:
                        # <----
                        if index_c == 0:
                            next_station = len(route_lines[i]) - 1
                        else:
                            next_station = index_c - 1
                    elif index_p < index_c:
                        # --->
                        if index_c == len(route_lines[i]) - 1:
                            next_station = 0
                        else:
                            next_station = index_c + 1
                    break

            for i in range(len(routes)):
                if routes[i].label == route_lines[next_line][next_station]:
                    return (self.scope.index(routes[i]), routes[i])
            return (-1, None)

        for i in range(8):  # 搜索所有连通当前站点且未走过的站点
            if (
                self.scope[i] is None
                or self.scope[i] in self.food_route
                or self.scope[i] in self.nest_route
            ):
                continue
            if self.square.label + self.scope[i].label in station_routes:
                routes.append(self.scope[i])
        routes_num = len(routes)
        if routes_num == 0:  # 没有可选线路，直接返回
            return (-1, None, "")
        if self.prev_direction != -1:
            (next_direction, square) = search_route_line()
        else:
            square = random.choice(routes)
            next_direction = self.scope.index(square)
        if random.random() <= Settings.PERTURBATION:  # 扰动规则，符合则进行随机扰动
            mistake = random.randint(0, 7)
            square = self.scope[mistake]
            for i in range(8):
                if square is not None and square not in routes:
                    square = None
            next_direction = mistake
        if square is None:
            line = ""
        else:
            line = self.square.label + square.label
        return (next_direction, square, line)

    def find_food_nest_rule(self):
        """
        觅食寻巢规则：在每只蚂蚁能感知的范围内寻找是否有食物，如果有就直接过去。否则看是否有外激素，
        并且比较在能感知的范围内哪一点的外激素最多，这样，它就朝外激素多的地方走，并且每只蚂蚁
        多会以小概率犯错误，从而并不是往外激素最多的点移动。蚂蚁找窝的规则和上面一样，只不过
        它对窝的外激素做出反应，而对食物外激素没反应。
        """

        def calculate_pheromone(square, behavior):  # 计算信息素量
            if behavior == "FOOD":
                return (
                    square.pheromone.food_pheromone * Settings.PHEROMONE_WEIGHT
                    + square.pheromone.nest_pheromone * (1 - Settings.PHEROMONE_WEIGHT)
                )
            if behavior == "NEST":
                return (
                    square.pheromone.food_pheromone * (1 - Settings.PHEROMONE_WEIGHT)
                    + square.pheromone.nest_pheromone * Settings.PHEROMONE_WEIGHT
                )

        station_routes = Settings.STATIONS_ROUTE_TIME.keys()
        routes = []
        bool_pheromone = False
        square = None
        next_direction = -1

        for i in range(8):  # 搜索所有连通当前站点且未走过的站点
            if (
                self.scope[i] is None
                or self.scope[i] in self.food_route
                or self.scope[i] in self.nest_route
            ):
                continue
            if self.square.label + self.scope[i].label in station_routes:
                routes.append(self.scope[i])

        routes_num = len(routes)
        if routes_num == 0:  # 没有可选线路，直接返回
            return (-1, None, "", False)

        for i in range(routes_num):
            if (
                routes[i].square_food.food_type is not None and self.behavior == "FOOD"
            ):  # 如果有食物，直接返回，不再考查激素
                return (i, routes[i], self.square.label + routes[i].label, True)
            if (
                routes[i].square_nest.nest_type is not None and self.behavior == "NEST"
            ):  # 如果有蚁巢，直接返回，不再考查激素
                return (i, routes[i], self.square.label + routes[i].label, True)  # 方向，方格，线路，是否找到
            if (routes[i].pheromone.food_pheromone != 0 and self.behavior == "FOOD") or (
                routes[i].pheromone.nest_pheromone != 0 and self.behavior == "NEST"
            ):
                bool_pheromone = True

        if not bool_pheromone:  # 如果没有激素，则返回执行后面的规则
            return (-1, None, "", False)
        else:  # 有激素，则获取所有可选方向
            pheromone_values = []
            pheromone_current_best_value = 0
            pheromone_value = 0
            line = ""
            for i in range(routes_num):
                pheromone_value = calculate_pheromone(routes[i], self.behavior)  # 计算方格值
                pheromone_values.append((pheromone_value, i))
                if pheromone_value > pheromone_current_best_value:
                    pheromone_current_best_value = pheromone_value
                    square = routes[i]
                    next_direction = self.scope.index(square)
                    next_index = i
                    line = self.square.label + square.label

            if len(pheromone_values) >= 2:  # 少于2个，不具有选错性
                if random.random() <= Settings.MISTAKE_RATE:  # 执行错误选择信息素方向规则，满足则去除最大值，从剩下的随机选取一个移动
                    values = [v[1] for v in pheromone_values if v[1] != next_index]
                    next_direction = random.choice(values)
                    (next_direction, square) = (
                        self.scope.index(routes[next_direction]),
                        routes[next_direction],
                    )
                    line = self.square.label + square.label

            return (next_direction, square, line, False)

    def avoid_obstacle_rule(self):
        """
        避障规则：如果蚂蚁要移动的方向有障碍物挡住，它会随机的选择另一个方向，并且有外激素指引的话，
        它会按照觅食的规则行为
        """
        station_routes = Settings.STATIONS_ROUTE_TIME.keys()
        routes = []

        for i in range(8):  # 搜索所有连通当前站点且未走过的站点
            if (
                self.scope[i] is None
                or self.scope[i] in self.food_route
                or self.scope[i] in self.nest_route
            ):
                continue
            if self.square.label + self.scope[i].label in station_routes:
                routes.append((i, self.scope[i], self.square.label + self.scope[i].label))
        routes_num = len(routes)
        if routes_num == 0:  # 没有可选线路，直接返回，设蚂蚁为死亡
            self.is_dead = True
            return (-1, None, "")
        else:
            return random.choice(routes)

    def spread_pheromone_rule(self, x, y):
        """
        播撒外激素规则：每只蚂蚁在刚找到食物或者窝的时候撒发的外激素最多，并随着它走远的距离，
        播撒的外激素越来越少。
        """
        # 释放信息素采用线性规则，即蚂蚁携带的激素总量除以当前路径损耗，可满足释放激素规则
        if self.behavior == "NEST":
            self.world.SQUARES[y][x].pheromone.food_pheromone += (
                1.0 * Settings.PHEROMONE / self.total_route_len
            )
        if self.behavior == "FOOD":
            self.world.SQUARES[y][x].pheromone.nest_pheromone += (
                1.0 * Settings.PHEROMONE / self.total_route_len
            )


# ####扩展类开始#############
class Food(object):  # 食物类
    def __init__(self, food_size, food_type):
        self.food_size = food_size
        self.food_type = food_type

    def get_food_size(self):
        return self.food_size

    def set_food_size(self, food_size):
        self.food_size = food_size

    def get_food_type(self):
        return self.food_type

    def set_food_type(self, food_type):
        self.food_type = food_type


class Nest(object):  # 蚁巢类
    def __init__(self, nest_type):
        self.nest_type = nest_type

    def get_nest_type(self):
        return self.nest_type

    def set_nest_type(self, nest_type):
        self.nest_type = nest_type


class Pheromone(object):  # 激素类
    def __init__(self, food_pheromone, nest_pheromone):
        self.food_pheromone = food_pheromone
        self.nest_pheromone = nest_pheromone

    def get_food_pheromone(self):
        return self.food_pheromone

    def set_food_pheromone(self, food_pheromone):
        self.food_pheromone = food_pheromone

    def get_nest_pheromone(self):
        return self.nest_pheromone

    def set_nest_pheromone(self, nest_pheromone):
        self.nest_pheromone = nest_pheromone


# ###########扩展类结束################


class Square(object):  # 范围类
    def __init__(self, x, y, label):
        self.square_food = Food(None, None)
        self.square_nest = Nest(None)
        self.pheromone = Pheromone(0.0, 0.0)
        self.x = x
        self.y = y
        self.label = label


class World(object):  # 世界类
    def __init__(self, world_map):
        self.food_pos = []
        self.nest_pos = []
        self.SQUARES = self.get_world_map_squares(world_map)

    def get_world_map_squares(self, world_map):  # 初始化世界地图
        width = len(world_map[0])
        height = len(world_map)
        SQUARES = [[None for x in range(width)] for y in range(height)]
        for y in range(height):
            for x in range(width):
                if world_map[y][x] == "1":
                    continue
                else:
                    SQUARES[y][x] = Square(x, y, world_map[y][x])
        return SQUARES

    def set_nest(self, label, *position):  # 设置窝点
        # (x,y)=position
        # self.SQUARES[x][y].square_nest.nest_type='NEST'
        #
        for y in range(len(self.SQUARES)):
            for x in range(len(self.SQUARES[0])):
                if self.SQUARES[y][x] is None:
                    continue
                if self.SQUARES[y][x].label == label:
                    self.SQUARES[y][x].square_nest.nest_type = "NEST"
                    self.nest_pos = [x, y]
                    return (x, y)

    def set_food(self, label, *position):  # 设置食物点
        # (x,y)=position
        # self.SQUARES[x][y].square_food.food_type='FOOD'
        # self.SQUARES[x][y].square_food.food_size='BIG'
        for y in range(len(self.SQUARES)):
            for x in range(len(self.SQUARES[0])):
                if self.SQUARES[y][x] is None:
                    continue
                if self.SQUARES[y][x].label == label:
                    self.SQUARES[y][x].square_food.food_type = "FOOD"
                    self.SQUARES[y][x].square_food.food_size = "BIG"
                    self.food_pos = [x, y]
                    return (x, y)

    def update_pheromone(self):  # 更新信息素规则(环境挥发)
        # 更新信息素,挥发
        for y in range(len(self.SQUARES)):
            for x in range(len(self.SQUARES[0])):
                if self.SQUARES[y][x] is not None:
                    self.SQUARES[y][x].pheromone.food_pheromone = (1 - Settings.RHO) * self.SQUARES[
                        y
                    ][x].pheromone.food_pheromone
                    self.SQUARES[y][x].pheromone.nest_pheromone = (1 - Settings.RHO) * self.SQUARES[
                        y
                    ][x].pheromone.nest_pheromone


def Init_ACO_K_ShortRoute():  # 算法函数
    # 加载地图
    world_map = []
    with open("map.txt") as f:
        for line in f:
            world_map.append("".join(line.split("\n")).split(","))
    # 迭代求解
    for od in Settings.O_D.keys():
        best_route = []  # 最短路径
        shortest_len = 10000  # 路径损耗(耗时)
        k_paths = []  # K短路径

        world = World(world_map)
        (x, y) = world.set_nest(od[0])
        world.set_food(od[1])
        ants = [Ant(world, world.SQUARES[y][x], i) for i in range(Settings.ANTS_NUM)]
        # 迭代
        for nc in range(Settings.MAX_NC):
            # 每只蚂蚁
            for ant_id in range(Settings.ANTS_NUM):
                # 觅食或找窝(移动)
                (route, route_len) = ants[ant_id].go_next_square()
                # 检查返回路线有效性，无效则返回
                if route is None or route[0] != od[0]:
                    continue
                # 返回路线是否大于当前最短10分钟，是则执行下轮操作
                if route_len - shortest_len > 10:
                    continue
                # 记录更新最短路线
                if route_len < shortest_len:
                    shortest_len = route_len
                    best_route = route
                # 第一条K短路线，则直接插入K短路线表，并执行下轮操作
                if len(k_paths) == 0:
                    k_paths.insert(0, (route, route_len))
                    continue
                # 如果该路线已经存在于K短路线表，则忽略，继续下轮操作
                if route in [k[0] for k in k_paths]:
                    continue
                else:
                    # 检查是否超过三次换乘
                    total_change = 0
                    for s in ["d", "f", "j", "n", "e", "k"]:
                        for r in route:
                            if s == r:
                                total_change += 1
                    if total_change > 3:
                        continue
                    # 查找路线应该插入K短路线表中的位置，如果到达K短路线表末尾，则直接加到末尾
                    tmp = len(k_paths)
                    for i in range(tmp):
                        if k_paths[i][1] >= route_len:
                            k_paths.insert(i, (route, route_len))
                            break
                        elif i == tmp - 1:
                            k_paths.append((route, route_len))
            # 更新信息素
            world.update_pheromone()
        # ###############迭代结束##############
        # 输出路线
        k_path_od = []
        if len(best_route) == 0:
            print("search K-path for OD:" + od + " fail,please run again any way.")
        else:
            for k_path in range(len(k_paths)):
                if k_paths[k_path][1] - shortest_len > 10:
                    break
                else:
                    print(
                        "The "
                        + str(k_path + 1)
                        + " K-path for OD: "
                        + od
                        + "->"
                        + "-".join(k_paths[k_path][0])
                        + " time:"
                        + str(k_paths[k_path][1])
                        + " min"
                    )
                    k_path_od.append((k_paths[k_path], Sx(k_paths[k_path][1], shortest_len)))

        # 断面流量处理
        # 只有一条路
        if len(k_path_od) == 1:
            # print("只有一条路")
            for i in range(len(k_path_od[0][0][0]) - 1):
                Settings.ALL_TWO[k_path_od[0][0][0][i] + k_path_od[0][0][0][i + 1]] += Settings.O_D[
                    od
                ]
        else:
            for k_path in range(len(k_path_od)):
                k_two = [k_path_od[i][0][1] for i in range(len(k_path_od))]
                k_two.insert(0, k_path_od[k_path][0][1])
                for i in range(len(k_path_od[k_path][0][0]) - 1):
                    Settings.ALL_TWO[
                        k_path_od[k_path][0][0][i] + k_path_od[k_path][0][0][i + 1]
                    ] += Settings.O_D[od] * Pk(k_two)


# 配送相关处理，求断面流量
def Sx(Ck, Cmin):
    return math.exp(-2 * (pow((Ck - Cmin), 2) / 25.0))


def Pk(Sk):
    return 1.0 * Sk[0] / (sum(Sk[1:]))


# 启动算法
if __name__ == "__main__":
    print("please waiting...")
    Init_ACO_K_ShortRoute()
    al = Settings.ALL_TWO.keys()
    for i in al:
        if int(Settings.ALL_TWO[i]) - Settings.ALL_TWO[i]:
            print(i, int(Settings.ALL_TWO[i]) + 1)
