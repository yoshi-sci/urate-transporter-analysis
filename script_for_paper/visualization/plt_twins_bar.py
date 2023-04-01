import matplotlib.pyplot as plt


def plt_twins_bar(database, option, reverse):
    plt.rcParams['font.family'] = 'Arial'

    fig = plt.figure(figsize=(4, 6), dpi=360)
    ax = fig.add_subplot(111, xlabel=database[0][0])
    if option == 1 or option == 3:
        plt.ylim(0, 100)

    l_label = []
    l_ITS = []
    l_ITS_DP = []
    l_ITS_DP_ITR = []
    l_ITS_DP_ITR_DN = []

    r_label = []
    r_ITS = []
    r_ITS_DP = []
    r_ITS_DP_ITR = []
    r_ITS_DP_ITR_DN = []

    # create data
    label = ['A1 A2', 'B1 B2', 'C1 C2', 'D1 D2']
    color = ['#5298FD', '#A16DFF', '#FD68A3', '#444444']
    patterns = ["/", "\\", "|", "-", "+", "x", "o", "O", ".", "*"]
    zorder = [4, 3, 2, 1]

    i = 0
    for data in database[1:]:
        # calculate total number divided by 100
        total = 0
        for str_num in data[1:]:
            total += int(str_num)
        if option == 1:
            total_dev_100 = total / 100
        elif option == 2:
            total_dev_100 = 1
        else:
            if i % 2 == 0:
                for str_num in database[i + 2][1:]:
                    total += int(str_num)
            else:
                for str_num in database[i][1:]:
                    total += int(str_num)
            total_dev_100 = total / 100

        # create data to plot
        if i % 2 == 0:
            l_label.append(data[0])
            l_ITS.append(int(data[1]) / total_dev_100)
            l_ITS_DP.append(int(data[2]) / total_dev_100 + l_ITS[-1])
            l_ITS_DP_ITR.append(int(data[3]) / total_dev_100 + l_ITS_DP[-1])
            l_ITS_DP_ITR_DN.append(int(data[4]) / total_dev_100 + l_ITS_DP_ITR[-1])
        else:
            r_label.append(data[0])
            r_ITS.append(int(data[1]) / total_dev_100)
            r_ITS_DP.append(int(data[2]) / total_dev_100 + r_ITS[-1])
            r_ITS_DP_ITR.append(int(data[3]) / total_dev_100 + r_ITS_DP[-1])
            r_ITS_DP_ITR_DN.append(int(data[4]) / total_dev_100 + r_ITS_DP_ITR[-1])
        i += 1

    if reverse:
        l_ITS, l_ITS_DP, l_ITS_DP_ITR, l_ITS_DP_ITR_DN, r_ITS, r_ITS_DP, r_ITS_DP_ITR, r_ITS_DP_ITR_DN, zorder = reverse_stack_bar_data(
            l_ITS, l_ITS_DP, l_ITS_DP_ITR, l_ITS_DP_ITR_DN, r_ITS, r_ITS_DP, r_ITS_DP_ITR, r_ITS_DP_ITR_DN, zorder)

    # plot bars in stack manner
    plot_twins_stack_bar(ax, color, database, fig, l_ITS, l_ITS_DP, l_ITS_DP_ITR, l_ITS_DP_ITR_DN, label, option,
                         patterns, r_ITS, r_ITS_DP, r_ITS_DP_ITR, r_ITS_DP_ITR_DN, zorder)


def plot_twins_stack_bar(ax, color, database, fig, l_ITS, l_ITS_DP, l_ITS_DP_ITR, l_ITS_DP_ITR_DN, label, option,
                         patterns, r_ITS, r_ITS_DP, r_ITS_DP_ITR, r_ITS_DP_ITR_DN, zorder):
    ax.bar(label, l_ITS, bottom=0, color=color[0], align='edge', width=-0.35, edgecolor='black', linewidth=1,
           zorder=zorder[0], hatch=patterns[0] * 2)
    ax.bar(label, l_ITS_DP, bottom=0, color=color[1], align='edge', width=-0.35, edgecolor='black', linewidth=1,
           zorder=zorder[1], hatch=patterns[9] * 2)
    ax.bar(label, l_ITS_DP_ITR, bottom=0, color=color[2], align='edge', width=-0.35, edgecolor='black', linewidth=1,
           zorder=zorder[2], hatch=patterns[4] * 2)
    ax.bar(label, l_ITS_DP_ITR_DN, bottom=0, color=color[3], align='edge', width=-0.35, edgecolor='black', linewidth=1,
           zorder=zorder[3], hatch=patterns[8] * 2)
    ax.bar(label, r_ITS, bottom=0, color=color[0], align='edge', width=0.35, edgecolor='black', linewidth=1,
           zorder=zorder[0], hatch=patterns[0] * 2)
    ax.bar(label, r_ITS_DP, bottom=0, color=color[1], align='edge', width=0.35, edgecolor='black', linewidth=1,
           zorder=zorder[1], hatch=patterns[9] * 2)
    ax.bar(label, r_ITS_DP_ITR, bottom=0, color=color[2], align='edge', width=0.35, edgecolor='black', linewidth=1,
           zorder=zorder[2], hatch=patterns[4] * 2)
    ax.bar(label, r_ITS_DP_ITR_DN, bottom=0, color=color[3], align='edge', width=0.35, edgecolor='black', linewidth=1,
           zorder=zorder[3], hatch=patterns[8] * 2)
    # setting legend here
    ax.legend(["ITS", "DP", "ITR", "DN"], bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=12)
    plt.subplots_adjust(right=0.7)
    fig.savefig('output/type' + str(option) + '/' + database[0][0],
                facecolor=fig.get_facecolor(), edgecolor=fig.get_edgecolor())
    print('figure type{}/{} saved successfully...'.format(str(option), database[0][0]))


def reverse_stack_bar_data(l_ITS, l_ITS_DP, l_ITS_DP_ITR, l_ITS_DP_ITR_DN, r_ITS, r_ITS_DP, r_ITS_DP_ITR,
                           r_ITS_DP_ITR_DN, zorder):
    l_DN = list(map(lambda x, y: x - y, l_ITS_DP_ITR_DN, l_ITS_DP_ITR))
    l_DN_ITR = list(map(lambda x, y: x - y, l_ITS_DP_ITR_DN, l_ITS_DP))
    l_DN_ITR_DP = list(map(lambda x, y: x - y, l_ITS_DP_ITR_DN, l_ITS))
    l_DN_ITR_DP_ITS = l_ITS_DP_ITR_DN
    r_DN = list(map(lambda x, y: x - y, r_ITS_DP_ITR_DN, r_ITS_DP_ITR))
    r_DN_ITR = list(map(lambda x, y: x - y, r_ITS_DP_ITR_DN, r_ITS_DP))
    r_DN_ITR_DP = list(map(lambda x, y: x - y, r_ITS_DP_ITR_DN, r_ITS))
    r_DN_ITR_DP_ITS = r_ITS_DP_ITR_DN
    l_ITS = l_DN_ITR_DP_ITS
    l_ITS_DP = l_DN_ITR_DP
    l_ITS_DP_ITR = l_DN_ITR
    l_ITS_DP_ITR_DN = l_DN
    r_ITS = r_DN_ITR_DP_ITS
    r_ITS_DP = r_DN_ITR_DP
    r_ITS_DP_ITR = r_DN_ITR
    r_ITS_DP_ITR_DN = r_DN
    zorder = [1, 2, 3, 4]
    return l_ITS, l_ITS_DP, l_ITS_DP_ITR, l_ITS_DP_ITR_DN, r_ITS, r_ITS_DP, r_ITS_DP_ITR, r_ITS_DP_ITR_DN, zorder
