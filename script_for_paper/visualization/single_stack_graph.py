# importing package
import os
from enum import Enum
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


# pattern_enum
class Bar(Enum):
    F2C = 1
    F2D = 2
    SF1B = 3
    SF5B = 4
    F3B = 5
    F5C = 6
    F3D = 7
    F4D = 8
    F4E = 9
    SF2F = 10
    SF2G = 11
    SF5A = 12


class UniversalColor(Enum):
    RED = '#FF4B00'
    BLUE = '#005AFF'
    GREEN = '#03AF7A'
    LIGHT_BLUE = '#4DC4FF'
    ORANGE = '#F6AA00'
    YELLOW = '#FFF100'
    BLACK = '#000000'
    PURPLE = '#990099'
    GRAY = '#84919E'
    INDIGO = '#332288'
    CYAN = '#88CCEE'
    PALE_GRAY = '#DDDDDD'


def get_color_list(pattern):
    uv = UniversalColor
    if pattern == Bar.F2C or pattern == Bar.SF2F:
        return ["#F5BDD3", uv.RED.value, "#F5BDD3"]
    if pattern == Bar.F2D or pattern == Bar.SF2G:
        # return [uv.INDIGO.value, uv.CYAN.value, uv.CYAN.value, uv.CYAN.value, uv.PALE_GRAY.value, uv.PALE_GRAY.value, uv.PALE_GRAY.value]
        return [uv.INDIGO.value, uv.LIGHT_BLUE.value, uv.LIGHT_BLUE.value, uv.LIGHT_BLUE.value, '#F2FEFF', '#F2FEFF', '#F2FEFF']
    if pattern == Bar.SF1B:
        return [uv.BLACK.value, uv.BLUE.value, uv.BLUE.value, uv.LIGHT_BLUE.value, uv.LIGHT_BLUE.value, uv.GREEN.value, uv.GREEN.value, uv.ORANGE.value, uv.RED.value]
    if pattern == Bar.SF5B:
        return [uv.GREEN.value, uv.ORANGE.value, uv.ORANGE.value, uv.ORANGE.value]
    if pattern == Bar.F3B:
        return [uv.BLACK.value, uv.RED.value, uv.GREEN.value, uv.BLUE.value]
    if pattern == Bar.F5C or pattern == Bar.SF5A:
        return [uv.RED.value, uv.GREEN.value, uv.BLUE.value]
    if pattern == Bar.F3D:
        return [uv.RED.value, uv.BLUE.value]
    if pattern == Bar.F4D:
        return ['#F2FEFF', uv.LIGHT_BLUE.value]
    if pattern == Bar.F4E:
        return [uv.RED.value, uv.GREEN.value, uv.BLUE.value, uv.BLACK.value]


    return ['darkgreen', 'tomato', 'blue', 'darkgreen', 'tomato', 'blue']


def get_hatch(pattern):
    if pattern == Bar.F2C or pattern == Bar.SF2F:
        return ['/', '', '\\']
    if pattern == Bar.F2D or pattern == Bar.SF2G:
        return ['', 'x', '/', '\\', '/', '\\', '']
    if pattern == Bar.SF1B:
        return ['', '', '\\', '', '\\', '', '\\', '', '']
    if pattern == Bar.SF5B:
        return ['', '/', '', '\\']
    return ['']


def gen_single_stacked_bar_graph(x_labels, data_labels, data, pattern, dir_name):
    # load dataset
    df = pd.DataFrame(index=x_labels)
    for i, label in enumerate(data_labels):
        df[label] = data[i]

    # view dataset
    print(df)

    # remove margin around the graph
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0

    # plot a Stacked Bar Chart using matplotlib
    ax = df.plot(
        kind='bar',
        stacked=True,
        # figsize=(4, 15),
        figsize=(1.25, 4.5),
        width=1,
        color=get_color_list(pattern),
        edgecolor=UniversalColor.BLACK.value, linewidth=1.5,
    )

    # hatching
    hatches = get_hatch(pattern)
    for i, patch in enumerate(ax.patches):
        patch.set_hatch(hatches[i % len(hatches)])

    # ax.yaxis.set_major_formatter(mtick.PercentFormatter())

    # plt.legend(data_labels, fontsize='xx-large', loc='upper left', bbox_to_anchor=(1, 1))
    plt.legend('', frameon=False)
    plt.tick_params(labelsize=12)
    plt.xticks(rotation=0)
    plt.xlim(-0.5, 0.5)
    plt.tight_layout()

    output_dir = "fig_out"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(output_dir + '/' + dir_name):
        os.makedirs(output_dir + '/' + dir_name)
    plt.savefig(output_dir + '/' + dir_name + '/' + x_labels[0] + '.tiff')

    plt.show()


def set_config():
    _d_name = input("directory name: ")
    _data_labels = input("data labels: ").split(',')
    _data_labels = [x for x in _data_labels if x != '']
    print(_data_labels[1:])
    return _data_labels[1:], _d_name


def set_data():
    _data_row = input("data row(type z to exit): ").split(',')
    if _data_row[0] == 'z':
        exit(0)
    _data_row = [x for x in _data_row if x != '']
    _data_row_int = [int(x) for x in reversed(_data_row[1:])]
    print(_data_row)
    return [_data_row[0]], _data_row_int


def to_percentage(data):
    total = sum(data)
    return list(map(lambda x: x / total * 100, data))


def main():
    # data_labels = ['A', 'B', 'C']
    # data_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
    # data_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    # data_labels = ['A', 'B', 'C', 'D']
    # data = [10, 20, 70]
    # x_label = ['one']

    data_labels, dir_name = set_config()
    pattern = Bar.F3D

    while True:
        x_label, data = set_data()
        data = to_percentage(data)
        gen_single_stacked_bar_graph(x_labels=x_label, data_labels=data_labels, data=data, pattern=pattern,
                                     dir_name=dir_name)


if __name__ == "__main__":
    main()
