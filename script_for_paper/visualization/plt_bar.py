from plt_quadruplets_bar import plt_quadruplets_bar
from plt_twins_bar import plt_twins_bar


def plt_bar(database, option, reverse, style):
    if style == 'twins':
        plt_twins_bar(database, option, reverse)
    elif style == 'quadruplets':
        plt_quadruplets_bar(database, option, reverse)
