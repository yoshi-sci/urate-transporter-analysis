# importing package
from csv_reader import csv_reader

from plt_bar import plt_bar


def main():
    # define file_path for the csv data
    file_path = 'data/Test-data.csv'

    # read csv file
    database = csv_reader(file_path)

    # generate figures
    reverse = True
    style = 'quadruplets'  # style = 'twins'
    for i in [1, 2, 3]:
        option = i
        plt_bar(database[0], option, reverse, style)
        plt_bar(database[1], option, reverse, style)
        plt_bar(database[2], option, reverse, style)
        plt_bar(database[3], option, reverse, style)


if __name__ == "__main__":
    main()
