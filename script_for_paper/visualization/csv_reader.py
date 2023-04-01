def csv_reader(file_path):
    import csv

    row_items = []
    with open(file_path) as f:
        reader = csv.reader(f)
        for row in reader:
            row_items.append(row[0:5])

    dataset = [row_items[0:9], row_items[10:19], row_items[20:29], row_items[30:39]]

    return dataset
