import math
import os
import numpy as np


class GetReference:
    def __init__(self, dir):
        self.dir = dir
        self.files = np.array([]).astype(FileWrap)
        self.quintiles = np.array([])
        self.get_library_files(dir)
        self.count_read_scaled_counts()

    def get_library_files(self, dir_name):
        for file_name in os.listdir(dir_name):
            full_path = os.path.join(dir_name, file_name)
            if not os.path.isfile(full_path):
                continue
            if file_name.endswith(".genes.results"):
                file = open(full_path, "r")
                self.files = np.append(self.files, np.array([FileWrap(file, file_name, dir_name)]))
        self.files = np.sort(self.files)

    def get_quintile(self, genes_count):
        sorted_genes_count = np.sort(genes_count)
        return sorted_genes_count[math.floor((0.75) * sorted_genes_count.shape[0]) - 1]

    def count_read_scaled_counts(self):
        for i, file in enumerate(self.files):
            if file.file_name.endswith(".genes.results"):
                genes_count = np.loadtxt(file.file_name, dtype=str, delimiter="\t")[1:, 4].astype(
                    np.float)
                col_sum = np.sum(genes_count)
                if col_sum > 0:
                    genes_count /= col_sum
                self.quintiles = np.append(self.quintiles, self.get_quintile(genes_count))

        quintiles_average = np.average(self.quintiles)
        diff = np.abs(quintiles_average - self.quintiles)
        min_diff = np.argmin(diff)
        return self.files[min_diff], self.files


class FileWrap:

    def __init__(self, file, file_name, dir):
        self.file = file
        self.file_name = dir + "\\" + file_name
        self.file_index = int(file_name.split(".")[1])

    def __gt__(self, other):
        return self.file_index > other.file_index


if __name__ == '__main__':
    # getReference = GetReference("C:\\Users\\Danielle\\Documents\\SharedUbuntu\\output")
    getReference = GetReference("C:\\Users\\Danielle\\Documents\\SharedUbuntu\\output\\test")
