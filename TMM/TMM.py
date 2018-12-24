import math
import sys
import pandas as pd
import numpy as np


def write_to_file(N_k_1, N_k_2, Y_k_1, Y_k_2, first_library_index, genes, second_library_index):
    out_file = open("libraries.{}.and.{}.results".format(first_library_index, second_library_index), "w")
    out_file.write("gene_id\texpected_count_library_{}\texpected_count_library_{}\n".format(first_library_index,
                                                                                            second_library_index))
    print("Writing to file....")
    for i in range(1, Y_k_1.shape[0]):
        print(i)
        out_file.write(genes[i] + "\t")
        out_file.write(str(float(Y_k_1[4][i]) / N_k_1) + "\t")
        out_file.write(str(float(Y_k_2[4][i]) / N_k_2) + "\n")


def get_M_g_v(N_k_1, N_k_2, Y_k_1, Y_k_2):
    M_g_v = np.zeros((Y_k_1.shape[0], 1))

    for i in range(1, Y_k_1.shape[0]):
        Y_g_k_1 = float(Y_k_1[4][i])
        Y_g_k_2 = float(Y_k_2[4][i])

        if Y_g_k_1 == 0 or Y_g_k_2 == 0:
            M_g = 0

        else:
            M_g = np.log(Y_g_k_1 / N_k_1) / np.log(Y_g_k_2 / N_k_2)
            W_g = (N_k_1 - Y_g_k_1) / (N_k_1 * Y_g_k_1) + (N_k_2 - Y_g_k_2) / (N_k_2 * Y_g_k_2)
            M_g = M_g * W_g

        M_g_v[i] = M_g

    return M_g_v



def get_files(input_file_1, input_file_2):
    file_1 = pd.read_csv(input_file_1, header=None, sep="\t")
    file_2 = pd.read_csv(input_file_2, header=None, sep="\t")

    df_file_1 = pd.DataFrame(file_1)
    df_file_2 = pd.DataFrame(file_2)

    return df_file_1, df_file_2


def trim(M_g_v):
    M_g_v_sort = np.sort(M_g_v, axis=None)
    M_g_v_sort_without_zeros = M_g_v_sort[np.nonzero(M_g_v_sort)]
    n = math.floor(0.7 * M_g_v_sort_without_zeros.shape[0])
    return M_g_v_sort_without_zeros[:n]


def sum_vector(Y_k):
    sum = 0
    for i in range(1, Y_k.shape[0] - 1):
        sum += float(Y_k[4][i])
    return sum


def get_sum_M_g_v_f(M_g_v_f):
    sum = 0
    for i in range(len(M_g_v_f)):
        sum += M_g_v_f[i]
    return sum


def get_log_TMM_k(Y_k_1, Y_k_2, N_k_1, N_k_2):
    log_TMM_k = 0

    for i in range(1, Y_k_1.shape[0]):
        Y_g_k_1 = float(Y_k_1[4][i])
        Y_g_k_2 = float(Y_k_2[4][i])
        if Y_g_k_1 != 0 and Y_g_k_2 != 0:
            log_TMM_k += (N_k_1 - Y_g_k_1) / (N_k_1 * Y_g_k_1) + (N_k_2 - Y_g_k_2) / (N_k_2 * Y_g_k_2)

    return log_TMM_k


def run(input_file_1, input_file_2):
    first_library_index = input_file_1.split(".")[1]
    second_library_index = input_file_2.split(".")[1]

    print("Collecting files....")
    df_file_1, df_file_2 = get_files(input_file_1, input_file_2)

    Y_k_1 = df_file_1.iloc[:, 4:5]
    N_k_1 = sum_vector(Y_k_1)

    Y_k_2 = df_file_2.iloc[:, 4:5]
    N_k_2 = sum_vector(Y_k_2)

    genes_ = df_file_2.iloc[:, 0:1]
    genes = genes_.values.flatten()

    M_g_vector = get_M_g_v(N_k_1, N_k_2, Y_k_1, Y_k_2)
    trimmed_M_g = trim(M_g_vector)

    sum_M_g_v_f = get_sum_M_g_v_f(trimmed_M_g)
    log_TMM_k = get_log_TMM_k(Y_k_1, Y_k_2, N_k_1, N_k_2)

    f_k = np.power(2, sum_M_g_v_f / log_TMM_k)
    f_k_s_r = np.sqrt(f_k)

    N_k_1 = N_k_1 / f_k_s_r
    N_k_2 = N_k_2 * f_k_s_r

    write_to_file(N_k_1, N_k_2, Y_k_1, Y_k_2, first_library_index, genes, second_library_index)


if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])

