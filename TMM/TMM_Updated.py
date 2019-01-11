import sys
import numpy as np

from GetReference import GetReference

NUM_OF_REPEATS = 4


def write_to_file_two_conditions_one_repeat(reference, new_sample, genes, reference_index, sample_index, input_dir):
    out_file = open(input_dir + "\\test\\" + "new.libraries.{}.and.{}.results".format(reference_index, sample_index),
                    "w")
    out_file.write("gene_id\texpected_reference_{}\texpected_count_library_{}\n".format(reference_index,
                                                                                        sample_index))
    print("Writing to file....")
    count = 0
    for i in range(reference.shape[0]):
        out_file.write(genes[i] + "\t")
        out_file.write(str(float(reference[i])) + "\t")
        out_file.write(str(float(new_sample[i])) + "\n")
    print("count:" + str(count))


def write_to_file_four_repeats(input_dir, normalized_samples, genes):
    for first_condition in range(7):
        for second_condition in range(first_condition + 1, 7):
            out_file = open(
                input_dir + "\\test\\" + "compare.condition.{}.and.condition{}.results".format(first_condition,
                                                                                               second_condition), "w")
            out_file.write(
                "gene_id\texpected_count_with_all_repeats_condition_{}\texpected_count_with_all_repeats_condition_{}\n".format(
                    first_condition,
                    second_condition))

            first_cond_repeats = np.zeros((NUM_OF_REPEATS, 0))
            for library in range(1, 29):
                if library % NUM_OF_REPEATS == first_condition:
                    first_cond_repeats = np.append(first_cond_repeats, library)

            second_cond_repeats = np.zeros((NUM_OF_REPEATS, 0))
            for library in range(1, 29):
                if library % NUM_OF_REPEATS == first_condition:
                    second_cond_repeats = np.append(second_cond_repeats, library)

            for gene_count in range(genes.shape[0]):
                out_file.write(genes[gene_count] + "\t")
                for repeat_index in range(NUM_OF_REPEATS):
                    out_file.write(str(float(normalized_samples[first_cond_repeats[repeat_index]][gene_count] + "\t")))
                    out_file.write(str(float(normalized_samples[second_cond_repeats[repeat_index]][gene_count] + "\t")))


def calculate_weights_and_factor(log_reference_sample, geo_mean_of_logs):
    sum = np.sum(geo_mean_of_logs)
    weights = np.zeros((geo_mean_of_logs.shape[0],))
    if sum != 0:
        weights = geo_mean_of_logs / np.sum(geo_mean_of_logs)
    scale_factor = np.sum(np.dot(log_reference_sample, weights))
    scale_factor = np.power(2, scale_factor)
    return scale_factor


def find_indeces(x, bottom, top):
    q1 = np.quantile(x, bottom)
    q2 = np.quantile(x, top)

    c1 = np.where(q1 < x)
    c2 = np.where(x < q2)
    c3 = np.where(x != float("-inf"))
    c4 = np.intersect1d(c1, c2)
    d = np.intersect1d(c4, c3)
    return d


def get_vectors_at_indexes(geo_mean_of_logs, log_reference_sample):
    first_indeces = find_indeces(log_reference_sample, 0.3, 0.7)
    second_indeces = find_indeces(geo_mean_of_logs, 0.05, 0.95)
    intersec_indeces = np.intersect1d(first_indeces, second_indeces)
    log_reference_sample = np.take(log_reference_sample, intersec_indeces)
    geo_mean_of_logs = np.take(geo_mean_of_logs, intersec_indeces)

    return geo_mean_of_logs, log_reference_sample


def delete_zero_cells(reference, sample):
    indexes = np.array([]).astype(int)
    for i in range(sample.shape[0]):
        if sample[i] == 0 and reference[i] == 0:
            indexes = np.append(indexes, np.array(i))

    sample = np.delete(sample, indexes)
    reference = np.delete(reference, indexes)

    return reference, sample


def calculate_new_sample(new_library_size, sample):
    sum = np.sum(sample)
    new_sample = np.zeros((sample.shape[0]))
    if sum != 0:
        new_sample = (sample / np.sum(sample)) * new_library_size
    return new_sample


def get_data(reference_library, sample):
    reference_index = reference_library.file_index
    reference = np.loadtxt(reference_library.file_name, dtype=str, delimiter="\t")[1:, 4].astype(np.float)
    sample_index = sample.file_index
    genes = np.loadtxt(sample.file_name, dtype=str, delimiter="\t")[1:, 0]
    sample = np.loadtxt(sample.file_name, dtype=str, delimiter="\t")[1:, 4].astype(np.float)
    return genes, reference, reference_index, sample, sample_index


def run(input_dir):
    get_reference = GetReference(input_dir)
    reference_library, samples = get_reference.count_read_scaled_counts()
    normalized_samples = np.zeros(28, samples.shape[0])

    for i, sample in enumerate(samples):
        if i <= 2:
            genes, reference, reference_index, sample, sample_index = get_data(reference_library, sample)
            reference, sample = delete_zero_cells(reference, sample)

            sample_sum = np.sum(sample)
            log_reference_sample = np.zeros((reference.shape[0]))
            if sample_sum != 0:
                log_reference_sample = np.log(reference / (np.sum(reference)) / (sample / np.sum(sample)))

            geo_mean_of_logs = (np.log(reference) + np.log(sample)) / 2
            geo_mean_of_logs, log_reference_sample = get_vectors_at_indexes(geo_mean_of_logs, log_reference_sample)

            scale_factor = calculate_weights_and_factor(log_reference_sample, geo_mean_of_logs)
            new_library_size = np.sum(sample) * scale_factor

            new_sample = calculate_new_sample(new_library_size, sample)
            normalized_samples = np.append(normalized_samples, np.array(new_sample), axis=i)

            write_to_file_two_conditions_one_repeat(reference, new_sample, genes, reference_index, sample_index,
                                                    input_dir)

            # import matplotlib.pyplot as plt
            # plt.plot(new_sample, reference, 'bs')
            # plt.title("Library{}_VS_Reference".format(sample_index))
            # plt.xlim([0, 10000])
            # plt.ylim([0, 10000])
            # plt.ylabel("Reference_{}".format(reference_index))
            # plt.xlabel("Library{}".format(sample_index))
            # plt.show()
            # plt.savefig(input_dir + "\\test\\" + "Library{}_VS_Reference".format(sample_index))

    # write_to_file_two_conditions_one_repeat(input_dir, normalized_samples, genes)
    print(1)


if __name__ == '__main__':
    run(sys.argv[1])
