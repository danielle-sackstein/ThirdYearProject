import collections
import os
import gzip
import shutil


class LibraryFile:
    def __init__(self, file_name, dir_name):
        self.file_name = file_name
        self.dir_name = dir_name

        results = file_name.split('_')
        self.lib_index = int(results[1][1:])
        self.part_index = int(results[2][1:])
        self.direction_index = int(results[3][1:])

    def __str__(self):
        return "{}, {}, {}: {}".format(self.lib_index, self.part_index, self.direction_index, self.file_name)


def get_library_files(dir_name):
    for file_name in os.listdir(dir_name):
        full_path = os.path.join(dir_name, file_name)
        if not os.path.isfile(full_path):
            continue
        if file_name.endswith(".fastq.gz"):
            yield LibraryFile(file_name, dir_name)


def get_libraries(dir_name):
    library_map = collections.defaultdict(list)
    for library_file in get_library_files(dir_name):
        # if library_file.direction_index == 2:
        if library_file.direction_index == 1:
            library_map[library_file.lib_index].append(library_file)
    return library_map


def concatenate(concat_library_file_name, library):
    with open(concat_library_file_name, 'w') as f_out:
        for library_file in library:
            print('\treading {}'.format(library_file))
            file_path = os.path.join(library_file.dir_name, library_file.file_name)
            with gzip.open(file_path, 'rt') as f_in:
                for line in f_in:
                    f_out.write(line)


def create_empty_dir(dir_name):
    if os.path.exists(dir_name):
        shutil.rmtree(dir_name)
    os.mkdir(dir_name)


def run():
    dir_name = 'C:\\Users\\Danielle\\Documents\\SharedUbuntu\\SeqData'

    concat_library_dir_name = '{}\\concat'.format(dir_name)
    # create_empty_dir(concat_library_dir_name)
    # TODO: 28 library and 9
    library_map = get_libraries(dir_name)
    for index in range(1, 29):
        concat_library_file_name = '{}\\library.{}.fasta'.format(concat_library_dir_name, index)
        print('Creating {}...'.format(concat_library_file_name))
        concatenate(concat_library_file_name, library_map[index])
        print('Finished {}'.format(concat_library_file_name))


if __name__ == '__main__':
    run()
    print('done')
