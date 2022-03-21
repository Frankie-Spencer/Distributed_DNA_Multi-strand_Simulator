import os


def convert_link_address(d):
    chars = ['\\', '/', '_', ':', '-', '.']
    dir = ''

    for c in d:
        if c.isalpha() or c.isdigit() or c in chars:
            dir += c
        else:
            dir += '^' + c
    dir_win = dir.replace('/', '\\')

    return dir_win


def read_file(file_path):
    with open(file_path, 'r') as f:
        all_lines = f.readlines()

        lines = []
        for l in all_lines:
            lines.append(str(l.strip()))

        f.close()

        return lines


def write_file(write_path, str_list):
    with open(write_path, 'w') as f:
        f.write("%s" % '\n')
        for str_line in str_list:
            f.write("%s" % str_line + '\n')

    f.close()


def delete_temp_files(root):
    for path, subdirs, files in os.walk(root):
        for name in files:
            c_dir = os.path.join(path, name)
            if c_dir.endswith('.xml'):
                os.remove(c_dir)
