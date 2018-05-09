import os


def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        parser.error('The file at %s does not exist' % arg)
    else:
        return arg


def is_valid_dir(parser, arg):
    if not os.path.isdir(arg):
        parser.error('The directory at %s does not exist' % arg)
    else:
        return arg
