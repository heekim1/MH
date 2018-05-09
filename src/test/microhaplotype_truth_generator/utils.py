import os


def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        parser.error('The file at %s does not exist' % arg)
    else:
        return arg


def is_valid_min_cov_value(parser, arg):
    try:
        if not 0.0 <= float(arg) <= 1.0:
            parser.error('Please enter a value between 0-1 for mincov. Given value %s' % arg)
        else:
            return float(arg)
    except ValueError:
        parser.error('Please enter a float value for mincov. Given value %s' % arg)
