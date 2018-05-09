import os


def is_valid_file(parser, arg):
    extension = arg.split('.')[-1]
    if not os.path.isfile(arg):
        parser.error('The file at %s does not exist' % arg)
    elif extension not in ['json', 'txt', 'xlsx']:
        parser.error('Unknown extension %s for file %s please give json, txt or xlsx file' % (extension, arg))
    else:
        return arg
