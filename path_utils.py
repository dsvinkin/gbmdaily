import os
import logging as log

#import config
#info = config.read_config('config.yaml')

#import setlog
#setlog.set_log()

def file_is_ok(file_name, min_size_kb):

    if file_name is None:
        return False

    size = 0
    if os.path.isfile(file_name): 
        size = os.path.getsize(file_name)
    else:
        log.info("File {:s} not found.".format(file_name))
        return False
    
    log.info("File {:s} has size: {:d} B".format(file_name, size))

    if (size >= min_size_kb * 1024):
        return True
    else:
        return False


def get_files(path, pattern='', prefix=True, all=False):
    """
    Find file with specific prefix (prefix=True) or suffix (prefix=False) 
    """

    list_files = os.listdir(path)

    if len(list_files) == 0:
        log.info("Directory {:s} is empty!".format(path))
        return None

    if prefix:
        file_folder = list(filter(lambda x: x.startswith(pattern), list_files))
    else:
        file_folder = list(filter(lambda x: x.endswith(pattern), list_files))

    if len(file_folder) == 0:
        log.info("No required file with pattern: {:s} in {:s}.".format(pattern, path))

    #print(file_folder)

    if all:
        return file_folder
    elif not all and len(file_folder) > 0:
        return file_folder[0]
    else:
        return None 