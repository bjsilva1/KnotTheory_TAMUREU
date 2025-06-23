from collections import namedtuple

LinkInfo = namedtuple('LinkInfo', ['num_components', 'alexander_nullity', 'jones_nullity', 'seifert_nullity', 'signature', 'linking_matrix_zero'])

ALL_LINKS = dict()

with open("KnotTheory_TAMUREU/link_info.txt") as infile:
    for line in infile.readlines():
        link_info = eval(line)
        ALL_LINKS[link_info[0]] = LinkInfo(link_info[1], link_info[2], link_info[3], link_info[4], link_info[5], link_info[6])