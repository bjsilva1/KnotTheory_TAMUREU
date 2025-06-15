from collections import namedtuple

LinkInfo = namedtuple('LinkInfo', ['name', 'alexander_nullity', 'jones_nullity', 'seifert_nullity', 'signature', 'linking_matrix_zero'])

ALL_LINKS = []

with open("KnotTheory_TAMUREU/links_info.txt") as infile:
    for line in infile.readlines():
        link_info = eval(line)
        ALL_LINKS.append(LinkInfo(link_info[0], link_info[1], link_info[2], link_info[3], link_info[4], link_info[5]))