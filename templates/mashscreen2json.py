#!/usr/bin/env python3

from statistics import median
import json

if __file__.endswith(".command.sh"):
    MASH_TXT = '$mashtxt'

def main(mash_output):
    '''converts top results to json

    Parameters
    ----------
    mash_output: str
        this is a string that stores the path to this file

    '''

    read_mash_output = open(mash_output)

    dic = {}
    median_list = []

    for line in read_mash_output:
        #print(line)
        tab_split = line.split("\t")
        identity = tab_split[0]
        #shared_hashes = tab_split[1]
        median_multiplicity = tab_split[2]
        #p_value = tab_split[3]
        query_id = tab_split[4]
        #query-comment should not exist here and it is irrelevant

        # here identity is what in fact interests to report to json but
        # median_multiplicity also is important since it gives an rough
        # estimation of the coverage depth for each plasmid.
        # Plasmids should have higher coverage depth due to their increased
        # copy number in relation to the chromosome.
        dic[query_id] = [identity, median_multiplicity]
        median_list.append(float(median_multiplicity))

    # median cutoff is twice the median of all median_multiplicity values
    # reported by mash screen. In the case of plasmids, since the database
    # has 9k entries and reads shouldn't have that many sequences it seems ok...
    median_cutoff = median(median_list)

    output_json = open(" ".join(mash_output.split(".")[:-1]) + ".json", "w")

    filtered_dic = {}
    for k,v in dic.items():
        # estimated copy number
        copy_number = int(float(v[1])/median_cutoff)
        # assure that plasmid as at least twice the median coverage depth
        if float(v[1]) > median_cutoff:
            filtered_dic["_".join(k.split("_")[0:3])] = [v[0], str(copy_number)]
    print("Exported dictionary has {} entries".format(len(filtered_dic)))
    output_json.write(json.dumps(filtered_dic))
    output_json.close()

if __name__ == "__main__":
    # a variable from nextflow process
    main(MASH_TXT)