import numpy as np
from numpy import dtype

# 1- compare matrices in several steps:
# a: element-wise comparison and detect those that have less than T distance
# b: Find shared core :
# b1: Consider matrices having less than T1 differnce in teir row numbers (n1-n2) < T1
# b2: 

dist_threshold = 0.11
shape_diff = 4

def load_matrices(infile):
    mtx_dict = {}
    mtx_file = open(infile)
    for line in iter(mtx_file):
        line = line.rstrip()
        if (line.startswith(">")):
            pwm_name = line[2:]
            mtx_arr = []
            continue
        if(line.startswith("=")):
            line_parts = line.split("\t")
            mtx_arr.append(line_parts[1:])
        else:
            arr = np.array(mtx_arr, dtype=float)
            mtx_dict[pwm_name] = arr
    return mtx_dict

def compare_matrices(mtx_dict):
    print("load")
    outdir = "sim_pwms/"
    #outfile = open("similar_pwms.txt", "w")
    for j in xrange(0,len(mtx_dict.keys())):
	name = mtx_dict.keys()[j]
        mtx1 = mtx_dict[mtx_dict.keys()[j]]
	outfile = outdir + name + ".txt"
	f = open(outfile, "w")
	f.write(mtx1)
	f.close()
        

if __name__ == "__main__":
    infile = "all_PWMs.mat"
    mtx_dict = load_matrices(infile)
    compare_matrices(mtx_dict)
