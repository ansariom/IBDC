#!/local/cluster/bin/python
import numpy as np
import pickle
import io
from sklearn import preprocessing
from sklearn import linear_model
import pandas as pd
import re
import datetime
import os
import argparse
import sys

# to help with raw-formatting help text...
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)



parser = argparse.ArgumentParser(description = "Runs 'in-silico knockdowns', taking a specification for which features to adjust and which TSSes to report results for. See http://megraw.cgrb.oregonstate.edu/wiki/index.php/IBDC_pipeline#File_format_suggestion for format details.", formatter_class = SmartFormatter)
parser.add_argument("--scaled_features", dest = "scaled_features_filename", required = True, help = "Scaled features file, tab separated.")
parser.add_argument("--diff_info", dest = "diff_info_filename", required = True, help = "Differential expression info file, comma-separated (Mitra's fault :-P).")
parser.add_argument("--feature_info", dest = "feature_info_filename", required = True, help = "Feature info file, tab-separated.")
parser.add_argument("--model", dest = "model_filename", required = True, help = "Sklearn model file.")
parser.add_argument("--knock_spec", dest = "knock_spec_filename", required = True, help = "R|A file (tab-separated) indicating which features to \n \
knockdown/adjust, and for which TSSes, and how. Example: \n \
    experiment_id   feature                                set_to \n \
    ex1             ^TATA_FWD_1$                           0 \n \
    ex2             ^TATA_FWD_1$                           0 \n \
    ex2             ^TATA_FWD_7$                           0 \n \
    ex3             ^TATA_(FWD|REV)_.$                     0.5 \n \
    ex4             ^TATA_(FWD|REV)_._OC_P_LEAF$           1 \n \
Here, ex1 just sets the TATA_FWD_1 feature to 0 for all TSSes.  \n \
 \n \
ex2 sets all the tss rows to 0 for the features TATA_FWD_1 and TATA_FWD_7.  \n \
 \n \
ex3 sets all features matching the given regex (in this case the TFBS  \n \
features for TATA, both forward and reverse) to 0.5. \n \
 \n \
ex4 sets all TSS rows to 1 for all the features matching the regex \n \
(here the OC features in leaf for TATA).") 
action = parser.add_mutually_exclusive_group(required = True)
action.add_argument("--tss_list", dest = "tss_list_filename", help = "Single column file (with header, ignored) indicating which TSS probability changes should be repored in the output.")
action.add_argument("--tss_percentile", dest = "tss_percentile", help = "Report the top percentile of TSSes, by absolute value of change. e.g. 0.05.")

args = parser.parse_args()

## Load TSS list
tss_list = None
tss_percentile = None
if args.tss_list_filename:
    tss_list_table = pd.read_csv(args.tss_list_filename, sep = "\t")
    tss_list = list(tss_list_table.loc[:, 'gene_id'])
else:
    tss_percentile = float(args.tss_percentile)
    assert tss_percentile > 0 and tss_percentile < 1, "--tss_percentile should be between 0 and 1 (e.g. 0.05)"



### read the data!
sys.stderr.write("Loading data... " + str(datetime.datetime.time(datetime.datetime.now())) + "\n")
all_features_scaled = pd.read_csv(args.scaled_features_filename, sep = "\t")
all_tss_diff_info = pd.read_csv(args.diff_info_filename, sep = ",")
feature_info = pd.read_csv(args.feature_info_filename, sep = "\t")
# load the model
logreg = pickle.load(io.open(args.model_filename, "r"))
# load spec files
knock_spec = pd.read_csv(args.knock_spec_filename, sep = "\t")



# get a numpy array of all the features, extract row and column names as lists
rownames = list(all_features_scaled['tss_name'])
colnames = all_features_scaled.columns.get_values()
# drop "tss_name" from the colnames
colnames = colnames[1:]

all_features_scaled = all_features_scaled.drop('tss_name', axis = 1)
features_array = all_features_scaled.as_matrix()




## Parse knock_spec to produce a set of experiments...
## which will be a dictionary mapping experiment ids to pandas data frames with the exp info
expids = set(list(knock_spec['experiment_id']))
experiments = dict()
for expid in expids:
    rows = knock_spec.loc[knock_spec['experiment_id'] == expid ,]
    experiments[expid] = rows


sys.stderr.write("Data Loaded! " + str(datetime.datetime.time(datetime.datetime.now())) + "\n")


# Pandas! Grrr. https://gist.github.com/hellpanderrr/599bce82ecc6934aa9e1 with modifications for sorting by column name
def sort_df(df, column_name, key):
    '''Takes dataframe, column index and custom function for sorting, 
    returns dataframe sorted by this column using this function'''
    
    col = df.loc[:,column_name]
    temp = np.array(col.values.tolist())
    order = sorted(range(len(temp)), key=lambda j: key(temp[j]))
    return df.iloc[order]



## scaled_features_np: numpy array of features
## rownames: rownames of the numpy array
## colnames: colnames of the numpy array
## tss_names: the TSSes we want to see the changes in
## knock_spec: pandas data frame with cols experiment_id, feature, tss, and set_to
## e.g.
#   experiment_id        feature tss set_to
# 1           ex2  TATAbox_FWD_1  .+      0
# 2           ex2  TATAbox_FWD_7  .+      0
## (the experiment_id column is ignored...)
## model: a sklearn model capable of .predict_proba()
## exp_id: an identifier for this knockdown experiment
## RETURNS: a pandas dataframe
def knockdown(scaled_features_np, rownames, colnames, tss_names, knock_spec, model, exp_id, tss_percentile):
    assert tss_names != None or tss_percentile != None, "Error: knockdown() called but both tss_names (list of string) and tss_percentile (float) are None. Exactly one must be set."
    # check the shape of this stuff: throw an error if the rownames and colnames arent the right size:
    if scaled_features_np.shape[0] != len(rownames) or scaled_features_np.shape[1] != len(colnames):
        sys.stderr.write("knockdown function called, but the number of rows and columns in the feature matrix is not equal to the length of rownames or colnames. Something is wrong. Checkit.\n")
        sys.stderr.write("Matrix size: " + str(scaled_features_np.shape) + "\n")
        sys.stderr.write("rowname and colnames lengths: " + str(len(rownames)) + "," + str(len(colnames)) + "\n")
        quit()

    # We gotta make a deep copy of the array before we go and modify it...	
    scaled_features_np = np.array(scaled_features_np)

    # run pre-knocked model (prob-1 values are at index 1 of the element lists)
    predicted_pre = model.predict_proba(scaled_features_np)
    one_class_predictions_pre = [x[1] for x in predicted_pre]

    # for each row in the knockdown specification:
    for index, row in knock_spec.iterrows():
        # extract the patterns and set_to
        col_regex = row['feature']
        set_to = float(row['set_to'])
        
        # consider every row and column of the matrix; if both names match,
        # set the entry to set_to
        #for row_index in xrange(0, scaled_features_np.shape[0]):
            #sys.stderr.write("Churning over columns for row " + str(row_index) + "\n")
            #if re.search(row_regex, rownames[row_index]):
        ## The above was slooooowww; no need to consider rows individually, since the model 
        ## is run on all rows (indpendently) anyway
        for col_index in xrange(0, scaled_features_np.shape[1]):
            if re.search(col_regex, colnames[col_index]):
                sys.stderr.write("  Setting col " + colnames[col_index] + " to " + str(set_to) + "\n")
                scaled_features_np[:, col_index] = set_to
                        
    
    # run post-knocked model
    predicted_post = model.predict_proba(scaled_features_np)
    one_class_predictions_post = [x[1] for x in predicted_post]

    # compute the differences as a list
    diff = [pre - post for pre, post in zip(one_class_predictions_pre, one_class_predictions_post)]

    results = pd.DataFrame({'experiment_id': exp_id,
                        'tss_name':rownames, 
                        'pre_knocked_prob1': one_class_predictions_pre,
                        'post_knocked_prob1': one_class_predictions_post,
                        'diff': diff})

    if tss_names:
        # make a dictionary of tss_names we want to report for (for speedy lookup)
        select_dict = dict()
        for tss in tss_names:
            select_dict[tss] = 1    
        
        # grab all the rows that are keys in in the select_dict
        selector = [bool(x in select_dict) for x in results['tss_name']]
        results_selected = results.loc[selector, :]
        return results_selected
        
    results = sort_df(results, 'diff', abs)

    rowcount = results.shape[0]
    num_to_keep = int(rowcount * tss_percentile)
    return results.iloc[(rowcount - num_to_keep):(rowcount),:]


i = 0
## run da experiments
for expid in experiments.keys():
    exp_df = experiments[expid]
    sys.stderr.write("Running experiment " + expid + " " + str(datetime.datetime.time(datetime.datetime.now())) + "\n")
    res_df = knockdown(features_array, rownames, colnames, tss_list, exp_df, logreg, expid, tss_percentile)

    if i == 0:
        res_df.to_csv(sys.stdout, header = True, sep = "\t", index = False)  # index = False means no row names
    else:
        res_df.to_csv(sys.stdout, header = False, sep = "\t", index = False)  # index = False means no row names

    i = i + 1



sys.stderr.write("All done! " + str(datetime.datetime.time(datetime.datetime.now())) + "\n")
