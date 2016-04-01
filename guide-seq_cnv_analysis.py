import pandas
import numpy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.stats as stats
#######################################
# USER: Set input paths here

# Guide-seq data
guide_seq_data_file = "data/crispr_guide_data.tsv"

# Cell line CNV data. Supplied as log10 ratios from GEO.
cnv_data_file = "data/GSM952580_full.txt"

# Manufacturer-supplied descriptor of CNV data array
array_index_file = "data/GPL9777-1843.txt"

# END parameter setting
#######################################
guide_seq_data = pandas.read_table(open(guide_seq_data_file, "r"),\
                          dtype={'gene_name': numpy.str})
cnv_data = pandas.read_table(open(cnv_data_file, "r"), converters={\
    'ID_REF': numpy.int},engine="python", skiprows = 4, skipfooter=1)
array_index = pandas.read_table(open(array_index_file ,"r"), dtype={'ID':\
    numpy.int}, skiprows = 12)
    #comment="#)

# Match up CNV data entries (which are supplied as numeric gene IDs)
# with the index describing how those IDs map to gene names
print("Matching CNV data to array key...")
cnv_ratios_and_genenames = pandas.concat([cnv_data, array_index], axis=1)

# A potential error is that missing data may cause erroneous matches between
# the array key and the CNV dataset. This block reports the first non-matching
# gene, or reports an confirmation if necessary.
for item in cnv_ratios_and_genenames.iterrows():
    if item[1]["ID_REF"] == item[1]["ID"]:
        None
    else:
        print(item[1])
        print("\n")
        print(repr(item[1]["ID_REF"]))
        print(repr(item[1]["ID"]))
        break
else:
    print("\tArray keys matched and verified.")

# Next take the guide data and extract that subset of the CNV data that overlaps
# with the guide data.

# At this point we need to take the CNV data that intersects with regions
# covered by guides. The ideal here would to be map each guide back to genomic
# coordinates and correlate those with the coordinates supplied in the CGH array
# data. However, for the sake of completing this in a reasonable amount of time,
# I'm going to trust that the `gene_name` in the guide data can be intersected
# with the `GENE_SYMBOL` in the CNV data.

# So, next we take the CNV data entries that have `GENE_SYMBOL` set (i.e. not an
# `NaN`).

cnv_ratios_and_genenames = cnv_ratios_and_genenames.dropna(\
    subset=["GENE_SYMBOL"])

# The next step relies of being able to intersect the two datasets; ideally all
# the gene names in the sgRNA data are represented in the CGH data. We're also
# going to proceed differently depending on whether genes in the CGH dataset are
# represented in more than one array spot.

# First task:

# The CNV data contains multiple array spots for different parts of each gene.
# We really want to reduce these data to a single median value per gene, as
# genes are the atomic unit in our sgRNA data. ("VALUE" is the column containing
# the log10 ratios of control to A375 spot intensity). log10 Cy5/Cy3, i.e.
# A375/control. Greater than 0 is enriched in A375; less than 1 is depleted in
# A375.

# Collapse the values for each gene into a single median gene value
cnv_grouped_by_gene = cnv_ratios_and_genenames[["VALUE", "GENE_SYMBOL"]]\
    .groupby('GENE_SYMBOL').agg([numpy.median, numpy.std, len])

# Second task:
# Figure out how much guides targeting each gene in the guide-seq data are
# depleted or enriched at each timepoint. This is a little more complicated as
# we first need to normalize the cell data against the starting plasmid
# concentration. Looking at the data, the minimum value in the guide-seq data is
# 1, indicating that the rest of the data is likely a multiple of that base
# value. Simply dividing the sequenced cell-guide abundance by the plasmid input
# data and log transforming for each gene should be sufficient.
guide_abundance = guide_seq_data.iloc[:,4:].div(guide_seq_data\
                                            ["norm_count_plasmid"],\
                                            axis='index').applymap(numpy.log)

# Merge the first column of the original data (i.e. the gene labels) with the
# newly normalized data
guide_abundance = pandas.merge(pandas.DataFrame(guide_seq_data.iloc[:,0]), \
    guide_abundance, left_index=True, right_index=True)

# Now summarize by gene name as was done with the CNV data
guide_abundance = guide_abundance.groupby('gene_name').agg([numpy.median, \
    numpy.std, len])

# Next concern is how well the datasets map onto each other. Are there many
# genes targeted by guides that for some reason aren't in the CNV data, or
# vice-versa?
print(str(len(cnv_grouped_by_gene)) + " genes found in the CNV dataset " + \
    "supplied.")
print(str(len(guide_abundance)) + " genes found in the guide-seq dataset " + \
    "supplied.")

# So, there are 17808 genes in the CNV dataset, but only 17419 in the sgRNA
# dataset. That's not unexpected; some sgRNAs are likely to be totally lost in
# library preparation, whereas all intended spots are imaged on the CNV array.

# Merge those rows containing genes in both CNV and guide-seq datasets:
cnv_and_guides = pandas.merge(cnv_grouped_by_gene, guide_abundance, \
    left_index=True, right_index=True, how="inner")
print(str(len(cnv_and_guides)) + " genes occur in both the CNV and " + \
    "guide-seq datasets. Analysis will continue with these genes.")

# Show how many unique genes are seen in both datasets:
print("(A total of " + str(len(pandas.merge(cnv_grouped_by_gene, \
    guide_abundance, left_index=True, right_index=True, how="outer"))) + \
    " unique genes were encountered across both datasets.)")

#Rename some columns to help keep track.
cnv_and_guides =cnv_and_guides.rename(columns={'VALUE': 'CNV'})

# There are 19,959 genes when you pool both datasets, but we only have both
# sgRNA abundance data and CNV data for 15,268 of them. Ultimately we don't have
# a choice but to work only with those 15,268, but it's worth seeing what the
# non-matching genes are to check if we've done something wrong. Uncomment the
# following lines to output a list of all genes.

# all_genes = pandas.merge(cnv_grouped_by_gene, guide_abundance, left_index=True,\
#                             right_index=True, how="outer")\
#                             .rename(columns={'VALUE': 'CNV'})
#
# # Write rows with null values in the CNV data to csv file
# all_genes[all_genes.isnull().CNV.iloc[:,0]].to_csv(open(\
#     "genes_found_only_in_guide-seq_data.csv", "w"))
#
# # Write rows with null values in the guide-seq data to csv file
# all_genes[all_genes.isnull().iloc[:,-1]].to_csv(open(\
#     "genes_found_only_in_CNV_data.csv", "w"))

# Anyway,  we want the inner join - genes that are in both sets:

# Get the median guide-seq data for each sample.
column_names = cnv_and_guides.columns.get_level_values(0)
sample_column_medians = [("CNV", "median"), ("CNV", "std")]

# Take every third item in the list that's not CNV data
sample_columns = [name for index, name in enumerate(column_names) if index in \
    range(0, len(column_names), 3) and name !="CNV"]
sample_column_medians.extend([(item, "median") for item in sample_columns])

# Trim the data to the medians of the guide-seq data for each gene
df = cnv_and_guides.loc[:,sample_column_medians]

# What is the actual mean CNV in the total dataset?
# Null hypothesis is that best linear fit of values will have y = 0x + mean(CNV)
print("Mean CNV in the dataset across guide-targeted regions is " + \
str(df.loc[:,("CNV", "median")].mean()))

# What guide-targeted regions of the genome exhibit the most extreme CNV?
df.loc[:,("CNV", "median")].sort_values().to_csv(open(\
    "ranked_cnv_by_gene.csv","w"))
print("Ranked list of CNV by guide-targeted genes output to" + \
    " ranked_cnv_by_gene.csv\n")
print("Top copy-number amplified genes\ngene\tfold enriched over control cells")
print(df.loc[:,("CNV", "median")].sort_values().tail(15).iloc[::-1])

print("\nTop copy-number depleted genes\ngene\tfold enriched over control" \
    + " cells")
print(df.loc[:,("CNV", "median")].sort_values().head(15))

# For this question, we're most interested in selection as a function of CNV.
# Whether this results in a growth advantage (leading to enrichment of a
# particular guide) or a growth disadvantage (leading to depletion of a
# particularl guide) is a peculiarity of the gene in question rather than the
# global effect of CNV on the ability of a guide to alter editing efficiency as
# a whole. In other words, the absolute factor of selection rather than the
# direction of selection is of most interest. So, I'm going to take the absolute
# value of the guide-seq data and fit that as a function of copy number.

# One hypothesis is that extra copies (high CNV) will result in less guide
# selection and that low copy would mean greater selective power of a given
# guide. That is, as CNV increases, guide selection index should decrease.

# The results of linear function fitting, parameters and plots are output to a
# PDF file "CNV vs Guide-Seq Abundance.pdf".

#http://www.wired.com/2011/01/linear-regression-with-pylab/
#https://stackoverflow.com/questions/2451264/creating-a-colormap-legend-in-
#matplotlib

font = {'size'   : 15}
matplotlib.rc('font', **font)
IndexSlice = pandas.IndexSlice
df.sortlevel(inplace=True, axis=1)
guide_enrichment = df.loc[:,IndexSlice[:,"median"]].iloc[:,1:].apply(abs)
figure, axes = plt.subplots(len(guide_enrichment.columns), 1, sharex=False, \
    sharey=False, figsize=(15,10*len(guide_enrichment.columns)))

cnv = df["CNV"]["median"]

for i, axes in enumerate(axes.ravel())  :
    y_axis = guide_enrichment.iloc[:,i]
    name = y_axis.name
    # Make a dummy plot to generate a color scale
    cbardummy = plt.hist2d(cnv, y_axis, bins = 200, cmap="Greys")
    figure.colorbar(cbardummy[3], ax=axes)

    # Generate the 2d histogram plot itself
    axes.hist2d(cnv, y_axis, bins = 200, cmap="Greys")

    # Fit a linear regression to the data
    slope, intercept, rvalue, pvalue, \
        stderr = scipy.stats.linregress(cnv, y_axis)

    # Display rsquared and slope of fit on the plot
    axes.annotate(str("p = " + str(pvalue) + "\nslope = " + str(slope)), \
                  xy = (0.05,0.95), color="k", va="top", ha="left", \
                  xycoords='axes fraction')

    # Add the fit line to the plot
    yp = np.polyval([slope,intercept], cnv)
    axes.plot(cnv, yp, color="r")

    # Set tile and axis text
    axes.set_title(name[0])
    axes.set_ylabel(\
        "Guide Selection Factor (absolute\nvs input plasmid representation)")
    axes.set_xlabel("Copy Number Variation Factor")

print("Plotted data and fit saved to CNV vs Guide Abundance.pdf.")
figure.savefig("CNV vs Guide Abundance.pdf", format='pdf')
