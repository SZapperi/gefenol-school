"""
GEO2pandas: Gene Expression Omnibus to Pandas.

GEO2pandas is a python package to retrieve
data and metadata from the Gene Expression Omnibus
and place it into pandas dataframes.

The interaction with GEO is handled with GEOparse:
LINK TO GEOPARSE

GEO2pandas adds some convenient wrappers and sanity checks.

Francesc Font-Clos
Center for Complexity & Biosystems
April 2018
"""

# python 2/3 compatibility
from __future__ import division
from __future__ import print_function
from builtins import str


import numpy as np
import pandas as pd
import GEOparse
import copy
import os


def get_meta_from_gse(gse, metadata_field="characteristics_ch1"):
    """
    Construct a metadata dataframe from a `gse` object and metadata field.

    Parameters
    -----------
    gse : GSE object
        The main GSE object holding all data.
    metadata_field: str or list[str]
        Metadata field(s) from which to parse metadata.

    Returns
    -------
    meta : pandas.DataFrame
        A dataframe holding the metadata, i.e., information on the samples of
        the experiment (group, condition, ...).

    """
    if type(metadata_field) is str:
        fields = [metadata_field]
    else:
        for x in metadata_field:
            assert type(x) is str
        fields = metadata_field

    for x in fields:
        assert x in list(gse.gsms.values())[0].metadata.keys()

    tmp = []
    for field in fields:
        t_meta = pd.DataFrame()
        for gsm_name, gsm in gse.gsms.items():
            series = _create_series_from_gsm(gsm_name, gsm,
                                             metadata_field=field)
            t_meta = t_meta.append(series)
        tmp.append(t_meta.T)

    meta = pd.concat(tmp).T

    # get rid of spaces in column names
    meta.columns = [x[1:] if x.find(" ") == 0 else x for x in meta.columns]
    meta.columns = [x.replace(" ", "_") for x in meta.columns]

    meta.sort_index(axis=0, inplace=True)
    return meta


def get_expr_from_gse(gse, gpl_gene_col=None, gpl=None,
                      value="VALUE", gpl_annot_col="ID",
                      ):
    """
    Construct expression dataframe from a `gse` object.

    Parameters
    -----------
    gse : GSE object
        The main GSE object holding all data.
    gpl_gene_col : str
        Field used to annotate. If not provided,
        tries to get Entrez labels, then GenBank.
    gpl : GEOparse.GPL
        Annotation platform. If not provided,
        first one is used.

    Returns
    -------
    expr : pandas.DataFrame
        A dataframe holding the raw gene expression data.

    See also
    --------
    GEO2pandas.clean_and_check_expr

    Other parameters
    ----------------
    value : str
        Column that holds the values. Defaults to "VALUE".
    gpl_annot_col : str
        Column used to link annotations. Defaults to "ID".

    """
    # chose gpl if not provided
    if gpl is None:
        gpl = list(gse.gpls.values())[0]
        if len(gse.gpls.keys()) != 1:
            print("# WARNING: this GSE has several gpls.")
            print("# using first one on the list: %s" % gpl.name)

    # choose annotation column
    if gpl_gene_col is None:
        for entrez in ["Entrez", "ENTREZ", "EntrezID",
                       "ENTREZ_GENE_ID", "Entrez_Gene_ID"]:
            if entrez in gpl.columns.index:
                gpl_gene_col = entrez
                print("")
                print("# WARNING: chose %s as gene name column" % gpl_gene_col)
                print("# Use gpl_gene_col keyword to choose yourself")
                break
        # print warning if there were no matches for entrez
        if gpl_gene_col is None:
            print("# WARNING: found no columns matching to Entrez ID")
            print("# Now looking for Genebank columns")

    # look for genebank only if no entrez matches!
    if gpl_gene_col is None:
        for genebank in ["GeneBank", "GB_ACC", "GB_LIST"]:
            if genebank in gpl.columns.index:
                gpl_gene_col = genebank
                print
                print("# WARNING: chose %s as gene name column" % gpl_gene_col)
                print("# Use gpl_gene_col keyword to choose yourself")
                break
        # print warning if there were no matches for Genebank either
        if gpl_gene_col is None:
            print("# WARNING: found no columns matching to Genebank")
            print("# You need to use the gpl_gene_col keyword")
            return None

    # pivot samples
    tmp_expr = gse.pivot_samples(value)

    # annotate expression dataframe
    tmp = (gpl.table[[gpl_annot_col, gpl_gene_col]])
    annot_dict = dict(zip(tmp.values.T[0], tmp.values.T[1]))

    def mymapper(x):
        return annot_dict[x] if x in annot_dict.keys() else x

    tmp_expr.index = tmp_expr.index.map(mymapper)
    expr = tmp_expr.T
    expr.columns.name = gpl_gene_col
    expr.sort_index(axis=0, inplace=True)
    expr.sort_index(axis=1, inplace=True)
    return expr


def get_gse(gse=None, filepath=None):
    """
    Get a GSE object.

    Wrapper around GEOparse.get_GEO that also checks if data
    is available with small caps. Example: if gse1234.soft.gz
    is present, but not GSE1234.soft.gz, GEO2parse redownloads it.

    Parameters
    ----------
    gse: str
        GEO accession code of the form 'GSE12345' or 'gse12345'.
    filepath : str
        Path to local copy of GSE, if present.

    Returns
    -------
    gse : GEOparse.GSE
        The GSE object.

    """
    if gse is None:
        print("# You MUST specify the GSE accession code!")
        print("")
    else:
        assert(type(gse) == str)
        _gse = gse.upper()

    # check if file is available
    filepath = None
    for filestr in ["GSE"+_gse[3:]+".soft.gz", "gse"+_gse[3:]+".soft.gz"]:
        if os.path.isfile(filestr):
            filepath = filestr

    if filepath is not None:
        gse = GEOparse.get_GEO(filepath=filepath)
    else:
        gse = GEOparse.get_GEO(geo=_gse)
    return gse


def clean_and_check_expr(expr=None,
                         average_out_duplicate_cols=True,
                         drop_ambiguous_cols=True,
                         drop_nan_cols=True,
                         verbose=False,
                         allow_nans=False,
                         allow_duplicate_rows=False,
                         remove_suffixes_from_genenames=True,
                         check_log=True):
    """
    Sanity checks for data retrieved from GEO.

    This function performs several sanity checks designed mostly
    for gene expression data retrieved from GEO via
    `GEO2pandas.get_expr_from_gse()`.

    Parameters
    ----------
    expr : pandas.DataFrame
        The gene expression data.
    average_out_duplicate_columns : bool
        Columns with the same name are averaged out.
    drop_ambiguous_cols : bool
        Drops all columns containing '/' or ',' characters in their names.
    drop_nan_cols : bool
        Drop any column names 'nan', 'NaN', np.nan.
    allow_nans : bool
        Allow some matrix entries to be nan, i.e., missing data.
        Otherwise it deletes all columns that contain missin values.
    allow_duplicate_rows : bool
        Allows non-unique labels in rows.
    remove_suffixes_from_genenames : bool
        Removes suffixes (such as .1) from genenames
    check_log : bool
        Tries to guess if data is in log space.

    Returns
    -------
    _expr : pandas.DataFramse
        A sanitized copy of `expr`.

    """
    if expr is None:
        print("# You MUST pass me a pandas dataframe holding gene"
              "expression values (samples in rows, genes in cols)")
        return None

    # define which strings mean NaN
    nanvalues = ["nan", "Nan", "NaN", "n/a"]

    # define what strings are a sign of ambiguity
    ambig_strings = ["/", "|", ","]

    # get a copy to work on
    _expr = copy.deepcopy(expr)

    if verbose:
        print("# initial shape")
        print(_expr.shape)
        print("")

    # cast column names to unicode
    if remove_suffixes_from_genenames:
        _expr.columns = [str(x).split(".")[0] for x in _expr.columns]
    else:
        _expr.columns = [str(x) for x in _expr.columns]

    # drop columns called nan
    if drop_nan_cols:
        _expr = _expr.T[[x is not np.nan for x in _expr.columns]].T
        for nan in nanvalues:
            _expr = _expr.T[[x.find(nan) == -1 for x in _expr.columns]].T

        if verbose:
            print("# after deleting nan columns")
            print(_expr.shape)
            print("")

    # drop ambiguous columns
    if drop_ambiguous_cols:
        for ambig in ambig_strings:
            _expr = _expr.T[[x.find(ambig) == -1 for x in _expr.columns]].T
        if verbose:
            print("# after dropping ambiguous columns")
            print(_expr.shape)
            print("")

    # drop columns with nans
    # check if there are nans indeed...
    dims = _expr.shape
    allvals = _expr.values.reshape(dims[0]*dims[1])
    l1 = len([x for x in allvals if (x in nanvalues)])
    l2 = len([x for x in allvals if x is np.nan])
    l3 = l1+l2
    # ... and map them to np.nan
    if l3 > 0:
        _expr = _expr.applymap(_nicenan)

    if allow_nans:
        if verbose:
            print("# you set allownan=True, up to you mate!")
            print("")
    else:
        _expr = _expr.T[_expr.isnull().sum(axis=0) == 0].T
        if verbose:
            print("# after dropping columns with nans")
            print(_expr.shape)
            print("")

    # average out duplicate columns
    if average_out_duplicate_cols:
        _expr = _expr.groupby(level=0, axis=1).mean()
        if verbose:
            print("# after averaging out duplicate columns")
            print(_expr.shape)
            print("")

    # warn about duplicate rows
    if allow_duplicate_rows is False:
        if np.unique(_expr.index.values).shape[0] < _expr.shape[0]:
            raise ValueError("There are non-unique labels in rows."
                             "Set allow_duplicate_rows to True to"
                             "avoid this error.")

    # make sure we are on logspace
    if check_log:
        assert _expr.max().max() < 100

    return _expr


def _create_series_from_gsm(gsm_name, gsm,
                            metadata_field="characteristics_ch1"):
    metadata = gsm.metadata[metadata_field]

    if len(metadata) == 1:
        num_commas = len(metadata[0].split(","))-1
        num_dospunts = len(metadata[0].split(":"))-1
        # many fields
        if (num_dospunts > 1) and (num_commas > 0) and\
           (num_dospunts == num_commas + 1):
            metadata = [x for x in metadata[0].split(",")]
        # a single field without name
        if num_dospunts == 0 and num_commas == 0:
            metadata = [metadata_field+":"+metadata[0]]

    return pd.Series(dict(zip(
        [x.split(":")[0] for x in metadata],
        [x.split(":")[1] for x in metadata])), name=gsm_name)


def _nicenan(x):
    nanvalues = ["nan", "Nan", "NaN", "n/a"]

    assert (x in nanvalues) or (type(x) is float) or (type(x) is np.float64)

    if x in nanvalues:
        return np.nan
    elif x is np.inf:
        return np.nan
    else:
        return x
