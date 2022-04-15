###############################################################################
#
# DESCRIPTION
#
# This code will parse a GSI diagnostic winds (diag_conv_uv) file to identify
# duplicate pairs of observations as defined by GSI/src/setupw.f90. When
# duplicate observations are discovered, observation errors can be inflated
# based on the time-difference between the pairs, which can affect assimilation
# penalty terms and statistics in the gsistat file.
#
# The purpose of this program is to explicitly find which observations may be
# flagged as duplicates in the GSI, so that an experiment and control cycle can
# be compared. It is possible, for example, for an experiment to make changes
# to one observation type/subtype, and then for penalty term statistics to
# change for some other type/subtype. This can happen if the observations that
# were changed in the experiment are now changing the distribution of
# duplications found by setupw.f90. So if you see the penalty terms change on
# the first assimilation cycle for an observation type/subtype you are not
# intending to change, you can compare an experiment and control diagnostic
# file with this program to see if the penalty term changes could be the result
# of a new distribution of duplications (as opposed to an error in the way the
# experiment is designed).
#
# YAML FORMAT
#
# dupsearch_yaml.yaml has two key entries:
#  1) 'diag': values defining the diag_conv file used
#         diag_filedir: full-path to diag_conv file (including "/" at the end)
#         diag_obvar: diag_conv file variable-type ('t','uv','q', or 'ps')
#         diag_type: diag_conv file reference-type ('ges' or 'anl')
#         diag_date: diag_conv file date (YYYYMMDDHH, or %Y%m%d%H format)
#
# 2) 'target': values defining type and subtype filters for target data
#         target_typ: filter target data by type, list of types, or None
#         target_styp: filter target data by subtype, list of subtypes, or None
#
# type and subtype filters for 'target' can either specify a single value, a
# list of multiple values (in which case an observation matching any element in
# list will be retained), or None (in which case the filter is bypassed
# allowing any observation)
###############################################################################

def parse_yaml(yaml_file):
    import yaml
    # YAML entries come in two types:
    # key='diag' : provides diag_conv file directory path, ob-var (t, uv, q,
    #              ps), diag-type (ges or anl), and analysis date (YYYYMMDDHH)
    # key='target' : provides type and subtype filters for target obs to search
    #                for duplicates (or None, for no filter on type or subtype)
    #
    with open(yaml_file, 'r') as stream:
        try:
            parsed_yaml = yaml.safe_load(stream)
        except yaml.YAMLError as YAMLError:
            parsed_yaml = None
            print(f'YAMLError: {YAMLError}')
        if parsed_yaml is not None:
            # Extract 'diag' data
            try:
                diag_filedir = parsed_yaml['diag']['filedir']
                diag_obvar = parsed_yaml['diag']['obvar']
                diag_type = parsed_yaml['diag']['type']
                diag_date = parsed_yaml['diag']['date']
            except KeyError as MissingDiagError:
                diag_filedir = None
                diag_obvar = None
                diag_type = None
                diag_date = None
                print(f'MissingDiagError: {MissingDiagError}')
            # Extract 'target' data
            try:
                target_typ = None if parsed_yaml['target']['typ'] == 'None'\
                    else parsed_yaml['target']['typ']
                target_styp = None if parsed_yaml['target']['styp'] == 'None'\
                    else parsed_yaml['target']['styp']
            except KeyError as MissingTargetError:
                target_typ = None
                target_styp = None
                print(f'MissingTargetError: {MissingTargetError}')
    return (diag_filedir, diag_obvar, diag_type, diag_date, target_typ,
            target_styp)


if __name__ == "__main__":
    from pyGSI.diags import Conventional
    import numpy as np
    import pandas as pd

    # Parsy YAML file for diag and target data
    (diag_filedir, diag_obvar, diag_type, diag_date, target_typ,
        target_styp) = parse_yaml('dupsearch_yaml.yaml')
    # Define diagnostic filename
    diagnosticFile = diag_filedir + 'diag_conv_' + diag_obvar + '_' +\
        diag_type+'.' + diag_date + '.nc4'
    # Load diagnosticFile, filter to include only assimilated (use==1) obs
    diag = Conventional(diagnosticFile)
    s = diag.data_df.loc[(slice(None), slice(None), slice(None), slice(None),
                          slice(None), slice(None), slice(None)), :]
    s = s[s.index.isin([1], level='Analysis_Use_Flag')]
    obs_use = []
    obs_pre = []
    obs_typ = []
    obs_styp = []
    obs_err2 = []  # adjusted obs error, not input error or final error
    obs_lat = np.asarray(s['latitude'])
    obs_lon = np.asarray(s['longitude'])
    obs_err2 = np.asarray(s['errinv_adjust'])
    # Error is read as inverse-error, recover error by inverting
    obs_err2 = obs_err2**-1
    for x in s.index:
        obs_typ.append(x[2])
        obs_styp.append(x[3])
        obs_pre.append(x[4])
        obs_use.append(x[6])
    obs_typ = np.asarray(obs_typ, dtype=np.float32)
    obs_styp = np.asarray(obs_styp, dtype=np.float32)
    obs_pre = np.asarray(obs_pre, dtype=np.float32)
    obs_use = np.asarray(obs_use, dtype=np.float32)
    n_obs = np.size(obs_lat)
    # Define targeted obs as a subset of all obs
    all_idx = np.arange(n_obs)  # vector of all observation indices
    # Ob-type filter
    if target_typ is not None:
        # Assert target_type as a list if it is not a list
        target_typ = [target_typ] if type(target_typ) != list else target_typ
        # Find indices of all obs with obs_typ in target_type list
        tgt_typ_idx = all_idx[np.isin(obs_typ, target_typ)]
    else:
        tgt_typ_idx = all_idx
    # Ob-subtype filter
    if target_styp is not None:
        # Assert target_subtype as a list if it is not a list
        target_styp = [target_styp] if type(target_styp) != list\
            else target_styp
        tgt_styp_idx = all_idx[np.isin(obs_styp, target_styp)]
    else:
        tgt_styp_idx = all_idx
    # Combine filters
    tgt_obs_idx = np.intersect1d(tgt_typ_idx, tgt_styp_idx)
    n_tgt = np.size(tgt_obs_idx)
    tgt_dup_idx = []
    obs_dup_idxs = []
    for i in range(n_tgt):
        xlat = obs_lat[tgt_obs_idx[i]]
        xlon = obs_lon[tgt_obs_idx[i]]
        xpre = obs_pre[tgt_obs_idx[i]]
        xerr2 = obs_err2[tgt_obs_idx[i]]
        if xerr2 < 1000.:
            dups = np.where((obs_lat == xlat) &
                            (obs_lon == xlon) &
                            (obs_pre == xpre) &
                            (obs_err2 < 1000.))[0]
            dups = dups[np.where(~np.isin(dups, tgt_obs_idx[i]))]
            if np.size(dups) > 0:
                tgt_dup_idx.append(tgt_obs_idx[i])
                obs_dup_idxs.append(tuple(dups))
    n_dup = len(tgt_dup_idx)
    # Report
    print(f'There are {n_dup} duplicate observations among selected targets:')
    for i in range(n_dup):
        print(f'Target       {tgt_dup_idx[i]}: Type {obs_typ[tgt_dup_idx[i]]},\
 SubType {obs_styp[tgt_dup_idx[i]]}')
        for d in obs_dup_idxs[i]:
            print(f'   Duplicate {d}: Type {obs_typ[d]}, SubType \
{obs_styp[d]}')
