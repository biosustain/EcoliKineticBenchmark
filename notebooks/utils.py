# -*- coding: utf-8 -*-
import scipy
import scipy.io as sio
import numpy as np
import pandas as pd


def loadmat(filename):
    """
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """

    def _check_keys(d):
        """
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        """
        for key in d:
            if isinstance(d[key], sio.matlab.mio5_params.mat_struct):
                d[key] = _todict(d[key])
        return d

    def _has_struct(elem):
        """Determine if elem is an array and if any array item is a struct"""
        return isinstance(elem, np.ndarray) and any(
            isinstance(e, scipy.io.matlab.mio5_params.mat_struct) for e in elem
        )

    def _todict(matobj):
        """
        A recursive function which constructs from matobjects nested dictionaries
        """
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, sio.matlab.mio5_params.mat_struct):
                d[strg] = _todict(elem)
            elif _has_struct(elem):
                d[strg] = _tolist(elem)
            else:
                d[strg] = elem
        return d

    def _tolist(ndarray):
        """
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        """
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, sio.matlab.mio5_params.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif _has_struct(sub_elem):
                elem_list.append(_tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list

    data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def get_khodayari_kos():
    return {
        "fbaA": "result_cont_Delta_fbaAB.mat",
        "fbaB": "result_cont_Delta_fbaAB.mat",
        "fbp": "result_cont_Delta_fbp.mat",
        "gnd": "result_cont_Delta_gnd.mat",
        "pfkA": "result_cont_Delta_pfkAB.mat",
        "pfkB": "result_cont_Delta_pfkAB.mat",
        "pgi": "result_cont_Delta_pgi.mat",
        "pgl": "result_cont_Delta_pgl.mat",
        "ppsA": "result_cont_Delta_ppsA.mat",
        "pts": "result_cont_Delta_pts.mat",
        "pykA": "result_cont_Delta_pykA.mat",
        "pykF": "result_cont_Delta_pykF.mat",
        "rpe": "result_cont_Delta_rpe.mat",
        "rpiA": "result_cont_Delta_rpiAB.mat",
        "rpiB": "result_cont_Delta_rpiAB.mat",
        "sdhCD": "result_cont_Delta_sdhCD.mat",
        "sucA": "result_cont_Delta_sucAB.mat",
        "talA": "result_cont_Delta_talAB.mat",
        "tktA": "result_cont_Delta_tktAB.mat",
        "tktB": "result_cont_Delta_tktAB.mat",
        "tpi": "result_cont_Delta_tpi.mat",
        "zwf": "result_cont_Delta_zwf.mat",
        "WT": "result_cont_WT.mat",
    }

def get_millard_kos():
    return {
        "fbaA": "Millard_result_Delta_fba.csv",
        "fbaB": "Millard_result_Delta_fba.csv",
        "fbp": "Millard_result_Delta_fbp.csv",
        "gnd": "Millard_result_Delta_gnd.csv",
        "pfkA": "Millard_result_Delta_pfk.csv",
        "pfkB": "Millard_result_Delta_pfk.csv",
        "pgi": "Millard_result_Delta_pgi.csv",
        "pgl": "Millard_result_Delta_pgl.csv",
        "ppsA": "Millard_result_Delta_pps.csv",
        "pts": "Millard_result_Delta_pts.csv",
        "pykA": "Millard_result_Delta_pyk.csv",
        "pykF": "Millard_result_Delta_pyk.csv",
        "rpe": "Millard_result_Delta_rpe.csv",
        "rpiA": "Millard_result_Delta_rpi.csv",
        "rpiB": "Millard_result_Delta_rpi.csv",
        "sdhCD": "Millard_result_Delta_sdh.csv",
        "talAB": "Millard_result_Delta_tal.csv",
        "tkt1": "Millard_result_Delta_tkt1.csv",
        "tkt2": "Millard_result_Delta_tkt2.csv",
        "tpi": "Millard_result_Delta_tpi.csv",
        "zwf": "Millard_result_Delta_zwf.csv",
        "sucA": "Millard_result_Delta_aAkgdh.csv",
        "WT": "Millard_result_WT.csv",
    }

def get_kurata_kos():
    return {
        "fbaA": "result_cont_Delta_fbaB.mat",
        "fbaB": "result_cont_Delta_fbaB.mat",
        "fbp": "result_cont_Delta_fbp.mat",
        "gnd": "result_cont_Delta_gnd.mat",
        "gpmA": "result_cont_Delta_gpmA.mat",
        "pfkA": "result_cont_Delta_pfkA.mat",
        "pfkB": "result_cont_Delta_pfkB.mat",
        "pgi": "result_cont_Delta_pgi.mat",
        "ppc": "result_cont_Delta_ppc.mat",
        "pgl": "result_cont_Delta_pgl.mat",
        "ppsA": "result_cont_Delta_ppsA.mat",
        "pts": "result_cont_Delta_pts.mat",
        "pykA": "result_cont_Delta_pykA.mat",
        "pykF": "result_cont_Delta_pykF.mat",
        "rpe": "result_cont_Delta_rpe.mat",
        "rpiA": "result_cont_Delta_rpiA.mat",
        "rpiB": "result_cont_Delta_rpiB.mat",
        "sdhCD": "result_cont_Delta_sdhC.mat",
        "sucA": "result_cont_Delta_sucAC.mat",
        "talA": "result_cont_Delta_talA.mat",
        "talB": "result_cont_Delta_talB.mat",
        "tktA": "result_cont_Delta_tktA.mat",
        "tktB": "result_cont_Delta_tktB.mat",
        "tpi": "result_cont_Delta_tpi.mat",
        "zwf": "result_cont_Delta_zwf.mat",
        "glk": "result_cont_Delta_glk.mat",
        "WT": "result_cont_WT.mat",
    }

def get_kotte_kos():
    return {
        "fbaAB": "Kotte_result_Delta_fba.csv",
        "fbp": "Kotte_result_Delta_fbp.csv",
        "pfkA": "Kotte_result_Delta_pfk.csv",       
        "ppsA": "Kotte_result_Delta_pps.csv",
        "pts": "Kotte_result_Delta_pts.csv",        
        "pykF": "Kotte_result_Delta_pyk.csv",
    }

def get_chassagnole_kos():
    return {
        "fbaA": "Chassagnole_result_Delta_fba.csv",
        "fbaB": "Chassagnole_result_Delta_fba.csv",        
        "gnd": "Chassagnole_result_Delta_gnd.csv",
        "pfkA": "Chassagnole_result_Delta_pfk.csv",
        "pfkB": "Chassagnole_result_Delta_pfk.csv",
        "pgi": "Chassagnole_result_Delta_pgi.csv",
        #"pts": "Chassagnole_result_Delta_pts.csv",
        "pykAF": "Chassagnole_result_Delta_pyk.csv",
        "rpe": "Chassagnole_result_Delta_rpe.csv",
        "rpiA": "Chassagnole_result_Delta_rpi.csv",
        "rpiB": "Chassagnole_result_Delta_rpi.csv",        
        "talA": "Chassagnole_result_Delta_tal.csv",
        "talB": "Chassagnole_result_Delta_tal.csv",
        "tkt1": "Chassagnole_result_Delta_tkt1.csv",
        "tkt2": "Chassagnole_result_Delta_tkt2.csv",
        "tpi": "Chassagnole_result_Delta_tpi.csv",
        "zwf": "Chassagnole_result_Delta_zwf.csv",
    }

def get_khodayari_zwf():
    samples = ["dzwf", "WT", "zwf(15)"]
    return {k:f"khodayari_zwf_sens_{k}.mat" for k in samples}

def get_kurata_zwf():
    samples = ["dzwf", "WT", "zwf(15)"]
    return {k:f"kurata_zwf_sens_{k}.mat" for k in samples}

def get_millard_zwf():
        samples = ["dzwf", "WT", "zwf(15)"]
        return {k:f"millard_zwf_sens_{k}.csv" for k in samples}

def get_khodayari_pgi():
    samples = ["dpgi", "pgi(0)", "pgi(20)", "pgi(50)", "pgi(100)", "WT"]
    return {k:f"khodayari_pgi_sens_{k}.mat" for k in samples}

def get_kurata_pgi():
    samples = ["dpgi", "pgi(0)", "pgi(20)", "pgi(50)", "pgi(100)", "WT"]
    return {k:f"kurata_pgi_sens_{k}.mat" for k in samples}  

def get_millard_pgi():
        samples = ["dpgi", "pgi(0)", "pgi(20)", "pgi(50)", "pgi(100)", "WT"]
        return {k:f"millard_pgi_sens_{k}.csv" for k in samples}

def get_khodayari_eno():
    samples = ["eno(0)", "eno(50)", "eno(200)", "eno(500)", "WT"]
    return {k:f"khodayari_eno_sens_{k}.mat" for k in samples}

def get_kurata_eno():
    samples = ["eno(0)", "eno(50)", "eno(200)", "eno(500)", "WT"]
    return {k:f"kurata_eno_sens_{k}.mat" for k in samples} 

def get_millard_eno():
        samples = ["eno(0)", "eno(50)", "eno(200)", "eno(500)", "WT"]
        return {k:f"millard_eno_sens_{k}.csv" for k in samples}

def relative_error(df, exp_name = None, measurement_type = "normalized_flux", threshold = 1000):
    """
    Calculate relative error in percentages compared with experimental dataset. 
    Supposed to be applied while grouping by BiGG_ID and sample_id.
    Example - all_data.groupby(["BiGG_ID", "sample_id"]).apply(relative_error, exp_name="Nicolas").reset_index().drop("level_2", axis = 1)
    Params:
    :df - group that would be passed by lambda
    :exp_name - value of "author" column that would be chosen as experimental data to be compared against
    :measurement_type - string that specifies on which kind of values to compute relative error
    :threshold - value that would be the maximum possible error
    """
    # find the experimental value 
    if exp_name is None:
        raise ValueError("please specify the name of experimental dataset")
    exp_value = df.loc[df["author"] == exp_name, measurement_type].values[0]
    if np.isnan(exp_value):
        exp_value = 0.0
    value = abs(abs(df[measurement_type] - exp_value) / df[measurement_type]) * 100
    
    # Check for errors which are more than threshold and change them to :trim_value
    if threshold is not None:
        value.loc[value > threshold] = threshold 
    
    # create new column for relative error
    new_df = df.assign(relative_error = value)
    new_df = new_df[["author", "flux", "normalized_flux", "relative_error"]]
    return new_df

def calculate_norm(x, exp_data):
    """
    Calculate ||x - exp_data || / ||exp_data|| in such way that NAs are discarded in both numerator and denominator
    Params:
    :x - pd.Series of measurement_type that was declared in normalized_error 
    :exp_data - pd.Sereis of experimental data with the same measurement_type    
    """    
    diff = x - exp_data
    nonna_ids = diff[~diff.isna()].index
    value = np.linalg.norm(diff.loc[nonna_ids]) / np.linalg.norm(exp_data.loc[nonna_ids])
    return value


def normalized_error(df, exp_name = None, measurement_type = "normalized_flux", threshold = 1000):
    """
    Calculate relative error in percentages compared with experimental dataset. 
    Supposed to be applied while grouping by BiGG_ID and sample_id.
    Example - selected_data.groupby(["sample_id"]).apply(normalized_error, exp_name="Ishii").reset_index()
    Params:
    :df - group that would be passed by lambda
    :exp_name - value of "author" column that would be chosen as experimental data to be compared against
    :measurement_type - string that specifies on which kind of values to compute normalized error
    """
    # find the experimental value 
    df = df.set_index("BiGG_ID")
    if exp_name is None:
        raise ValueError("please specify the name of experimental dataset")
    exp_values = df.loc[df["author"] == exp_name, measurement_type]

    # find the error estimate tgat is norm of (vexp - vsim) / norm of vexp. 
    # https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003580    
    
    #diff = df.groupby("author").apply(lambda x: np.linalg.norm((exp_values - x[measurement_type]).dropna()) / np.linalg.norm(exp_values))
    diff = df.groupby("author").apply(lambda x: calculate_norm(x = x[measurement_type], exp_data = exp_values))

    return diff

# Set of routines to load data for various models 
# and present them as pandas dataframe

def load_khodayari(sample_names, load_path, id_df, files = None):
    """ Will return single dataframe with columns:
    - author=Khodayari, 
    - sample_id corresponding to relevant sample
    - flux with correspoding value in original units
    - ID showing original id
    - BiGG ID with BiGG identifier
    params:
    :sample_names - list of sample names or "all"
    :load_path - Path object where the samples are located
    :id_df - pd.DataFrame with conversion Model ID -> BiGG ID
    :files - dict where key is sample name and value is filename
    """

    if files is None:
        raise ValueError("files dictionary is not specified")

    if type(sample_names) is not list:
        sample_names = [sample_names]

    # check mismatch with sample names in inputs
    unknown_ids = [elem for elem in sample_names if elem not in files.keys()]

    if sample_names == ["all"]:
        sample_names = files.keys()
    elif unknown_ids:
        error_msg = ", ".join(unknown_ids)
        raise ValueError(f"Unable to find relevant data for {error_msg}")

    res_df = pd.DataFrame()
    for sample_id in sample_names:
        file_name = files[sample_id]
        data = loadmat(load_path / file_name)
        data_shape = data["Vnet"].shape
        print(
            f"Loaded data file for sample {sample_id} which has flux matrix of {data_shape}"
        )

        # data['Vnet'][:,data_shape[1]-1] is the last column of integration, ideally it should be closer to steady state
        # khod_rxn_ids[455] is the index of 'Biomass' flux, the last flux id
        df = pd.DataFrame(
            {
                "flux": data["Vnet"][0:457, data_shape[1] - 1],
                "ID": id_df["ID"],
                "BiGG_ID": id_df["BiGG ID"],
            }
        )
        df = df.assign(author="Khodayari", sample_id=sample_id)
        
        # normalize to the glucose uptake
        glucose_uptake = df[df['ID'] == 'EX_glc(e)']['flux'].values[0]
        df = df.assign(normalized_flux = lambda x: x.flux * 100 / glucose_uptake)
        
        # update
        res_df = pd.concat([res_df, df])

    return res_df

def load_kurata(sample_names, load_path, id_df, files = None):
    """ Will return single dataframe with columns:
    - author=Kurata, 
    - sample_id corresponding to relevant sample
    - flux with correspoding value in original units
    - ID with BiGG identifier
    params:
    :sample_names - list of sample names or "all"
    :load_path - Path object where the samples are located
    :id_df - pd.DataFrame with conversion Model ID -> BiGG ID
    :files - dict where key is sample name and value is filename
      """
    if files is None:
        raise ValueError("files dictionary is not specified")

    if type(sample_names) is not list:
        sample_names = [sample_names]

    # check mismatch with sample names in inputs
    unknown_ids = [elem for elem in sample_names if elem not in files.keys()]

    if sample_names == ["all"]:
        sample_names = files.keys()
    elif unknown_ids:
        error_msg = ", ".join(unknown_ids)
        raise ValueError(f"Unable to find relevant data for {error_msg}")

    res_df = pd.DataFrame()
    for sample_id in sample_names:
        file_name = files[sample_id]
        data = loadmat(load_path / file_name)
        data_shape = data["FLUX"].shape
        print(
            f"Loaded data file for sample {sample_id} which has flux matrix of {data_shape}"
        )

        # this weird construction helps to deal with multiple IDs corresponding to Gapdh reaction
        # which in Kuratas model is a sum of gapA, tpiA, gpmA or gpmM, eno, pgk
        # resulting id would be GAPD.
        kurata_ids = id_df[["ID", "BiGG ID"]].drop_duplicates(subset="ID")

        df = pd.DataFrame(
            {
                "flux": data["FLUX"][2100,],
                "ID": kurata_ids["ID"].values,
                "BiGG_ID": kurata_ids["BiGG ID"].values,
            }
        )
        df = df.assign(author="Kurata", sample_id=sample_id)

        # add results from Gapdh reaction
        flux_val = df[df["BiGG_ID"] == "GAPD"]["flux"].values[0]
        add_fluxes_df = pd.DataFrame(
            {
                "flux": flux_val,
                "ID": "Gapdh",
                "BiGG_ID": ["TPI", "PGM", "ENO", "PGK"],
                "author": "Kurata",
                "sample_id": sample_id,
            }
        )
        df = df.append(add_fluxes_df, sort=False)
        
        # calculate normalized fluxes with respect to Glucose consumption
        glucose_uptake = df[df['ID'] == 'vPts4']['flux'].values[0] + df[df['ID'] == 'vNonpts']['flux'].values[0]
        df = df.assign(normalized_flux = lambda x: x.flux * 100 / glucose_uptake)
        
        # update
        res_df = pd.concat([res_df, df])
    
    return res_df

def load_millard(sample_names, load_path, id_df, files = None):
    """ Will return single dataframe with columns:
    - author=Millard, 
    - sample_id corresponding to relevant sample
    - flux with correspoding value in original units
    - ID with BiGG identifier
    params:
    :sample_names - list of sample names or "all"
    :load_path - Path object where the samples are located
    :id_df - pd.DataFrame with conversion Model ID -> BiGG ID
    :files - dict where key is sample name and value is filename
      """
    if files is None:
        raise ValueError("files dictionary is not specified")

    if type(sample_names) is not list:
        sample_names = [sample_names]

    # check mismatch with sample names in inputs
    unknown_ids = [elem for elem in sample_names if elem not in files.keys()]

    if sample_names == ["all"]:
        sample_names = files.keys()
    elif unknown_ids:
        error_msg = ", ".join(unknown_ids)
        raise ValueError(f"Unable to find relevant data for {error_msg}")

    res_df = pd.DataFrame()
    for sample_id in sample_names:
        file_name = files[sample_id]

        data = pd.read_csv(load_path / file_name)
        data_shape = data["ID"].shape
        print(
            f"Loaded data file for sample {sample_id} which has flux matrix of {data_shape}"
        )

        df = pd.merge(
            left=data.drop("Unnamed: 0", axis=1),
            right=id_df[["ID", "BiGG ID"]],
            how="left",
            on="ID",
        )

        df = df.assign(author="Millard", sample_id=sample_id)
        df = df.rename({"BiGG ID": "BiGG_ID", "Value": "flux"}, axis=1)
        
        # Set MDH reaction to be the difference between MQO and MDH flux
        mdh_flux = df.loc[df.ID == 'MQO', 'flux'].values[0] - df.loc[df.ID == 'MDH', 'flux'].values[0]
        df.loc[df.BiGG_ID == 'MDH', 'flux'] = mdh_flux  
        
        # Calculate normalized fluxes
        glucose_uptake = df[df['ID'] == 'XCH_GLC']['flux'].values[0]
        df = df.assign(normalized_flux = lambda x: x.flux * 100 / glucose_uptake)
        

        
        # update
        res_df = pd.concat([res_df, df])

    return res_df

def load_kotte(sample_names, load_path, id_df, files = None):

    """ Will return single dataframe with columns:
    - author=Millard, 
    - sample_id corresponding to relevant sample
    - flux with correspoding value in original units
    - ID with BiGG identifier
    params:
    :sample_names - list of sample names or "all"
    :load_path - Path object where the samples are located
    :id_df - pd.DataFrame with conversion Model ID -> BiGG ID
    :files - dict where key is sample name and value is filename
      """
    if files is None:
        raise ValueError("files dictionary is not specified")

    if type(sample_names) is not list:
        sample_names = [sample_names]

    # check mismatch with sample names in inputs
    unknown_ids = [elem for elem in sample_names if elem not in files.keys()]

    if sample_names == ["all"]:
        sample_names = files.keys()
    elif unknown_ids:
        error_msg = ", ".join(unknown_ids)
        raise ValueError(f"Unable to find relevant data for {error_msg}")

    res_df = pd.DataFrame()
    for sample_id in sample_names:
        file_name = files[sample_id]

        data = pd.read_csv(load_path / file_name)
        data_shape = data["ID"].shape
        print(
            f"Loaded data file for sample {sample_id} which has flux matrix of {data_shape}"
        )

        df = pd.merge(
            left=data.drop("Unnamed: 0", axis=1),
            right=id_df[["ID", "BiGG ID"]],
            how="left",
            on="ID",
        )

        df = df.assign(author="Kotte", sample_id=sample_id)
        df = df.rename({"BiGG ID": "BiGG_ID", "Value": "flux"}, axis=1)
        # update
        res_df = pd.concat([res_df, df])

    return res_df

def load_chassagnole(sample_names, load_path, id_df, files = None):
    """ Will return single dataframe with columns:
    - author=Millard, 
    - sample_id corresponding to relevant sample
    - flux with correspoding value in original units
    - ID with BiGG identifier
    params:
    :sample_names - list of sample names or "all"
    :load_path - Path object where the samples are located
    :id_df - pd.DataFrame with conversion Model ID -> BiGG ID
    :files - dict where key is sample name and value is filename
      """
    if files is None:
        raise ValueError("files dictionary is not specified")

    if type(sample_names) is not list:
        sample_names = [sample_names]

    # check mismatch with sample names in inputs
    unknown_ids = [elem for elem in sample_names if elem not in files.keys()]

    if sample_names == ["all"]:
        sample_names = files.keys()
    elif unknown_ids:
        error_msg = ", ".join(unknown_ids)
        raise ValueError(f"Unable to find relevant data for {error_msg}")

    res_df = pd.DataFrame()
    for sample_id in sample_names:
        file_name = files[sample_id]

        data = pd.read_csv(load_path / file_name)
        data_shape = data["ID"].shape
        print(
            f"Loaded data file for sample {sample_id} which has flux matrix of {data_shape}"
        )

        df = pd.merge(
            left=data.drop("Unnamed: 0", axis=1),
            right=id_df[["ID", "BiGG ID"]],
            how="left",
            on="ID",
        )

        df = df.assign(author="Chassagnole", sample_id=sample_id)
        df = df.rename({"BiGG ID": "BiGG_ID", "Value": "flux"}, axis=1)
        
        # Calculate normalized fluxes
        glucose_uptake = df[df['ID'] == 'vPTS']['flux'].values[0]
        df = df.assign(normalized_flux = lambda x: x.flux * 100 / glucose_uptake)

        # update
        res_df = pd.concat([res_df, df])

    return res_df