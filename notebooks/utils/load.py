import numpy as np
import pandas as pd
import xarray as xr

from pathlib import Path
from contextlib import redirect_stdout
import io

from .utils import (
    get_khodayari_kos,
    get_millard_kos,
    get_kurata_kos,
    get_chassagnole_kos,
    get_khodayari_dilutions,
    get_kurata_dilutions,
    get_millard_dilutions,
    get_khodayari_zwf,
    get_kurata_zwf,
    get_millard_zwf,
    get_khodayari_pgi,
    get_kurata_pgi,
    get_millard_pgi,
    get_khodayari_eno,
    get_kurata_eno,
    get_millard_eno,
    load_khodayari,
    load_kurata,
    load_millard,
    load_chassagnole,
    loadmat,
)


# Set up paths
data_path = Path("../data")
path_to_results = data_path / "simulation_results"

# Load ID dataframes
khod_idf = pd.read_csv(data_path / "khodayari_id.csv")
millard_idf = pd.read_csv(data_path / "millard_id.csv")
kurata_idf = pd.read_csv(data_path / "kurata_id.csv")
chassagnole_idf = pd.read_csv(data_path / "chassagnole_id.csv")


def _load_kinetic_ko_sims():
    """
    Load the data from kinetic models simulations
    """
    khodayari_results = load_khodayari(
        sample_names="all",
        load_path=(path_to_results / "Khodayari"),
        id_df=khod_idf,
        files=get_khodayari_kos(),
    )

    kurata_results = load_kurata(
        sample_names="all",
        load_path=(path_to_results / "Kurata"),
        id_df=kurata_idf,
        files=get_kurata_kos(),
    )

    millard_results = load_millard(
        sample_names="all",
        load_path=(path_to_results / "Millard"),
        id_df=millard_idf,
        files=get_millard_kos(),
    )

    chassagnole_results = load_chassagnole(
        sample_names="all",
        load_path=(path_to_results / "Chassagnole"),
        id_df=chassagnole_idf,
        files=get_chassagnole_kos(),
    )

    simulation_results = pd.concat(
        [khodayari_results, kurata_results, millard_results, chassagnole_results],
        sort=False,
    )
    return simulation_results


def _load_experimental_ko_data():
    """
    Load Ishii data
    """
    df = pd.read_csv("../data/datasets/ishii2007_tidy.csv")
    # this regexp matches deletions starting with d like dpgi
    df["sample_id"] = df.Genotype.str.extract(r"d(\w+)")
    df.loc[df.Genotype == "WT", "sample_id"] = "WT"

    df = df.assign(author="Ishii")
    df = df.rename(
        {
            "Measurement_ID": "BiGG_ID",
            "Original_Value": "normalized_flux",
            "Value": "flux",
            "Original_ID": "ID",
        },
        axis=1,
    )
    df = df[df["Measurement_Type"] == "flux"]
    df.loc[df["BiGG_ID"] == "PYKF", "BiGG_ID"] = "PYK"

    df = df[["flux", "ID", "BiGG_ID", "author", "sample_id", "normalized_flux"]]

    """
    douglas_flux_data = pd.read_csv("../../../DataAnalysis/DouglasKineticData/data/flux_data_processed.csv", index_col=0)
    douglas_sample_names = {
        "Evo04": "WT",
        "Evo04gnd": "gnd",
        "Evo04pgi": "pgi",
        "Evo04sdhCB": "sdh",
        "Evo04tpiA": "tpi",
    }
    douglas_exp_df = douglas_flux_data.query(
    "sample_name in @douglas_sample_names.keys()"
    )
    douglas_exp_df["sample_id"] = douglas_exp_df.sample_name.apply(
        lambda x: douglas_sample_names[x]
    )
    douglas_exp_df["author"] = "McCloskey"
    douglas_exp_df = douglas_exp_df.rename(
        {"rxn_id": "BiGG_ID", "sampling_median": "flux"}, axis=1
    ).drop(["sampling_var", "sampling_min", "sampling_max", "sample_name"], axis=1)
    x_doug = douglas_exp_df.set_index(["author", "sample_id", "BiGG_ID"]).to_xarray()
    x_doug['normalized_flux'] = 100*(x_doug.flux / x_doug.sel(BiGG_ID="GLCptspp").flux)
    mccloskey_results = x_doug.to_dataframe().reset_index()
    mccloskey_results.head()
    """

    return df


def _load_cobra_ko_sims():
    """
    Load simulations from iML1515, ECC2 and iML1515 and ECC2 conditioned on experimental data
    """
    iml_results = pd.read_csv(
        path_to_results / "COBRA" / "iML1515" / "knockouts_all.csv", index_col=0
    )
    ecc_results = pd.read_csv(
        path_to_results / "COBRA" / "ECC2" / "knockouts_all.csv", index_col=0
    )

    exp_iml_results = pd.read_csv(
        path_to_results / "COBRA" / "Exp_iML1515" / "knockouts_all.csv", index_col=0
    )
    exp_ecc_results = pd.read_csv(
        path_to_results / "COBRA" / "Exp_ECC2" / "knockouts_all.csv", index_col=0
    )
    
    df = pd.concat([iml_results, exp_iml_results])
    # Fix direction to match experimental data
    # Only for iML1515
    # fix PGM direction
    df.loc[df["ID"] == "PGM", "flux"] = -1*df.loc[df["ID"] == "PGM", "flux"].values[0]
    df.loc[df["ID"] == "PGM", "normalized_flux"] = -1*df.loc[df["ID"] == "PGM", "normalized_flux"].values[0]
    # fix PGK direction
    df.loc[df["ID"] == "PGK", "flux"] = -1*df.loc[df["ID"] == "PGK", "flux"].values[0]
    df.loc[df["ID"] == "PGK", "normalized_flux"] = -1*df.loc[df["ID"] == "PGK", "normalized_flux"].values[0]

    df = pd.concat([df, ecc_results, exp_ecc_results])
    # For both iML1515 and ECC2
    # fix RPI direction
    df.loc[df["ID"] == "RPI", "flux"] = -1*df.loc[df["ID"] == "RPI", "flux"].values[0]
    df.loc[df["ID"] == "RPI", "normalized_flux"] = -1*df.loc[df["ID"] == "RPI", "normalized_flux"].values[0]

    return df


def load_ko_data():
    """
    Load all simulations,
    """
    with io.StringIO() as buf, redirect_stdout(buf):
        simulation_data = _load_kinetic_ko_sims()
        cobra_data = _load_cobra_ko_sims()
        exp_data = _load_experimental_ko_data()
        file_info = buf.getvalue()
    return pd.concat([simulation_data, cobra_data, exp_data], sort=False), file_info


def _load_kinetic_dilution_sims():
    khodayari_dil = load_khodayari(
        sample_names="all",
        load_path=(path_to_results / "Khodayari" / "dilutions"),
        id_df=khod_idf,
        files=get_khodayari_dilutions(),
    )
    kurata_dil = load_kurata(
        sample_names="all",
        load_path=(path_to_results / "Kurata" / "dilutions"),
        id_df=kurata_idf,
        files=get_kurata_dilutions(),
    )
    millard_dil = load_millard(
        sample_names="all",
        load_path=(path_to_results / "Millard" / "dilutions"),
        id_df=millard_idf,
        files=get_millard_dilutions(),
    )
    return pd.concat([khodayari_dil, kurata_dil, millard_dil], sort=False)


def _load_experimental_dilution_data():
    yao_df = pd.read_csv(data_path / "datasets" / "yao2011_tidy.csv")
    consumption_rates = yao_df.query('Measurement_Type == "consumption_rate"')
    yao_fluxes = yao_df.query('Measurement_Type == "flux"')

    def normalize_to_uptake(group):
        consumption_rate = consumption_rates.loc[
            consumption_rates.Dilution == group.name, "Value"
        ].values[0]
        print(f"Consumption rate for D {group.name} is {consumption_rate}")
        group = group.assign(normalized_flux=lambda x: x.Value / consumption_rate * 100)
        return group

    df = (
        yao_fluxes.groupby("Dilution").apply(normalize_to_uptake).reset_index(drop=True)
    )

    df = df.assign(author="Yao")
    df = df.rename(
        {
            "Measurement_ID": "BiGG_ID",
            "Value": "flux",
            "Original_ID": "ID",
            "Dilution": "sample_id",
        },
        axis=1,
    )
    df = df[df["Measurement_Type"] == "flux"]

    df = df[["flux", "ID", "BiGG_ID", "author", "sample_id", "normalized_flux"]]
    df.sample_id = df.sample_id.apply(str)
    return df


def load_dilution_data():
    with io.StringIO() as buf, redirect_stdout(buf):
        simulation_data = _load_kinetic_dilution_sims()
        exp_data = _load_experimental_dilution_data()
        file_info = buf.getvalue()
    return pd.concat([simulation_data, exp_data], sort=False), file_info


def _load_experimental_sensitivity_data():
    """
    Load experimental results and return a tuple of 3 dataframes each corresponding to 
    genes zwf, pgi and eno
    """

    """
    Nicloas, 2007 data, for zwf knockout
    """
    df = pd.read_csv("../data/datasets/nicolas2007_tidy.csv")

    df = df.assign(author="Nicolas")
    df = df.rename(
        {
            "Measurement_ID": "BiGG_ID",
            "Original_Value": "normalized_flux",
            "Value": "flux",
            "Original_ID": "ID",
            "Genotype": "sample_id",
        },
        axis=1,
    )
    df = df[df["Measurement_Type"] == "flux"]

    df = df[["flux", "ID", "BiGG_ID", "author", "sample_id", "normalized_flux"]]
    exp_results_zwf = df

    """
    Usui, 2012 data for pgi and eno data
    """
    df = pd.read_csv("../data/datasets/usui2012_tidy.csv")
    df = df.assign(author="Usui")
    df = df.rename(
        {
            "Measurement_ID": "BiGG_ID",
            "Original_Value": "normalized_flux",
            "Value": "flux",
            "Original_ID": "ID",
            "Genotype": "sample_id",
        },
        axis=1,
    )
    df = df[df["Measurement_Type"] == "flux"]

    df = df[["flux", "ID", "BiGG_ID", "author", "sample_id", "normalized_flux"]]
    exp_results_pgi = df
    exp_results_eno = df
    return (exp_results_zwf, exp_results_pgi, exp_results_eno)


def _load_kinetic_sensitivity_sims():
    khodayari_zwf = load_khodayari(
        sample_names="all",
        load_path=(path_to_results / "Khodayari" / "zwf_sensitivity"),
        id_df=khod_idf,
        files=get_khodayari_zwf(),
    )
    khodayari_pgi = load_khodayari(
        sample_names="all",
        load_path=(path_to_results / "Khodayari" / "pgi_sensitivity"),
        id_df=khod_idf,
        files=get_khodayari_pgi(),
    )
    khodayari_eno = load_khodayari(
        sample_names="all",
        load_path=(path_to_results / "Khodayari" / "eno_sensitivity"),
        id_df=khod_idf,
        files=get_khodayari_eno(),
    )
    kurata_zwf = load_kurata(
        sample_names="all",
        load_path=(path_to_results / "Kurata" / "zwf_sensitivity"),
        id_df=kurata_idf,
        files=get_kurata_zwf(),
    )
    kurata_pgi = load_kurata(
        sample_names="all",
        load_path=(path_to_results / "Kurata" / "pgi_sensitivity"),
        id_df=kurata_idf,
        files=get_kurata_pgi(),
    )
    kurata_eno = load_kurata(
        sample_names="all",
        load_path=(path_to_results / "Kurata" / "eno_sensitivity"),
        id_df=kurata_idf,
        files=get_kurata_eno(),
    )
    millard_zwf = load_millard(
        sample_names="all",
        load_path=(path_to_results / "Millard" / "zwf_sensitivity"),
        id_df=millard_idf,
        files=get_millard_zwf(),
    )
    millard_pgi = load_millard(
        sample_names="all",
        load_path=(path_to_results / "Millard" / "pgi_sensitivity"),
        id_df=millard_idf,
        files=get_millard_pgi(),
    )
    millard_eno = load_millard(
        sample_names="all",
        load_path=(path_to_results / "Millard" / "eno_sensitivity"),
        id_df=millard_idf,
        files=get_millard_eno(),
    )

    simulation_zwf = pd.concat([khodayari_zwf, kurata_zwf, millard_zwf], sort=False)
    simulation_pgi = pd.concat([khodayari_pgi, kurata_pgi, millard_pgi], sort=False)
    simulation_eno = pd.concat([khodayari_eno, kurata_eno, millard_eno], sort=False)
    return (simulation_zwf, simulation_pgi, simulation_eno)


def load_sensitivity_data():
    with io.StringIO() as buf, redirect_stdout(buf):
        simulation_data_zwf, simulation_data_pgi, simulation_data_eno = (
            _load_kinetic_sensitivity_sims()
        )
        exp_data_zwf, exp_data_pgi, exp_data_eno = _load_experimental_sensitivity_data()
        file_info = buf.getvalue()
    return (
        (
            pd.concat([simulation_data_zwf, exp_data_zwf], sort=False),
            pd.concat([simulation_data_pgi, exp_data_pgi], sort=False),
            pd.concat([simulation_data_eno, exp_data_eno], sort=False),
        ),
        file_info,
    )

