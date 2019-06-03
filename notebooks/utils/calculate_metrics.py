import numpy as np
import pandas as pd
import xarray as xr


from pathlib import Path
from .load import _load_experimental_ko_data


# Set up paths
data_path = Path("../data")
path_to_results = data_path / "simulation_results"

# Load ID dataframes
khod_idf = pd.read_csv(data_path / "khodayari_id.csv")
millard_idf = pd.read_csv(data_path / "millard_id.csv")
kurata_idf = pd.read_csv(data_path / "kurata_id.csv")
chassagnole_idf = pd.read_csv(data_path / "chassagnole_id.csv")



def _get_common_fluxes(data, author):
    """
    Find such fluxes that are common between all datasets
    """
    common_fluxes = set(khod_idf["BiGG ID"].unique()).intersection(
        set(millard_idf["BiGG ID"].unique()),
        set(kurata_idf["BiGG ID"].unique()),
        set(data.query("author == @author").BiGG_ID.unique()),
    )
    common_fluxes = {x for x in common_fluxes if pd.notna(x)}
    return common_fluxes


def process_data(data, author):
    """ 
    Subselect and check if the data is alright. Fixes some issues which can lead to numerical troubles.
    """
    selected_fluxes = _get_common_fluxes(data, author)
    selected_data = (
        data.query("BiGG_ID in @selected_fluxes")
        .groupby(["BiGG_ID", "sample_id", "author"])
        .median()
        .reset_index()
    )
    xdf = selected_data.set_index(["sample_id", "author", "BiGG_ID"]).to_xarray()

    abs_flux = xdf.flux
    # trim extremely small values
    xdf["normalized_flux"] = xdf.normalized_flux.where(
        abs(xdf.normalized_flux) > 1e-1, 0.0
    )
    xdf["normalized_flux"] = xdf.normalized_flux.where(~abs_flux.isnull(), np.NaN)
    return xdf


def relative_errors(xdata, author=None):
    """
    Calculates error metrics. Supposed to be run after check_data
    """
    xdf = xdata
    abs_flux = xdf.flux

    nm_flux = xdf.normalized_flux
    exp_flux = xdf.sel(author=author).normalized_flux

    # Mean absolute percent error (MAPE)
    xdf["relative_error"] = abs(nm_flux - exp_flux) / abs(nm_flux) * 100
    # If either predicted OR original are zeros then error is 100%
    xdf["relative_error"] = xdf["relative_error"].where(
        ((nm_flux != 0) & (exp_flux != 0)), 100
    )
    # If original data and predicted data were both zeros then error is zero as well
    xdf["relative_error"] = xdf["relative_error"].where(
        ~((nm_flux == 0) & (exp_flux == 0)), 0
    )
    # Put back NaNs where they were originally
    xdf["relative_error"] = xdf["relative_error"].where(~abs_flux.isnull(), np.NaN)

    # Calculate symmetric MAPE
    xdf["symm_relative_error"] = (
        100 * abs(nm_flux - exp_flux) / ((abs(nm_flux) + abs(exp_flux)) / 2)
    )
    return xdf


def summary_errors(xdata, author=None):
    """
    Calculates summary error per each sample_id for each model. Supposed to be run after check_data.
    Returns xarray DataSet with both normalized and non-normalized error.
    Normalized error is L2norm(pred-exp) divided by
    """

    def vector_norm(x, dim, ord=None):
        return xr.apply_ufunc(
            np.linalg.norm, x, input_core_dims=[[dim]], kwargs={"ord": ord, "axis": -1}
        )

    xdf = xdata
    nm_flux = xdf.normalized_flux
    exp_flux = xdf.sel(author=author).normalized_flux

    diff_norm = vector_norm(nm_flux - exp_flux, dim="BiGG_ID", ord=2)
    ishii_norm = vector_norm(exp_flux, dim="BiGG_ID", ord=2)
    x_norm_error = (
        (diff_norm / ishii_norm)
        .rename("normalized_error")
        .where(~diff_norm.isnull(), np.NaN)
    )
    diff_norm = diff_norm.rename("unnormalized_error")

    return xr.Dataset(
        data_vars={"normalized_error": x_norm_error, "unnormalized_error": diff_norm}
    )

