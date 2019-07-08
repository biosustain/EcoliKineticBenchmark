# -*- coding: utf-8 -*-
import pandas as pd
import xarray as xr

import altair as alt
from altair.expr import datum


def heatmap(xdf, author=None, sample_id=None):
    """
    Function to create relative error heatmaps showing missing values as grey, 
    can work with heatmap per 1 sample or per 1 author
    """
    data = xdf["relative_error"].to_dataframe().reset_index()
    if sample_id:
        data = data.query(f"sample_id == '{sample_id}'")
        title = f"Heatmap for sample_id {sample_id}"
    else:
        data = data.query(f"author == '{author}'")
        title = f"Heatmap for author {author}"

    base = (
        alt.Chart(data, title=title)
        .mark_rect()
        .encode(
            y=alt.Y("BiGG_ID:N"),
            color=alt.condition(
                "datum.relative_error === null",
                alt.ColorValue("lightgrey"),
                alt.Color(
                    "relative_error",
                    scale=alt.Scale(
                        domain=[0, 25, 50, 75, 100, 200, 300],
                        type="threshold",
                        scheme="greenblue",
                    ),
                    title="relative error (clipped at 300)",
                ),
            ),
            tooltip=["sample_id", "BiGG_ID", "relative_error"],
        )
    )
    if sample_id:
        chart = base.encode(x=alt.X("author:N"))
    else:
        chart = base.encode(x=alt.X("sample_id:N"))

    return chart.configure(invalidValues=None)


def summary_chart(norm_error, author=None, title=None, sort_list=None):
    """
    Draws summary chart where one sample is represented by one dot. 
    Supposed to be used on data where errors are "summarized" 
    """
    selector = alt.selection_single(empty="all", fields=["sample_id"])
    color = alt.condition(
        selector, alt.Color("author:N", sort=sort_list), alt.ColorValue("lightgray")
    )

    opacity = alt.condition(selector, alt.OpacityValue(0.4), alt.OpacityValue(1.0))

    size = alt.condition(selector, alt.SizeValue(100), alt.SizeValue(40))

    base = (
        alt.Chart(
            norm_error.to_dataframe().reset_index().query(f"author != '{author}'"),
            title=title,
        )
        .mark_circle()
        .encode(
            y=alt.Y("normalized_error", title="Normalized error"),
            tooltip=["author", "sample_id", "normalized_error"],
            x=alt.X("author", sort=sort_list),
            size=size,
            opacity=opacity,
        )
    )

    errors = (
        base.encode(color=color)
        .add_selection(selector)
        .transform_filter("datum.normalized_error !== null")
    )

    # na_vals = base.encode(color=alt.value("lightgrey")).transform_filter("datum.normalized_error === null")

    chart = (
        errors.properties(width=700, height=600)
        .configure_axis(labelFontSize=24, titleFontSize=24)
        .configure_legend(labelFontSize=16, titleFontSize=20)
        .configure(invalidValues=None)
    )

    return chart
