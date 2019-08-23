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
    data = xdf.to_dataframe().reset_index()
    if sample_id:
        source = data.query(f"sample_id == '{sample_id}'")
        title = f"Heatmap for sample_id {sample_id}"
    else:
        source = data.query(f"author == '{author}'")
        title = f"Heatmap for author {author}"

    base = (
        alt.Chart(source, title=title)
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
            tooltip=["sample_id", "BiGG_ID", "relative_error", "normalized_flux"],
        )
    )
    if sample_id:
        chart = base.encode(x=alt.X("author:N"))
    else:
        chart = base.encode(x=alt.X("sample_id:N"))

    return chart.configure(invalidValues=None)


def _get_colors():
    domain = [
        "Khodayari",
        "Millard",
        "Kurata",
        "iML1515",
        "Chassagnole",
        "Exp_iML1515",
        "Exp_ECC2",
        "Long",
    ]
    range_ = [
        "#1b9e77",
        "#d95f02",
        "#7570b3",
        "#e7298a",
        "#66a61e",
        "#e6ab02",
        "#a6761d",
        "#666666",
    ]
    return domain, range_


def summary_chart(norm_error, author=None, title=None, sort_list=None):
    """
    Draws summary chart where one sample is represented by one dot. 
    Supposed to be used on data where errors are "summarized" 
    """
    source = norm_error.to_dataframe().reset_index().query(f"author != '{author}'")

    selector = alt.selection_single(empty="all", fields=["sample_id"])
    color = alt.condition(
        selector, alt.Color("author:N", sort=sort_list), alt.ColorValue("lightgray")
    )
    opacity = alt.condition(selector, alt.OpacityValue(0.4), alt.OpacityValue(1.0))
    size = alt.condition(selector, alt.SizeValue(100), alt.SizeValue(40))

    base = (
        alt.Chart(source, title=title)
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
        .configure_axis(labelFontSize=20, titleFontSize=20)
        .configure_legend(labelFontSize=20, titleFontSize=20)
        .configure_header(titleFontSize=24)
        .configure(invalidValues=None)
    )

    return chart


def jitter_summary_chart(
    norm_error,
    author=None,
    title=None,
    sort_list=None,
    opacity=False,
    color_scheme="category10",
):
    source = norm_error.to_dataframe().reset_index()
    if type(author) == list:
        source = source.query(f"author not in @author")
    elif type(author) == str:
        source = source.query(f"author != '{author}'")

    selector = alt.selection_single(empty="none", fields=["sample_id"])
    if opacity:
        opacity = alt.condition(selector, alt.OpacityValue(1.0), alt.OpacityValue(0.5))
    else:
        opacity = alt.OpacityValue(1.0)

    size = alt.condition(selector, alt.SizeValue(150), alt.SizeValue(60))

    # modify that for changes in coloring. Make sure that domain is correct
    domain, range_ = _get_colors()

    stripplot = (
        alt.Chart(source, width=100, height=600)
        .mark_circle()
        .encode(
            x=alt.X(
                "jitter:Q",
                title=None,
                axis=alt.Axis(values=[0], ticks=True, grid=False, labels=False),
                scale=alt.Scale(),
            ),
            y=alt.Y("normalized_error:Q", title="Normalized error"),
            color=alt.Color(
                "author:N", legend=None, scale=alt.Scale(domain=domain, range=range_)
            ),
            column=alt.Column(
                "author:N",
                header=alt.Header(
                    labelAngle=-90,
                    titleOrient="top",
                    labelOrient="bottom",
                    labelAlign="right",
                    labelPadding=3,
                    labelFontSize=20,
                    title=title,
                ),
                sort=sort_list,
            ),
            size=size,
            opacity=opacity,
            tooltip=["author", "sample_id", "normalized_error"],
        )
        .transform_calculate(
            # Generate Gaussian jitter with a Box-Muller transform
            jitter="sqrt(-2*log(random()))*cos(2*PI*random())"
        )
        .configure_facet(spacing=5)
        .configure_view(stroke=None)
        .add_selection(selector)
        .transform_filter("datum.normalized_error !== null")
    )

    return stripplot.configure_axis(
        labelFontSize=20, titleFontSize=20
    ).configure_header(titleFontSize=24)


def boxplot(
    norm_error,
    author=None,
    title=None,
    sort_list=None,
    opacity=False,
    color_scheme="category10",
):
    source = norm_error.to_dataframe().reset_index()
    if type(author) == list:
        source = source.query(f"author not in @author")
    elif type(author) == str:
        source = source.query(f"author != '{author}'")

    selector = alt.selection_single(empty="none", fields=["sample_id"])
    if opacity:
        opacity = alt.condition(selector, alt.OpacityValue(1.0), alt.OpacityValue(0.5))
    else:
        opacity = alt.OpacityValue(1.0)

    size = alt.condition(selector, alt.SizeValue(150), alt.SizeValue(60))

    # modify that for changes in coloring. Make sure that domain is correct
    domain, range_ = _get_colors()

    stripplot = (
        alt.Chart(width=100, height=600)
        .mark_circle()
        .encode(
            x=alt.X(
                "jitter:Q",
                title=None,
                axis=alt.Axis(values=[0], ticks=True, grid=False, labels=False),
                scale=alt.Scale(),
            ),
            y=alt.Y("normalized_error:Q", title="Normalized error"),
            color=alt.Color(
                "author:N", legend=None, scale=alt.Scale(domain=domain, range=range_)
            ),
            size=size,
            opacity=opacity,
            tooltip=["author", "sample_id", "normalized_error"],
        )
        .add_selection(selector)
    )

    boxplot = (
        alt.Chart()
        .mark_boxplot(size=95)
        .encode(
            y=alt.Y("normalized_error:Q", title="Normalized error"),
            color=alt.Color(
                "author:N", legend=None, scale=alt.Scale(domain=domain, range=range_)
            ),
            opacity=alt.OpacityValue(0.7),
        )
    )

    layer = (
        alt.layer(boxplot, stripplot, data=source)
        .transform_calculate(
            # Generate Gaussian jitter with a Box-Muller transform
            jitter="sqrt(-2*log(random()))*cos(2*PI*random())"
        )
        .facet(
            column=alt.Column(
                "author:N",
                header=alt.Header(
                    labelAngle=-90,
                    titleOrient="top",
                    labelOrient="bottom",
                    labelAlign="right",
                    labelPadding=3,
                    labelFontSize=20,
                    title=title,
                ),
                sort=sort_list,
            )
        )
        .configure_facet(spacing=5)
        .configure_view(stroke=None)
        .transform_filter("datum.normalized_error !== null")
    )

    return layer.configure_axis(labelFontSize=20, titleFontSize=20).configure_header(
        titleFontSize=24
    )
