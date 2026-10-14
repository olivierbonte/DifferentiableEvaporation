# %% Imports
import pandas as pd
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def knn_gapfill(
    df: pd.DataFrame, target: str, predictors: list[str], k=5
) -> tuple[pd.Series, pd.Series]:
    """
    Fill the gaps of `target` in `df` with a K-nearest neighbour regression on `predictors`.

    Follows `fill_netrad` of FluxDataKit
    (https://github.com/geco-bern/FluxDataKit/blob/43fb25847bb99ab8b02abf943fa82933c5a5315b/R/helper_functions.R#L61):
    the predictors are centred and scaled, and a gap is filled with the mean of the `k`
    nearest neighbours, trained on all rows where target and predictors are present.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing the target and predictors.
    target : str
        Name of the target column to fill.
    predictors : list[str]
        List of predictor column names to use for the KNN regression.
    k : int, optional
        Number of nearest neighbours to use for the regression. Default is 5.

    Returns
    -------
    filled : pd.Series
        The filled target column.
    to_fill : pd.Series
        A boolean mask of the rows that were filled.
    """
    complete_predictors = df[predictors].notna().all(axis=1)
    train = complete_predictors & df[target].notna()
    to_fill = complete_predictors & df[target].isna()
    filled = df[target].copy()
    if to_fill.any():
        knn = make_pipeline(StandardScaler(), KNeighborsRegressor(n_neighbors=k))
        knn.fit(df.loc[train, predictors], df.loc[train, target])
        filled.loc[to_fill] = knn.predict(df.loc[to_fill, predictors])
    return filled, pd.Series(to_fill, index=df.index)
