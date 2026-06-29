"""Pure stability classification kernel — testable without scanpy."""
import pandas as pd

from pipeline import config


def classify(df):
    short_max = config.STABILITY_CLASSES["short_max_h"]
    long_min = config.STABILITY_CLASSES["long_min_h"]
    out = df.copy()

    def _c(h):
        if pd.isna(h):
            return "unknown"
        if h < short_max:
            return "short"
        if h < long_min:
            return "medium"
        return "long"

    out["stability_class"] = out["half_life_h"].map(_c)
    return out
