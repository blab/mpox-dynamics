import pandas as pd


def clean_owid_cases(
    raw_cases, filter_locations, min_date=None, max_date=None
):
    """
    Filter raw OWID cases to locations and time period of interest.
    """
    rc = raw_cases.copy()
    rc = rc[["location", "date", "new_cases"]]

    # Filter to locatins of interest
    rc = rc[rc.location.isin(filter_locations)]
    rc = rc.rename(columns={"new_cases": "cases"})

    # Define and filter time period for analysis
    min_date = rc.date.min() if min_date is None else min_date
    max_date = rc.date.max() if max_date is None else max_date
    rc = rc[(rc.date > min_date) & (rc.date < max_date)]

    return rc


OWID_MONKEYPOX_CASES_URL = "https://raw.githubusercontent.com/owid/monkeypox/main/owid-monkeypox-data.csv"
FILTER_LOCATIONS = [
    "United States",
    "United Kingdom",
    "Spain",
    "Portugal",
    "Germany",
    "France",
    "Colombia",
    "World",
]
MIN_DATE, MAX_DATE = "2022-05-01", "2022-10-31"

if __name__ == "__main__":
    raw_cases = pd.read_csv(OWID_MONKEYPOX_CASES_URL)
    cleaned_cases = clean_owid_cases(
        raw_cases, FILTER_LOCATIONS, MIN_DATE, MAX_DATE
    )
    cleaned_cases.to_csv("./monkeypox-cases-counts.tsv", sep="\t", index=False)
