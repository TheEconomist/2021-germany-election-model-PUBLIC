#! /usr/bin/env python3

# %% imports and setup
from datetime import date
import requests
import html
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup as Soup

URL_BASE = "https://www.wahlrecht.de/umfragen"
POLLSTERS = [
    "allensbach",
    "emnid",
    "forsa",
    "politbarometer",
    "gms",
    "dimap",
    "insa",
    "yougov",
]

#%% testing
poll = "insa"
page = requests.get("{url}/{poll}.htm".format(url=URL_BASE, poll=poll))
soup = Soup(page.text, features="html.parser")
table = soup.select_one("table.wilko")
# Remove the footers
[e.extract() for e in table.select(".foot")]
data = pd.read_html(str(table), parse_dates=[0])[0]
data = (
    data.rename(columns={"Unnamed: 0": "Datum"})
    .dropna(axis="columns", how="all")
    .dropna(axis="rows", how="all")[:-2]
)
num_cols = data.columns.drop(["Datum", "Befragte", "Zeitraum"])
data[num_cols] = (
    data[num_cols]
    .apply(
        lambda x: x.apply(
            lambda y: y if type(y) is float else y.rstrip("%").replace(",", ".")
        )
    )
    .apply(pd.to_numeric, errors="coerce")
    / 100
)
data["Befragte"] = (
    data["Befragte"]
    .apply(lambda x: x.replace(".", ""))
    .apply(pd.to_numeric, errors="coerce")
)
data["poll"] = poll
data.head()

# %%
def parsePoll(poll):
    # poll = POLLSTERS[0]
    print(poll)
    page = requests.get("{url}/{poll}.htm".format(url=URL_BASE, poll=poll))
    soup = Soup(page.text, features="html.parser")
    table = soup.select_one("table.wilko")
    # Remove the footers
    [e.extract() for e in table.select(".foot")]
    data = pd.read_html(str(table), parse_dates=[0])[0]
    data = (
        data.rename(columns={"Unnamed: 0": "Datum"})
        .dropna(axis="columns", how="all")
        .dropna(axis="rows", how="all")[:-2]
    )
    num_cols = data.columns.drop(["Datum", "Befragte", "Zeitraum"])
    data[num_cols] = (
        data[num_cols]
        .apply(
            lambda x: x.apply(
                lambda y: y if type(y) is float else y.rstrip("%").replace(",", ".")
            )
        )
        .apply(pd.to_numeric, errors="coerce")
        / 100
    )
    data["Befragte"] = (
        data["Befragte"]
        .apply(lambda x: x.replace(".", ""))
        .map(lambda x: x.split(" • ")[-1])
        .apply(pd.to_numeric, errors="coerce")
    )
    data["poll"] = poll
    # data[num_cols] = data[num_cols].apply(pd.to_numeric, errors='coerce')
    return data


# %%
cols_to_use = [
    "Datum",
    "CDU/CSU",
    "SPD",
    "GRÜNE",
    "FDP",
    "LINKE",
    "AfD",
    "Sonstige",
    "Befragte",
    "Zeitraum",
    "poll",
]
all_polls = pd.concat([parsePoll(p) for p in POLLSTERS], axis="rows")[cols_to_use]

# %% Process and consolidate
all_polls["Datum"] = pd.to_datetime(all_polls["Datum"], dayfirst=True)
all_polls = all_polls[
    all_polls.Datum > "2017-09-24"
]  # throw out polls from before the last election
all_polls = all_polls[all_polls["Zeitraum"] != "Bundestagswahl"]
all_polls["year"] = all_polls["Datum"].dt.to_period("Y")
all_polls["start"] = all_polls.Zeitraum.str.split("–").map(lambda x: x[0])
all_polls["end"] = all_polls.Zeitraum.str.split("–").map(lambda x: x[-1])
all_polls["start"] = pd.to_datetime(
    all_polls.apply(
        lambda row: row.start
        + str(
            row.year - 1
            if (row.Datum.month == 1) & (row.start[3:5] == "12")
            else row.year
        ),
        axis=1,
    ),
    format="%d.%m.%Y",
    errors="coerce",
)
all_polls["end"] = pd.to_datetime(
    all_polls.apply(
        lambda row: row.end
        + str(
            row.year - 1
            if (row.Datum.month == 1) & (row.end[3:5] == "12")
            else row.year
        ),
        axis=1,
    ),
    format="%d.%m.%Y",
    errors="coerce",
)
melt_polls = all_polls.melt(
    id_vars=("Datum", "poll", "year", "start", "end", "Befragte", "Zeitraum"),
    var_name="party",
    value_name="share",
).rename(columns={"Datum": "date", "Befragte": "n"})
melt_polls = melt_polls[melt_polls.columns.drop(["year", "Zeitraum"])]
melt_polls


# %% Save out
melt_polls.to_csv("./output-data/all-polls-2021.csv")
# %%
