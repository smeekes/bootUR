# Macroeconomic Time Series

Macroeconomic data from Eurostat on GDP, consumption, inflation and
unemployment for Belgium, Germany, France, the Netherlands and the
United Kingdom.

## Usage

``` r
MacroTS
```

## Format

A time series object containing 20 macroeconomic seasonally adjusted
time series, quarterly observed from 1992-2019 for Belgium (BE), Germany
(DE), France (FR), the Netherlands (NL) and the United Kingdom (UK).

- `GDP_BE`:

  Gross domestic product at market prices (index, 2015=100) for Belgium.

- `GDP_DE`:

  Gross domestic product at market prices (index, 2015=100) for Germany.

- `GDP_FR`:

  Gross domestic product at market prices (index, 2015=100) for France.

- `GDP_NL`:

  Gross domestic product at market prices (index, 2015=100) for the
  Netherlands.

- `GDP_UK`:

  Gross domestic product at market prices (index, 2015=100) for the
  United Kingdom.

- `CONS_BE`:

  Final consumption expenditure (index, 2015=100) for Belgium.

- `CONS_DE`:

  Final consumption expenditure (index, 2015=100) for Germany.

- `CONS_FR`:

  Final consumption expenditure (index, 2015=100) for France.

- `CONS_NL`:

  Final consumption expenditure (index, 2015=100) for the Netherlands.

- `CONS_UK`:

  Final consumption expenditure (index, 2015=100) for the United
  Kingdom.

- `HICP_BE`:

  Harmonised Indices of Consumer Prices (annual rate of change,
  2015=100) for Belgium.

- `HICP_DE`:

  Harmonised Indices of Consumer Prices (annual rate of change,
  2015=100) for Germany.

- `HICP_FR`:

  Harmonised Indices of Consumer Prices (annual rate of change,
  2015=100) for France.

- `HICP_N`:

  Harmonised Indices of Consumer Prices (annual rate of change,
  2015=100) for the Netherlands.

- `HICP_UK`:

  Harmonised Indices of Consumer Prices (annual rate of change,
  2015=100) for the United Kingdom.

- `UR_BE`:

  Unemployment rate (percentage of the active population) for Belgium.

- `UR_DE`:

  Unemployment rate (percentage of the active population) for Germany.

- `UR_FR`:

  Unemployment rate (percentage of the active population) for France.

- `UR_NL`:

  Unemployment rate (percentage of the active population) for the
  Netherlands.

- `UR_UK`:

  Unemployment rate (percentage of the active population) for the United
  Kingdom.

## Source

https://ec.europa.eu/eurostat/data/database

## Note

- Unemployment rates are seasonally but not calendar adjusted, all other
  series are both seasonally and calendar adjusted.

- Quarterly inflation rates are sampled from Eurostat's monthly series
  with annual rates of change as the final month of the respective
  quarter.

- The unemployment rate for France excludes overseas territories
  ('France continental' in the Eurostat database).
