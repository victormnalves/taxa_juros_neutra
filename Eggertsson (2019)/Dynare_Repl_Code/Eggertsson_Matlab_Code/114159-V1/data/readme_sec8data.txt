Source data, Section 8 for Eggertsson, Mehrotra, & Robbins, "A Model of Secular Stagnation: Theory and Quantitative Evaluation"


DATA SOURCES

This folder contains the underlying data used to calibrate the quantitative secular stagnation model. The model requires a data series for US government debt (including state and local government debt), age-specific survival rates, total factor productivity, total fertility rate, and the age distribution of the population.

1. TOTAL FACTOR PRODUCTIVITY: (Column A, jr_tfp.xlsx) HP-filtered annual utilization adjusted total factor productivity (dtfp_util) from John G. Fernald, “A Quarterly, Utilization Adjusted Series on Total Factor Productivity.” Available from <https://www.frbsf.org/economic-research/indicators-data/total-factor-productivity-tfp/>. This series is modified in two ways - the initial level of TFP in 1970 is set at the average from 1948-1974 and the final value of TFP is set at the average from 1975-1993. The optimistic series simply sets the terminal level of TFP to its 1948-1974 average.

2. GOVERNMENT DEBT: (Column B, jr_gdebt.xlsx) 5-year moving average of federal plus state and local government debt (% of GDP). Federal debt is FRED series, Gross Federal Debt as a Percent of Gross Domestic Product (GFDGDPA188S). State and local debt is imputed from 2012-13 data on state and local debt obtained from US Census Bureau’s Survey of State and Local Governments. Available from <https://www.census.gov/govs/local/index.html>. The 2012 and 2013 average of state/local government debt is assumed to be a constant multiple of gross federal debt. A five-year centered moving average is taken to obtain the series in Column B.

3. TOTAL FERTILITY RATE: (Column A, jr_tfr.xlsx) 5-year centered moving average of 0.5 x US total fertility rate (births per woman). US total fertility rate from 1960-2015 obtained from World Bank World Development Indicator, series SP.DYN.TFRT.IN. Available from <https://data.worldbank.org/indicator/SP.DYN.TFRT.IN?locations=US>. US total fertility rate from 1945-1960 obtained from Centers for Disease Control Vital Statistics. Available from Table 1-7: <https://www.cdc.gov/nchs/data/statab/t991x07.pdf>.

4. SURVIVAL RATES: (Columns A-BD, jr_survival.csv) One-year age-specific survival probabilities by year (rows, 1970-2007) and columns are from age 25-80. Survival rates are obtained as 1-mortality rate within age-specific bins: 25-34, 35-44, 45-54, 55-64, 65-74, 75-84. For each year, survival rates by age are interpolated within age bins. Mortality (per 100,000 persons) within age bins by year reported by Centers for Disease Control, Vital Statistics. Available from <https://www.cdc.gov/nchs/nvss/mortality/gmwk23r.htm> and <https://www.cdc.gov/nchs/nvss/mortality/hist290.htm>.

5. US POPULATION PYRAMID: (Columns A-B, us_census_pop.xlsx) US population by age 25-81 in 1970 and 2015. Full population pyramid estimates by year available from US Census Bureau: <https://www.census.gov/data/tables/time-series/demo/popest/pre-1980-national.html>. 


