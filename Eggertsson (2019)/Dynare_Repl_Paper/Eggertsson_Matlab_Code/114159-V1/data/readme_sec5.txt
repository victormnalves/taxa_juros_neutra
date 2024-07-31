Instructions for replication of Figure 5 and Appendix H in Eggertsson, Mehrotra, & Robbins, "A Model of Secular Stagnation: Theory and Quantitative Evaluation"


DATA

The underlying data in Figure 5 (blue starred lines) can be found in Figure 5.xlsx which contains separate spreadsheets for the US, Japan, and the Eurozone. The underlying data sources along with websites for downloading updated series are documented at the bottom of the spreadsheet.

MODEL

The key moments from the model used to plot the impulse response of the 3 period secular stagnation model are highlighted in light blue in the spreadsheets of Figure 5.xlsx. These moments are obtained from the matlab programs stag_US.m, stag_JP.m, and stag_EU.m. Each of these m-files numerically calibrate and solve the log-linearized three equation secular stagnation model with hysteresis. The details of the calibration are discussed in Appendix H. These m-files call the programs gx_hx.m and ir.m available from Stephanie Schmitt-Grohe and Martin Uribe <http://www.columbia.edu/~mu2166/1st_order/1st_order.htm>. The program computes the impulse response for output, inflation, and total factor productivity for each country.

The calculations in Figure 5.xlsx take these underlying moments for log-linear deviations from steady state and recover projections for output and inflation from the model.

