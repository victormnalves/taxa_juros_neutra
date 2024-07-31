library(readxl)
sheets <- excel_sheets("census_pop_estimates_age_sex.xlsx")
for (sheet in sheets) {
    data <- read_excel("census_pop_estimates_age_sex.xlsx", sheet = sheet)
    write.csv(data, file = paste0("census_pop_estimates_age_sex_", sheet, ".csv"), na = "\"NaN\"")
}
