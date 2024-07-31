# This function is used to read the demographic data.
using DelimitedFiles
function read_data(path_csv::String)
    data = data_t(readdlm(string(path_csv, "marriageAgeDiff_Q.csv"), ','),
              readdlm(string(path_csv, "interp_death_rate_1900_2220_Q.csv"), ','),
              readdlm(string(path_csv, "births_Q_annualized.csv"), ','),
              readdlm(string(path_csv, "population_1899_Q.csv"), ','),
              readdlm(string(path_csv, "share_births_mothers_Q.csv"), ','),
              readdlm(string(path_csv, "netmigration_Q.csv"), ','),
              readdlm(string(path_csv, "parent_child_1899_Q.csv"), ','),
              readdlm(string(path_csv, "dependents_1899_Q.csv"), ','),
              readdlm(string(path_csv, "epr_trend_1900_2100_Q.csv"), ',')/100.0,
              readdlm(string(path_csv, "parent_fertility_Q.csv"), ','),
              )
    return data
end
