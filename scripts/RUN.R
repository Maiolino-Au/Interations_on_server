source(file = "/scripts/functions_on_server.R")

# timepoint can be ("23days", "1month", "1.5month", "2month", "3month", "4month", "5month", "6month")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Please provide a timepoint as the first argument (e.g. '23days')")
}
timepoint <- args[1]

tmp <- capture.output(
    iterations_list <- iterations(
        timepoints = timepoint,
        housekeeping_genes = c("ACTB", "DLG4"),
        genes_of_interest = c("SRCIN1", "KIAA1217", "CIT"),
        path_to_data = "/sharedFolder/Data/",
        res_list = seq(0.4, 1, by = 0.1),
        dim_list = c(20, 30, 40),
        top_n_genes = 50,
        check_for_saved_plot = T,
        output = F, 
        reduced.output = T, 
        print.plot = F,
        name_save = timepoint,
        exec.enrichR = T
    ),
    file = paste0("/sharedFolder/output_", timepoint, ".txt")
)