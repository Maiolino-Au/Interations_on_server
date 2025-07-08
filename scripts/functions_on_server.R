# LOAD LIBRARIES
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(presto))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(tictoc))

suppressPackageStartupMessages(library(enrichR))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))

suppressPackageStartupMessages(library(SingleR))

suppressPackageStartupMessages(library(GPTCelltype))
suppressPackageStartupMessages(library(openai))

total_time <- function(seconds) {
    d <- seconds %/% (86400)
    h <- (seconds %% 86400) %/% 3600
    m <- (seconds %% 3600) %/% 60
    s <- seconds %% 60
    
    cat(sprintf("Total Time: %d Days, %d Hours, %d Minutes and %f Seconds\n", d, h, m, s))
}

load.data <- function(
    file_name,
    data_path = path_to_data,
    output = F,
    reduced.output = T 
) {
    if (output | reduced.output) {
        print(paste("Loading data for time point:", file_name))
    }

    # Load the data
    sc_data <- Read10X(data.dir = paste(data_path, "expression_", file_name, sep = ""), gene.column = 1)

    # Create Seurat object
    sc_data <- CreateSeuratObject(counts = sc_data, min.cells = 3, min.features = 500, project = file_name, names.delim = "-", names.field = 2)

    # Normalize the data
    sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = output)

    # Find variable features
    sc_data <- FindVariableFeatures(sc_data, selection.method = "mvp", nfeatures = 2000, verbose = output)

    # Scale the data
    sc_data <- ScaleData(sc_data, verbose = output)

    return(sc_data)
}

PCA.cluster <- function(
    data = sc_data_scaled, 
    file_name = timepoints[time_point], 
    res = 1, 
    n_dim = 40, 
    save = F,
    output = F,
    reduced.output = T    
) {  
    if (output | reduced.output) {
        print(paste("Running PCA and clustering for time point:", file_name))
        print(paste("- Resolution:", res))
        print(paste("- Dimensions:", n_dim))
    }
        
    # PCA
    data <- RunPCA(data, npcs = n_dim, verbose = output)
    #print(ElbowPlot(object = data, ndims = 50))

    # Cluster the cells
    data <- FindNeighbors(data, dims = 1:n_dim, verbose = output)
    data <- FindClusters(data, resolution = res, verbose = output)
    
    #print(table(Idents(data)))

    # Save the parial
    if (save) {
        name_new_dir <- paste("Partial/", file_name, "/cluster", param, sep="")
        if (!dir.exists(name_new_dir)) {dir.create(name_new_dir)} 
    
        print(paste("Saving PCA for time point", file_name, "in", name_new_dir))
        save(data, file = paste(name_new_dir, "/PCA_res_", res, "_dim_", n_dim, "_", file_name, ".Robj", sep=""))
    }
    return(data)
}

UMAP.plot <- function(
    data = sc_data,
    file_name = timepoints[time_point],
    n_dim = 40,
    param = "",
    output = F,
    reduced.output = T,
    print_plot = F
) {
    if (output | reduced.output) {
        print("Making UMAP")
    }
    
    sc_data_UMAP <- RunUMAP(data, dims = 1:n_dim, verbose = output)
    
    # Visualization of clusters   
    plot <- DimPlot(sc_data_UMAP, reduction = "umap", label = TRUE, pt.size = 1) + 
        ggtitle(paste("UMAP of Clusters -",file_name, "-", gsub("_", " ", param)))

    if (print_plot) {print(plot)}
    
    return(plot)
}

# FIND ALL MARKERS
cluster.markers <- function(
    data, 
    file_name = timepoints[time_point],
    output = F,
    reduced.output = T
) {
    if (output | reduced.output) {
        print(paste("Finding all markers for time point:", file_name))
    }

    # Find all markers for every cluster compared to all remaining cells
    markers <- FindAllMarkers(data,
                              only.pos = TRUE,   # Considera solo i marker espressi positivamente
                              min.pct = 0.25,    # Percentuale minima di espressione nelle cellule del cluster
                              logfc.threshold = 0.25,  # Soglia minima di LogFC
                              verbose = output)
        
    return(markers)
}

annotation.enrichR <- function (
    top_genes,
    database = "Allen_Brain_Atlas_10x_scRNA_2021",
    output = F,
    reduced.output = T
) {
    gene_symbols <- top_genes$gene
    # entrez_ids <- suppressMessages(
    #     mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    # )
    # # Add Entrez IDs to your data.frame
    # top_genes$entrez <- entrez_ids
    # # Remove rows with NA Entrez IDs
    # top_genes <- top_genes[!is.na(top_genes$entrez), ]

    annotation_list <- list()
    
    for (cl in unique(top_genes$cluster)) {
        genes_cluster <- top_genes %>% filter(cluster == cl) %>% pull(gene)#entrez)
        # Optionally, use Entrez IDs instead:
        # genes_cluster <- top_genes %>% filter(cluster == cl) %>% pull(entrez)
        
        # Perform enrichment analysis
        a <- capture.output(
            enriched <- enrichr(genes_cluster, databases = "Allen_Brain_Atlas_10x_scRNA_2021")
        )

        annotation <- enriched$Allen_Brain_Atlas_10x_scRNA_2021 %>% as.data.table()
        annotation$cluster <- as.numeric(cl)

        annotation_list[[cl]] <- annotation[grepl("Human", annotation$Term, ignore.case = TRUE) & annotation$Adjusted.P.value < 5e-2]
    }


    return(annotation_list)
}

iterations <- function(
    timepoints = c("23days", "1month", "1.5month", "2month", "3month", "4month", "5month", "6month"),
    housekeeping_genes = c("ACTB", "DLG4"),
    genes_of_interest = c("SRCIN1", "KIAA1217", "CIT"),
    path_to_data = "/Data/",
    res_list = seq(0.1, 1, by = 0.1),
    dim_list = c(20, 30, 40, 50),
    top_n_genes = 50,
    check_for_saved_plot = F,
    output = F,
    reduced.output = T,
    print.plot = F,
    name_save = "",
    exec.enrichR = T
) {
    print(paste("Path to data:", path_to_data))
    print(paste("Timepoints:", timepoints))
    print(paste("Housekeeping genes:", housekeeping_genes))
    print(paste("Genes of interest:", genes_of_interest))
    print("")

    print(paste("Resolutions:", res_list))
    print(paste("Number of dimension:", dim_list))
    print("")
    
    # Total number of 
    cycles_per_timepoint <- length(res_list)*length(dim_list)
    print(paste("Number of cycles per timepoint:", cycles_per_timepoint))
    tot_cycles <- cycles_per_timepoint*length(timepoints)
    print(paste("Total number of cycles:", tot_cycles))
    print("")
    
    dir_iterations <- "/Results/Iterations"
    if (!dir.exists(dir_iterations)) {dir.create(dir_iterations)} 
    
    N_cycle <- 0
    
    print("Starting")
    
    timings <- data.frame()
    n_of_clusters <- data.frame()
    
    tic("Iterations")
    
    for (f_name in timepoints) {
        tic(paste("Loading", f_name))
    
        dir_timepoint <- paste(dir_iterations, "/", f_name, sep="")     
        if (!dir.exists(dir_timepoint)) {dir.create(dir_timepoint)}
        
        # Load data
        sc_data_scaled <- load.data(file_name = f_name, data_path = path_to_data, output = output, reduced.output = reduced.output)
    
        toc(quiet = T)
    
        for (clustering_resolution in res_list) {
            for (n_of_dimesnions in dim_list) {
                # Progression Message
                N_cycle <- N_cycle + 1
                print(paste(N_cycle, "of", tot_cycles)) 
                print(paste(f_name, "- Res:", clustering_resolution, "- Dim:", n_of_dimesnions))
    
                # param
                char_res <- as.character(clustering_resolution)
                if (grepl("\\.", char_res)) {char_res <- gsub("\\.", "_", char_res)}
                param <- paste0("_res_", char_res, "_dim_", n_of_dimesnions)
    
                # Create subfolder
                new_dir <- paste(dir_timepoint, "/", f_name, param, sep="")     
                if (!dir.exists(new_dir)) {dir.create(new_dir)}
    
                # start time
                tic(paste(f_name, "- Res:", clustering_resolution, "- Dim:", n_of_dimesnions))
    
                # Clusterize
                sc_data <- PCA.cluster(
                    data = sc_data_scaled, 
                    file_name = f_name, 
                    res = clustering_resolution, 
                    n_dim = n_of_dimesnions,
                    output = output, 
                    reduced.output = reduced.output
                )
    
                n_of_clusters[paste0("res_", clustering_resolution, "_dim_", n_of_dimesnions), f_name] <- length(table(Idents(sc_data)))
    
                plot_UMAP <- suppressMessages(suppressWarnings(UMAP.plot(
                    data = sc_data, 
                    file_name = f_name, 
                    n_dim = n_of_dimesnions, 
                    param = param,
                    output = output, 
                    reduced.output = reduced.output,
                    print_plot = print.plot
                )))

                # Save UMAP
                umap_file_name <- paste0(new_dir, "/UMAP", f_name, param, ".png")
                if (!(check_for_saved_plot & exists(umap_file_name))) {
                    print("Saving UMAP")
                    ggsave(
                        umap_file_name, 
                        plot=plot_UMAP, 
                        width = 1920*2, height = 1980*2, units = "px"
                    )
                }
    
                # Find Markers
                cluster_markers <- cluster.markers(
                    data = sc_data, 
                    file_name = f_name, 
                    output = output, 
                    reduced.output = reduced.output
                )
    
                write.csv(cluster_markers, file = paste0(
                    new_dir, "/",
                    "cluster_markers_",
                    f_name,
                    param,
                    ".csv"
                ))

                # de.genes
                de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest)
                if (output | reduced.output) {print(paste("de.genes nrow:", nrow(de_genes)))}

                write.csv(de_genes, file = paste0(
                    new_dir, "/",
                    "de_genes",
                    f_name,
                    param,
                    ".csv"
                ))

                # Top top_n_genes DE genes
                top_genes_by_cluster <- cluster_markers %>% group_by(cluster) %>% top_n(n = top_n_genes, wt = avg_log2FC) %>% as.data.frame()

                write.csv(top_genes_by_cluster, file = paste0(
                    new_dir, "/",
                    "top_genes_by_cluster_",
                    f_name,
                    param,
                    ".csv"
                ))
    
                # Annotation - enrichR
                if (exec.enrichR) {
                    enrichR_list <- annotation.enrichR(
                        top_genes = top_genes_by_cluster,
                        database = "Allen_Brain_Atlas_10x_scRNA_2021"
                    )
    
                    write.csv(bind_rows(enrichR_list), file = paste0(
                        new_dir, "/",
                        "enrichR_annotation_",
                        f_name,
                        param,
                        ".csv"
                    ), row.names = FALSE)
                }
    
                # Annotation - singelR
    
                # Annotation - CPTCellType
    
                # end time
                elapsed <- toc(log = TRUE, quiet = reduced.output)
                timings[f_name, paste0("res_", clustering_resolution, "_dim_", n_of_dimesnions)] <- elapsed$toc - elapsed$tic
            }
        }

        print("") # an empti line for better visualization 
    }
    
    timings$tot_time <- rowSums(timings)
    
    n_tocken_gpt <- n_of_clusters*837.5
    n_tocken_gpt$total <- rowSums(n_tocken_gpt)

    if (name_save != F) {
        if (name_save != "") name_save <- paste0("_", name_save)
        write.csv(timings, file = paste0(dir_iterations, "/tmings", name_save, ".csv"))
        write.csv(n_tocken_gpt, file = paste0(dir_iterations, "/number_of_tokens", name_save, ".csv"))
    }
    
    tot_elapsed <- toc(log = TRUE, quiet = F)
    print("END")
    return(list(Time = timings, N_of_Clusters = n_of_clusters))
}