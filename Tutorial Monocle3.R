# tutorial on Monocle 3
library(monocle3)
library(dplyr)

#load the RDS parallel function, very useful (cool-salmon repo)

# Load the data
expression_matrix <- readRDS("/Users/carloleonardi/Downloads/cao_l2_expression.rds")
cell_metadata <- readRDS("/Users/carloleonardi/Downloads/cao_l2_colData.rds")
gene_annotation <- readRDS("/Users/carloleonardi/Downloads/cao_l2_rowData.rds")

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)

plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds)

plot_cells(cds)

plot_cells(cds, color_cells_by="cao_cell_type")

# take away batch effects
cds = align_cds(cds, num_dim = 100, alignment_group = "plate")
cds = reduce_dimension(cds)
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)

cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds)

marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)

colData(cds)$assigned_cell_type <- as.character(partitions(cds))

colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$assigned_cell_type,
                                                "1"="Germline",
                                                "2"="Body wall muscle",
                                                "3"="Unclassified neurons",
                                                "4"="Vulval precursors",
                                                "5"="Failed QC",
                                                "6"="Seam cells",
                                                "7"="Pharyngeal epithelia",
                                                "8"="Coelomocytes",
                                                "9"="Am/PH sheath cells",
                                                "10"="Failed QC",
                                                "11"="Touch receptor neurons",
                                                "12"="Intestinal/rectal muscle",
                                                "13"="Pharyngeal neurons",
                                                "14"="NA",
                                                "15"="flp-1(+) interneurons",
                                                "16"="Canal associated neurons",
                                                "17"="Ciliated sensory neurons",
                                                "18"="Other interneurons",
                                                "19"="Pharyngeal gland",
                                                "20"="Failed QC",
                                                "21"="Ciliated sensory neurons",
                                                "22"="Oxygen sensory neurons",
                                                "23"="Ciliated sensory neurons",
                                                "24"="Ciliated sensory neurons",
                                                "25"="Ciliated sensory neurons",
                                                "26"="Ciliated sensory neurons",
                                                "27"="Oxygen sensory neurons",
                                                "28"="Ciliated sensory neurons",
                                                "29"="Unclassified neurons",
                                                "30"="Socket cells",
                                                "31"="Failed QC",
                                                "32"="Pharyngeal gland",
                                                "33"="Ciliated sensory neurons",
                                                "34"="Ciliated sensory neurons",
                                                "35"="Ciliated sensory neurons",
                                                "36"="Failed QC",
                                                "37"="Ciliated sensory neurons",
                                                "38"="Pharyngeal muscle")

plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type")

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition",  colorRampPalette(brewer.pal(8, "Set2"))(nb.cols))
                                 
                                 
library(RColorBrewer)
# Define the number of colors you want
  nb.cols <- 9999999
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
# Create a ggplot with 18 colors 
  # Use scale_fill_manual
  ggplot(df) 
  +     geom_col(aes(name, Sepal.Length, fill = factor(Sepal.Length))) +
  +     scale_fill_manual(values = mycolors) +
  +     theme_minimal() +
  +     theme(legend.position = "top")

