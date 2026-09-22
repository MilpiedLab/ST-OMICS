
## @knitr prepare_data

cat(sprintf("%s **OVERVIEW AND OBJECTIVES**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")

cat("\n\n")
cat("Mutually Exclusive Co-expression Rate (MECR) is a metric that can be used 
to assess the quality of cell segmentation in imaging spatial transcriptomics data.

MECR is computed by dividing the fraction of cells that express both markers by 
the fraction of cells that express at least one of the markers, for each pair of 
mutually exclusive genes (Satija et al., bioRxiv, 2024).

Computing this metric requires access to a list of mutually exclusive gene pairs 
for the dataset being analyzed. These correspond to genes that are not expected to 
be co-expressed within the same cell type (i.e., genes whose expression is typically 
restricted to distinct cell types).

The aim of this analysis step is therefore to retrieve a list of mutually 
exclusive gene pairs from a reference single-cell dataset. The cell types present 
in this reference are expected to also be represented in the spatial dataset. 
The gene pairs identified in this analysis will be used in subsequent steps to 
assess the quality of cell segmentation.")
cat("\n\n")

################################################################################
#'                1. LOAD AND PREPARE REFERENCE SINGLE-CELL DATASET           
#'                
#' Load the reference single-cell dataset and manually complete the cell type 
#' annotations in the corresponding Seurat object metadata.        
#'        
################################################################################

# Load reference single cell data as Seurat object
reference_sc_data = readRDS(file = REFERENCE_SC_SEURAT_INPUT)

# Set Normalized data as default assay for DEG analysis
DefaultAssay(reference_sc_data) = "RNA"
Idents(reference_sc_data) = "celltype"

# Attributing hierarchical level to cell types targeted by the gene panel
reference_sc_data@meta.data = reference_sc_data@meta.data %>% 
  mutate(Celltype_Level_4 = celltype,
         Celltype_Level_3 = ifelse(Celltype_Level_4 %in% c("Activated.MBC", "MYC.activated.B"),"Non_malignant_activated_B.cells",
                                   ifelse(Celltype_Level_4 %in% c("MBC", "FCRL4.MBC"),"Non_malignant_Mem_B.cells",
                                          ifelse(Celltype_Level_4 %in% c("Up_valve", "Down_valve"),"valveLEC.cell",
                                                 ifelse(Celltype_Level_4 %in% c("cLEC", "collectLEC"),"c.collect.LEC.cell",
                                                        ifelse(Celltype_Level_4 %in% c("SFRP2_SC", "TNF_SC"),"NE_TNF_SC.cell",
                                                               ifelse(Celltype_Level_4 %in% c("SFRP4_SC", "C7_SC"),"NE_C7_SC.cell",Celltype_Level_4)))))),
         Celltype_Level_2 = ifelse(Celltype_Level_3 %in% c("CD14moMac_0", "S100AmoMac_5",  "CD16moMac_10"), "moMac.cell",
                                   ifelse(Celltype_Level_3 %in% c("CD4.Central.Naive", "CD4.TREG", "CD4.TFH", "CD4.Effector", "CD4.CTL"), "CD4_T.cell",
                                          ifelse(Celltype_Level_3 %in% c("CD8.Central", "CD8.Effector", "CD8.Effector.IFNG", "CD8.Naive", "CD8.Effector.CXCL13", "CD8.Effector.ENTPD1"),"CD8_T.cell",
                                                 ifelse(Celltype_Level_3 %in% c("Non_malignant_activated_B.cells", "PC.PB", "Naive.B", "GC.B", "Non_malignant_Mem_B.cells"),"Non_malignant_B.cells",
                                                        ifelse(Celltype_Level_3 %in% c("Mem_like", "PC_like", "Myc_act_like", "Myc_NFKB_like", "GC_like"),"Malignant_B.cell",
                                                               ifelse(Celltype_Level_3 %in% c("hHEV", "aHEV", "VBEC", "CXCL10_HEV"),"VenousBEC.cell",
                                                                      ifelse(Celltype_Level_3 %in% c("tBEC", "C_aHEV", "cBEC"), "capillaryBEC.cell",
                                                                             ifelse(Celltype_Level_3 %in% c("ABEC", "caBEC", "aBEC"),"arterialBEC.cell",
                                                                                    ifelse(Celltype_Level_3 %in% c("valveLEC.cell", "msLEC","c.collect.LEC.cell", "bLEC"),"otherLEC.cell",
                                                                                           ifelse(Celltype_Level_3 %in% c("NE_TNF_SC.cell", "NE_C7_SC.cell", "advSC"),"NE_SC.cell",
                                                                                                  ifelse(Celltype_Level_3 %in% c("ATF3lo_SMC", "ATF3hi_SMC","TRC", "Pericyte"),"NE_SMC.cell",
                                                                                                         ifelse(Celltype_Level_3 %in% c("FDC", "MRC"),"NE_FSC.cell",Celltype_Level_3)))))))))))),
         Celltype_Level_1 = ifelse(Celltype_Level_2 %in% c("moMac.cell", "SIGLEC6cDC_4", "CLEC9AcDC1_8","CD1EcDC2_3","CREMcDC2_2","M1_1","MacC1Q_7", "MNDAcDC2_6"),"Myeloid.cell",
                                   ifelse(Celltype_Level_2 == "pDC","pDC", 
                                          ifelse(Celltype_Level_2 == "LAMP3cDC_9","LAMP3cDC_9", 
                                                 ifelse(Celltype_Level_2 %in% c("CD4_T.cell","ISG", "CD8_T.cell", "DN.CD81"),"T.cells",
                                                        ifelse(Celltype_Level_2 %in% c("NK","ILC"),"NK_ILC.cells",
                                                               ifelse(Celltype_Level_2 %in% c("Non_malignant_B.cells","Malignant_B.cell"),"B.cell",
                                                                      ifelse(Celltype_Level_2 %in% c("VenousBEC.cell","capillaryBEC.cell","arterialBEC.cell"),"BEC.cell",
                                                                             ifelse(Celltype_Level_2 %in% c("fLEC_pfsLEC", "otherLEC.cell"),"LEC.cell",
                                                                                    ifelse(Celltype_Level_2 %in% c("NE_SC.cell", "NE_SMC.cell","NE_FSC.cell","AGT_SC"), "NESC.cell",Celltype_Level_2))))))))))

################################################################################
#'                        2. LOAD XENIUM SPATIAL DATASET                      
#'                        
#' Load the Xenium dataset to be analyzed in this experiment. From the Xenium 
#' Seurat object, retrieve the genes of interest (i.e., the gene panel) used to 
#' generate the data, in order to restrict the search for mutually exclusive 
#' gene pairs to this specific set of genes. Constraining the exclusive gene 
#' pairs search to the Xenium gene panel reduces computation time.    
#'                    
################################################################################

# Load Xenium data as Seurat object
Xenium_data = LoadXenium_new(XENIUM_DATA_PATH[["TMA002"]])

# Retrieve gene panel
gene_panel_480 = rownames(Xenium_data)
gene_panel_df = data.frame("gene" = gene_panel_480)

# Clean previous data
rm("Xenium_data")
gc()

################################################################################
#'                        3. CELL TYPE TREE VIZUALISATION                         
#'                           
#' Build a hierarchical tree of cell types based on the annotations in the 
#' Seurat object metadata. This hierarchical tree will be displayed in the HTML 
#' report corresponding to this analysis step. 
#'                          
################################################################################

# Extracting hierarchy df
hierarchy_type = reference_sc_data@meta.data %>%
  dplyr::select(Celltype_Level_1, Celltype_Level_2, Celltype_Level_3, Celltype_Level_4) %>%
  distinct() %>%
  as_tibble() %>%
  dplyr::arrange(Celltype_Level_1,Celltype_Level_2,Celltype_Level_3, Celltype_Level_4)

# Add level 0
hierarchy_type$Celltype_Level_0 = "cell_type"

# Create a copy of hierarchy type
data = hierarchy_type

# Create a column with hierarchical path
data$pathString = paste("Celltypes", data$Celltype_Level_1, data$Celltype_Level_2, 
                        data$Celltype_Level_3, data$Celltype_Level_4, sep = "/")

# Convert to data.tree structure
cell_tree = as.Node(data)

cat("\n\n")
cat(sprintf("%s **Hierarchy overview**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Hierarchy of cell types found in the reference single cell data set:")
cat("\n\n")

# Print the cell type tree
cat("\n\n")
collapsibleTree(cell_tree, root = "Celltypes")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s **Gene panel**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Gene panel used to generate the Xenium data:")
cat("\n\n")

cat("\n\n")
print(htmltools::tagList(DT::datatable(gene_panel_df)))
cat("\n\n")
