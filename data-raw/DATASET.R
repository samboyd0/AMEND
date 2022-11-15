## code to prepare `DATASET` dataset goes here

# load packages
library(GEOquery); library(limma); library(ECEA); library(oligo)
library(dplyr); library(mgu74av2.db); library(biomaRt); library(igraph); library(data.table)

#================
# Glut4 Data Set ----
#================
gse1 <- GEOquery::getGEO("GSE35378", GSEMatrix = TRUE)
data <- Biobase::exprs(gse1[[1]])

# Background-subtracting and normalizing with oligo package
re <- oligo::basicRMA(data, row.names(data))

## checking for genes with zero sample variance
s <- matrix()
r <- matrix()
for(i in 1:nrow(re)){
  t <- sd(re[i,])
  if(t == 0){
    s <- c(s, t)
    r <- c(r, i)
  }
}
var_check <- data.frame(var = s[-1], row = r[-1], name = row.names(re)[r[-1]])

# Design matrix
groups = c(rep("OX", 3), rep("OX_ctrl", 3), rep("KO", 3), rep("KO_ctrl", 3))
X = stats::model.matrix(~0+groups)
colnames(X) = sort(unique(groups))

# Fitting a linear model
fit <- limma::lmFit(re[-var_check$row,], design = X)

# Getting contrast coefficents for KO vs. Ctrl 1
cnt_KO = limma::makeContrasts(contrasts = "KO-KO_ctrl", levels = colnames(X))
fit_KO = limma::contrasts.fit(fit, contrasts = cnt_KO)
fit_KO <- limma::eBayes(fit_KO)
tt_KO <- limma::topTable(fit_KO, number = 100000, adjust.method = "BH")
esize_KO <- tt_KO$logFC
pval_KO <- tt_KO$P.Value
names(esize_KO) <- row.names(tt_KO)
names(pval_KO) <- row.names(tt_KO)

# Getting contrast coefficents for OX vs. Ctrl 2
cnt_OX = limma::makeContrasts(contrasts = "OX-OX_ctrl", levels = colnames(X))
fit_OX = limma::contrasts.fit(fit, contrasts = cnt_OX)
fit_OX <- limma::eBayes(fit_OX)
tt_OX <- limma::topTable(fit_OX, number = 100000, adjust.method = "BH")
esize_OX <- tt_OX$logFC
pval_OX <- tt_OX$P.Value
names(esize_OX) <- row.names(tt_OX)
names(pval_OX) <- row.names(tt_OX)

# Getting ECI values for each gene
eci <- ECEA::getECI(smd1 = esize_OX, smd2 = esize_KO, p1 = pval_OX, p2 = pval_KO)
eci <- data.frame(ECI = eci, row.names = names(eci))

# Mapping probe identifier to gene symbol & aggregating ECI values and logFC values for KO vs. Ctrl
datFC = data.frame(logFC = tt_KO$logFC * (1 - tt_KO$P.Value), row.names = row.names(tt_KO))
dat <- eci
sym_map <- mgu74av2.db::mgu74av2SYMBOL

# Map the probes
genes <- unlist(as.list(sym_map[rownames(dat)]))

# Remove probes that didn't map
genes <- na.omit(genes)
common <- intersect(names(genes), rownames(dat))
genes <- genes[common]
dat <- dat[common,]
datFC = datFC[common,]

# Aggregate probes... take median ECI if multiple probes map to same gene
dat <- aggregate(dat, by = list(genes), FUN = median)
names(dat) <- c("Gene", "ECI")
row.names(dat) <- dat$Gene
dat <- dat %>% dplyr::select(-Gene)

# Aggregate probes... take median logFC if multiple probes map to same gene
datFC <- aggregate(datFC, by = list(genes), FUN = median)
names(datFC) <- c("Gene", "logFC")
row.names(datFC) <- datFC$Gene
datFC <- datFC %>% dplyr::select(-Gene)


#=========================#
# Corresponding PPI Network ----
#=========================#
# Getting mapping between Mouse Gene Symbols <--> Ensembl Peptide IDs
gene_symbols <- row.names(dat)
ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mapping <- biomaRt::getBM(attributes=c("mgi_symbol", "ensembl_peptide_id"),
                          filters = "mgi_symbol",
                          values = gene_symbols,
                          mart = ensembl)
mapping <- mapping %>% dplyr::filter(ensembl_peptide_id != "")

# Function for removing the duplicate edges from String PPI network (PPIN)
remove_duplicate_edges = function(x){
  g = igraph::graph_from_edgelist(as.matrix(x[,1:2]), directed = F)
  igraph::E(g)$weight = x[,3] / 1000
  adjM = as.matrix(igraph::as_adjacency_matrix(g, attr = "weight", sparse = T)) / 2
  g = igraph::graph_from_adjacency_matrix(adjM, mode = "undirected", weighted = TRUE)
  el = igraph::as_edgelist(g)
  return(cbind(el, igraph::E(g)$weight))
}

# Loading in mus musculus PPIN from String v11.0... can be downloaded from STRING website, but make sure you're on the right version
file_path = "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/10090.protein.links.v11.0.txt"
ppi = data.table::fread(file = file_path,
                        header = TRUE, sep = " ") %>%
  dplyr::filter(combined_score >= 800) %>%
  dplyr::mutate(protein1 = do.call(c, lapply(strsplit(protein1, "\\."), function(x) x[2])),
         protein2 = do.call(c, lapply(strsplit(protein2, "\\."), function(x) x[2])))  %>%
  as.data.frame() %>%
  remove_duplicate_edges()

# Keeping only genes that map to a protein in PPIN
mapping = mapping %>%
  dplyr::mutate(ppi = ensembl_peptide_id %in% ppi[,1] | ensembl_peptide_id %in% ppi[,2]) %>%
  dplyr::filter(ppi) %>%
  dplyr::group_by(mgi_symbol) %>%
  dplyr::mutate(N = sum(ppi)) %>%
  dplyr::filter(row_number() == 1) %>%
  as.data.frame()

# Keeping only proteins that map to a gene in the experiment
ppi = ppi[ppi[,1] %in% mapping$ensembl_peptide_id & ppi[,2] %in% mapping$ensembl_peptide_id,]

# Getting largest connected component of PPIN as an igraph object
g = igraph::graph_from_edgelist(as.matrix(ppi[,1:2]), directed = FALSE)
igraph::E(g)$weight = as.numeric(ppi[,3])
igraph::V(g)$symbol = mapping$mgi_symbol[match(V(g)$name, mapping$ensembl_peptide_id)]
l = igraph::decompose(g)
glut4_graph = l[[which.max(do.call(c, lapply(l, igraph::vcount)))]]

# Setting vertex attribute ECI and logFC (logFC for KO vs. Ctrl)
igraph::V(glut4_graph)$ECI = dat$ECI[match(igraph::V(glut4_graph)$symbol, row.names(dat))]
igraph::V(glut4_graph)$logFC = datFC$logFC[match(igraph::V(glut4_graph)$symbol, row.names(datFC))]

# Running the Louvain clusering algorithm to get a graph that is more manageable and quicker to run for examples
glut4_clusters = igraph::cluster_louvain(glut4_graph, weights = NULL)

# Taking largest cluster
glut4_graph = igraph::induced_subgraph(glut4_graph, glut4_clusters$membership == as.numeric(names(sort(table(glut4_clusters$membership), decreasing = T))[1]))

# Getting adjacency matrix from graph
glut4_adjM = igraph::as_adjacency_matrix(glut4_graph, attr = "weight", sparse = T)

# Creating vector of gene-wise experimental scores
eci_scores = igraph::V(glut4_graph)$ECI
logFC_KO = igraph::V(glut4_graph)$logFC

# write to package
usethis::use_data(glut4_graph, glut4_adjM, eci_scores, logFC_KO, overwrite = TRUE)

