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
# experiment 1: GLUT4 overexpressor vs control
s <- matrix()
r <- matrix()
for(i in 1:nrow(re)){
  t <- sd(re[i, 1:6])
  if(t == 0){
    s <- c(s, t)
    r <- c(r, i)
  }
}
df1 <- data.frame(var = s[-1], row = r[-1], name = row.names(re)[r[-1]])
# experiment 2: GLUT4 knockout vs control
s <- matrix()
r <- matrix()
for(i in 1:nrow(re)){
  t <- sd(re[i, 7:12])
  if(t == 0){
    s <- c(s, t)
    r <- c(r, i)
  }
}
df2 <- data.frame(var = s[-1], row = r[-1], name = row.names(re)[r[-1]])


## Differential Expression analysis using limma
design <- cbind(ctrl = 1, trt = c(0, 0, 0, 1, 1, 1))

# logFC for over-expressor experiment
fit1 <- limma::lmFit(re[-df1$row,1:6], design)
fit1 <- limma::eBayes(fit1)
tt1 <- limma::topTable(fit1, coef = 2, number = dim(re)[1], adjust = "BH")
esize1 <- tt1$logFC
pval1 <- tt1$P.Value
names(esize1) <- row.names(tt1)
names(pval1) <- row.names(tt1)

# logFC for knockout experiment
fit2 <- limma::lmFit(re[-df2$row,7:12], design)
fit2 <- limma::eBayes(fit2)
tt2 <- limma::topTable(fit2, coef = 2, number = dim(re)[1], adjust = "BH")
esize2 <- tt2$logFC
pval2 <- tt2$P.Value
names(esize2) <- row.names(tt2)
names(pval2) <- row.names(tt2)

# Getting ECI values for each gene
eci <- ECEA::getECI(smd1 = esize1, smd2 = esize2, p1 = pval1, p2 = pval2)
eci <- data.frame(ECI = eci, row.names = names(eci))

# mapping probe identifier to gene symbol & aggregating ECI values
dat <- eci
sym_map <- mgu74av2.db::mgu74av2SYMBOL

# Map the probes
genes <- unlist(as.list(sym_map[rownames(dat)]))

# Remove probes that didn't map
genes <- na.omit(genes)
common <- intersect(names(genes), rownames(dat))
genes <- genes[common]
dat <- dat[common,]

# Aggregate probes... take median ECI if multiple probes map to same gene
dat <- aggregate(dat, by = list(genes), FUN = median)
names(dat) <- c("Gene", "ECI")
r <- dat$Gene
dat <- dat %>% dplyr::select(-Gene)
row.names(dat) <- r


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

# Loading in mus musculus PPIN from String v11.0
ppi = data.table::fread(file = "~/10090.protein.links.v11.0.txt",
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

# Getting largest connected component of PPIN as an igraph object with vertex attribute "ECI"
g = igraph::graph_from_edgelist(as.matrix(ppi[,1:2]), directed = FALSE)
igraph::E(g)$weight = as.numeric(ppi[,3])
igraph::V(g)$symbol = mapping$mgi_symbol[match(V(g)$name, mapping$ensembl_peptide_id)]
l = igraph::decompose(g)
glut4_graph = l[[which.max(do.call(c, lapply(l, igraph::vcount)))]]
igraph::V(glut4_graph)$ECI = dat$ECI[match(igraph::V(glut4_graph)$symbol, row.names(dat))]

# Running the Louvain clusering algorithm to get a graph that is more manageable and quicker to run for examples
glut4_clusters = igraph::cluster_louvain(glut4_graph, weights = NULL)

# Taking largest cluster
glut4_graph = igraph::induced_subgraph(glut4_graph, glut4_clusters$membership == as.numeric(names(sort(table(glut4_clusters$membership), decreasing = T))[1]))

# Getting adjacency matrix from graph
glut4_adjM = igraph::as_adjacency_matrix(glut4_graph, attr = "weight", sparse = T)

# Creating vector of gene-wise experimental scores, in this case the ECI values
eci_scores = igraph::V(glut4_graph)$ECI

# write to package
usethis::use_data(glut4_graph, glut4_adjM, eci_scores, overwrite = TRUE)

