#### This file contains the script for stepwise analysis of data, utilized for figures
```
# =========================================
# Xenium S. aureus niche mapping (CELL-BASED, r=20, k=5)
# + consistent "low → high" ordering
# + efferocytosis module score
# =========================================
# R 4.3+ ; Seurat v5 compatible

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(writexl)
})

set.seed(1)
sf::sf_use_s2(FALSE)  # planar ops; Xenium coords are in µm

# ---------------------------
# Paths & output folders
# ---------------------------
base_dir   <- "/path/xenium/MARCO_xenium_data"
plots_dir  <- file.path(base_dir, "results", "plots")
tables_dir <- file.path(base_dir, "results", "tables")
rds_dir    <- file.path(base_dir, "results", "rds")
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir,    recursive = TRUE, showWarnings = FALSE)

# Xenium output folder (edit if needed)
xenium_path <- "/path/xenium/20250207__180856__MAC-BFR-Roy-1_02072025/output-XETG00281__0033560__EDT-073_SV1__20250207__180935"

# ---------------------------
# Load Xenium into Seurat
# ---------------------------
x <- LoadXenium(xenium_path, fov = "fov", assay = "Xenium")
x <- subset(x, subset = nCount_Xenium > 0)

# ---------------------------
# Helpers
# ---------------------------
LOW_HIGH <- c("low","high")

save_and_print <- function(p, file, w=8, h=8, dpi=300) {
  ggsave(file.path(plots_dir, file), p, width=w, height=h, dpi=dpi)
  print(p)
}

detect_xy <- function(df) {
  x_opts <- c("x","x_location","global_x","imagecol","pxl_col_in_fullres","xc","X","x_pos","xpos","xcoord","x_coord")
  y_opts <- c("y","y_location","global_y","imagerow","pxl_row_in_fullres","yc","Y","y_pos","ypos","ycoord","y_coord")
  xcol <- intersect(x_opts, colnames(df))[1]
  ycol <- intersect(y_opts, colnames(df))[1]
  if (is.na(xcol) || is.na(ycol)) stop("No x/y-like columns found.")
  if (!"x" %in% colnames(df)) df[["x"]] <- df[[xcol]]
  if (!"y" %in% colnames(df)) df[["y"]] <- df[[ycol]]
  list(df=df, x="x", y="y")
}

plot_transcripts_points <- function(df, xcol, ycol, out_prefix) {
  df <- df %>% filter(is.finite(.data[[xcol]]), is.finite(.data[[ycol]]))
  stopifnot(nrow(df) > 0)
  p_pts <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(size = 0.35, alpha = 0.7) +
    coord_fixed() + theme_void() +
    ggtitle("S. aureus transcripts (points)")
  save_and_print(p_pts, paste0(out_prefix, "_points.png"), 8, 7)
  
  p_den <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    stat_density_2d(geom = "raster", aes(fill = after_stat(density)), contour = FALSE) +
    coord_fixed() + theme_void() +
    ggtitle("S. aureus transcripts (density)")
  save_and_print(p_den, paste0(out_prefix, "_density.png"), 8, 7)
}

# Collapse case-insensitive duplicate genes (e.g., Pcbp1 + PCBP1 -> PCBP1) and
# drop panel/bacterial IDs like "400-SA-16S". Do this on the macrophage object.
collapse_case_dups <- function(obj, assay="Xenium") {
  suppressPackageStartupMessages(library(Matrix))
  DefaultAssay(obj) <- assay
  M <- GetAssayData(obj, slot="counts", assay=assay)
  genes <- rownames(M)
  grp <- toupper(genes)
  M_sum <- rowsum(as.matrix(M), group = grp)
  rownames(M_sum) <- unique(grp)
  keep <- !grepl("^[0-9]+-", rownames(M_sum))
  M_sum <- M_sum[keep, , drop=FALSE]
  obj[[assay]] <- CreateAssayObject(counts = Matrix::Matrix(M_sum, sparse = TRUE))
  obj
}

# ---------------------------
# Read raw transcripts and filter S. aureus
# ---------------------------
suppressWarnings({
  files <- list.files(
    xenium_path,
    pattern = "(transcript|molecule).*\\.(parquet|csv(\\.gz)?)$",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
})
stopifnot(length(files) > 0)
parquet_file <- files[grepl("\\.parquet$", files, ignore.case = TRUE)][1]
csv_file     <- files[grepl("\\.csv(\\.gz)?$", files, ignore.case = TRUE)][1]

if (!is.na(parquet_file)) {
  suppressPackageStartupMessages(library(arrow))
  df <- read_parquet(parquet_file); message("Loaded transcripts from: ", parquet_file)
} else {
  suppressPackageStartupMessages(library(data.table))
  df <- data.table::fread(csv_file); message("Loaded transcripts from: ", csv_file)
}
df <- as.data.frame(df)

gene_col <- intersect(c("feature_name","gene","target","Gene","name"), colnames(df))[1]
stopifnot(!is.na(gene_col))
df[[gene_col]] <- gsub("_","-", df[[gene_col]])

# S. aureus probes (e.g., 400-SA-16S, 401-SA-agrB, …)
df_sa <- dplyr::filter(df, grepl("^[0-9]+-SA[-_]", .data[[gene_col]]))
stopifnot(nrow(df_sa) > 0)

xy_sa <- detect_xy(df_sa)
df_sa <- xy_sa$df

# Quick point & density plots
plot_transcripts_points(df_sa, "x", "y", "00B_transcripts_SA_fromDisk")

# ==========================
# CELL-BASED NICHES (no KDE/polygons)
# Rule: HIGH if (in-cell SA ≥1) OR (SA within r μm neighborhood ≥ k)
# ==========================
coords_df <- as.data.frame(GetTissueCoordinates(x))
xy_candidates <- intersect(
  c("x","y","imagecol","imagerow","pxl_col_in_fullres","pxl_row_in_fullres"),
  colnames(coords_df)
)
stopifnot(length(xy_candidates) >= 2)
xcol <- xy_candidates[1]; ycol <- xy_candidates[2]

if (!is.null(rownames(coords_df)) && length(intersect(rownames(coords_df), colnames(x))) > 0) {
  cells_df <- tibble::tibble(cell = rownames(coords_df),
                             x = coords_df[[xcol]],
                             y = coords_df[[ycol]])
} else if (nrow(coords_df) == ncol(x)) {
  cells_df <- tibble::tibble(cell = colnames(x),
                             x = coords_df[[xcol]],
                             y = coords_df[[ycol]])
} else {
  stop("Could not reconcile cell IDs between Seurat and coordinates.")
}

# Count SA per cell via segmentation
sa_by_cell <- df_sa %>%
  dplyr::filter(!is.na(cell_id) & cell_id != "") %>%
  dplyr::mutate(cell_id = as.character(cell_id)) %>%
  dplyr::count(cell_id, name = "sa_in_cell")

cells_df <- cells_df %>%
  dplyr::left_join(sa_by_cell, by = c("cell" = "cell_id")) %>%
  dplyr::mutate(sa_in_cell = dplyr::coalesce(sa_in_cell, 0L))

# Neighborhood SA within r μm
sa_pts_sf <- sf::st_as_sf(df_sa, coords = c("x","y"), crs = NA)
cells_sf  <- sf::st_as_sf(cells_df, coords = c("x","y"), crs = NA)

r_um <- 20
k_neigh <- 5
cells_buf  <- sf::st_buffer(sf::st_geometry(cells_sf), dist = r_um)
cells_df$sa_neigh_r <- lengths(sf::st_intersects(cells_buf, sf::st_geometry(sa_pts_sf)))

# Call niches
cells_df$niche <- ifelse(cells_df$sa_in_cell >= 1 | cells_df$sa_neigh_r >= k_neigh, "high", "low")
cells_df$niche <- factor(cells_df$niche, levels = LOW_HIGH)  # enforce order

# Write to Seurat meta
x$sa_in_cell  <- cells_df$sa_in_cell[match(colnames(x), cells_df$cell)]
x$sa_neigh_r  <- cells_df$sa_neigh_r[match(colnames(x), cells_df$cell)]
x$niche       <- cells_df$niche[match(colnames(x), cells_df$cell)]
x$niche       <- factor(x$niche, levels = LOW_HIGH)          # enforce order in metadata

message(sprintf("Niche rule: HIGH if SA_in_cell ≥1 OR SA_within_%dum ≥ %d", r_um, k_neigh))
print(table(x$niche, useNA="ifany"))

# Visual QC (legend & mapping fixed to low→high)
p_cells_overlay <- ggplot() +
  geom_point(data = cells_df, aes(x = x, y = y, color = niche),
             size = 0.35, alpha = 0.65) +
  scale_color_manual(breaks = LOW_HIGH, values = c(low = "#2ca02c", high = "#d62728")) +
  geom_point(data = df_sa, aes(x = x, y = y), color = "black", size = 0.15, alpha = 0.4) +
  coord_fixed() + theme_void() +
  ggtitle(sprintf("Cell-based niches (r=%d μm, in-cell ≥1 OR neigh ≥%d)", r_um, k_neigh))
save_and_print(p_cells_overlay, "40_cell_based_niches_r20k5.png", 8, 8, 300)

p_rule <- ggplot(cells_df, aes(sa_neigh_r, fill=niche)) +
  geom_histogram(position="identity", alpha=0.6, bins=40) +
  geom_vline(xintercept = k_neigh, linetype=2) +
  scale_fill_manual(breaks = LOW_HIGH, values=c(low="#2ca02c", high="#d62728")) +
  theme_bw() + xlab(sprintf("Staph within %d μm", r_um)) + ylab("Cells") +
  ggtitle("Niche calling rule (neighborhood component)")
save_and_print(p_rule, "40C_rule_hist_neigh.png", 6, 4, 300)

# Summary table (cells)
cell_summary <- cells_df %>%
  count(niche, name = "n_cells") %>%
  mutate(frac = n_cells / sum(n_cells), pct = 100*frac)
write_xlsx(list("cell_based_counts" = cell_summary),
           file.path(tables_dir, "40B_cell_based_counts.xlsx"))

# ==========================
# Macrophage subset (CD68+ & LYZ+) with low→high ordering
# ==========================
DefaultAssay(x) <- "Xenium"
expr_mat <- GetAssayData(x, slot="counts", assay="Xenium")

cd68_pos <- if ("CD68" %in% rownames(expr_mat)) expr_mat["CD68", ] > 0 else rep(FALSE, ncol(x))
lyz_pos  <- if ("LYZ"  %in% rownames(expr_mat)) expr_mat["LYZ",  ] > 0 else rep(FALSE, ncol(x))
mph_cells <- colnames(x)[cd68_pos & lyz_pos & !is.na(x$niche)]
mph <- subset(x, cells = mph_cells)

# Clean duplicates / panel IDs inside macrophage object
mph <- collapse_case_dups(mph, assay="Xenium")

# Standard workflow
DefaultAssay(mph) <- "Xenium"
mph <- NormalizeData(mph)
mph <- FindVariableFeatures(mph, nfeatures=2000)
mph <- ScaleData(mph)
mph <- RunPCA(mph, npcs=30)
mph <- FindNeighbors(mph, dims=1:20)
mph <- FindClusters(mph, resolution=0.2)
mph <- RunUMAP(mph, dims=1:20)

# Ensure niche is ordered low→high everywhere
mph$niche <- factor(mph$niche, levels = LOW_HIGH)

p_umap_split <- DimPlot(mph, reduction="umap", group.by="seurat_clusters",
                        split.by="niche", ncol=2) +
  ggtitle("Macrophages UMAP (split by niche, res=0.2)")
save_and_print(p_umap_split, "50_mph_umap_split_by_niche_r20k5.png", 10, 5)

# ==========================
# MARCO / MERTK readouts (with fixed order)
# ==========================
mat_mph <- GetAssayData(mph, slot="counts", assay="Xenium")
have <- rownames(mat_mph)
CD68_pos_m <- if ("CD68" %in% have) mat_mph["CD68",] > 0 else rep(FALSE, ncol(mph))

MARCO_pos <- if ("MARCO" %in% have) mat_mph["MARCO", ] > 0 else rep(FALSE, ncol(mph))
MERTK_pos <- if ("MERTK" %in% have) mat_mph["MERTK", ] > 0 else rep(FALSE, ncol(mph))

mph$MARCO_pos <- MARCO_pos & CD68_pos_m
mph$MERTK_pos <- MERTK_pos & CD68_pos_m

# Dot-free violins with thin box overlays; split panels ordered by factor(mph$niche)
plots <- list()
if ("MARCO" %in% have) {
  plots$MARCO_in_MERTK <- VlnPlot(
    mph, features="MARCO", group.by="MERTK_pos", split.by="niche",
    pt.size=0, raster=TRUE
  ) + geom_boxplot(width=0.1, outlier.shape=NA, alpha=0.35) +
    ggtitle("MARCO expression in MERTK+ vs MERTK− (by niche)")
}
if ("MERTK" %in% have) {
  plots$MERTK_in_MARCO <- VlnPlot(
    mph, features="MERTK", group.by="MARCO_pos", split.by="niche",
    pt.size=0, raster=TRUE
  ) + geom_boxplot(width=0.1, outlier.shape=NA, alpha=0.35) +
    ggtitle("MERTK expression in MARCO+ vs MARCO− (by niche)")
}
for (nm in names(plots)) save_and_print(plots[[nm]], paste0("54_", nm, ".png"), 9, 4)

# UMAP highlighting MARCO+ in red, split by niche (facet order = low→high)
mph$MARCO_pos_fac <- ifelse(mph$MARCO_pos, "MARCO+", "MARCO−")
p_umap_marco <- DimPlot(
  mph, reduction="umap", group.by="MARCO_pos_fac", split.by="niche", ncol=2,
  cols = c("MARCO−" = "grey80", "MARCO+" = "#d62728"), pt.size = 0.5, shuffle = TRUE
) + ggtitle("MARCO+ cells highlighted (split by niche)")
save_and_print(p_umap_marco, "50B_mph_umap_MARCOpos_split_r20k5.png", 10, 5)

# % MARCO+ and % MERTK+ per niche (+ prop.test), with fixed bar/legend order
frac_by_marker <- function(marker_logical, marker_name) {
  tab <- table(mph$niche, marker_logical)
  if (!"TRUE" %in% colnames(tab)) tab <- cbind(tab, `TRUE`=0)
  tab <- tab[,c("FALSE","TRUE")]
  totals    <- rowSums(tab)
  successes <- tab[,"TRUE"]
  pct       <- 100 * successes / totals
  pt <- prop.test(x = c(successes["low"], successes["high"]),
                  n = c(totals["low"], totals["high"]), correct = FALSE)
  tibble::tibble(
    marker = marker_name,
    niche  = factor(c("low","high"), levels = LOW_HIGH),
    n_cells = c(totals["low"], totals["high"]),
    n_pos   = c(successes["low"], successes["high"]),
    pct_pos = c(pct["low"], pct["high"]),
    p_value = pt$p.value
  )
}
marco_tbl <- frac_by_marker(mph$MARCO_pos, "MARCO")
mertk_tbl <- frac_by_marker(mph$MERTK_pos, "MERTK")
frac_stats <- dplyr::bind_rows(marco_tbl, mertk_tbl)
write_xlsx(list("MARCO_MERTK_percent_by_niche" = frac_stats),
           file.path(tables_dir, "53C_fraction_percentages.xlsx"))

p_pct_marco <- marco_tbl %>%
  ggplot(aes(x = niche, y = pct_pos, fill = niche)) +
  geom_col(width=0.6) + ylim(0,25) +
  scale_x_discrete(limits = LOW_HIGH) +
  scale_fill_manual(breaks = LOW_HIGH, values=c(low="#1f77b4", high="#d62728")) +
  theme_bw() + labs(title="MARCO+ macrophages (%)", y="% of macrophages", x=NULL)
save_and_print(p_pct_marco, "53A_MARCO_percent.png", 5, 4)

p_pct_mertk <- mertk_tbl %>%
  ggplot(aes(x = niche, y = pct_pos, fill = niche)) +
  geom_col(width=0.6) + ylim(0,25) +
  scale_x_discrete(limits = LOW_HIGH) +
  scale_fill_manual(breaks = LOW_HIGH, values=c(low="#1f77b4", high="#d62728")) +
  theme_bw() + labs(title="MERTK+ macrophages (%)", y="% of macrophages", x=NULL)
save_and_print(p_pct_mertk, "53B_MERTK_percent.png", 5, 4)


# ===== Nearest Staph distance (μm) — MARCO+ macrophages =====
# Make sure sf planar ops are used
sf::sf_use_s2(FALSE)

# 1) Build macrophage coordinates in **the same order as colnames(mph)**
mph_coord_df <- cells_df %>%
  dplyr::filter(cell %in% colnames(mph)) %>%
  dplyr::slice(match(colnames(mph), cell))

stopifnot(nrow(mph_coord_df) == ncol(mph))
mph_sf <- sf::st_as_sf(mph_coord_df, coords = c("x","y"), crs = NA)

# 2) SA points (reuse if already present)
if (!exists("sa_pts_sf")) {
  sa_pts_sf <- sf::st_as_sf(df_sa, coords = c("x","y"), crs = NA)
}

# 3) For each macrophage, find nearest SA dot and distance (μm)
nearest_idx <- sf::st_nearest_feature(mph_sf, sa_pts_sf)
mph$nearest_sa_um <- as.numeric(
  sf::st_distance(mph_sf, sa_pts_sf[nearest_idx,], by_element = TRUE)
)

# 4) Violin for MARCO+ macrophages, ordered low → high with fixed colors
LOW_HIGH <- c("low","high")
mph$niche <- factor(mph$niche, levels = LOW_HIGH)

d_marco <- FetchData(mph, vars = c("niche","nearest_sa_um","MARCO_pos")) %>%
  dplyr::as_tibble() %>% dplyr::filter(MARCO_pos)

pval <- wilcox.test(nearest_sa_um ~ niche, data = d_marco)$p.value
stars <- ifelse(pval < 1e-3,"***", ifelse(pval < 1e-2,"**", ifelse(pval < 0.05,"*","ns")))
ymax  <- max(d_marco$nearest_sa_um, na.rm = TRUE)

p_near <- ggplot(d_marco, aes(x = niche, y = nearest_sa_um, fill = niche)) +
  geom_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
  scale_x_discrete(limits = LOW_HIGH) +
  scale_fill_manual(breaks = LOW_HIGH, values = c(low = "#2ca02c", high = "#d62728")) +
  theme_bw() +
  labs(title = "Mph–MARCO+  nearest staph distance (μm)",
       x = NULL, y = "nearest staph distance (μm)") +
  annotate("segment", x = 1, xend = 2, y = ymax*1.03, yend = ymax*1.03) +
  annotate("text", x = 1.5, y = ymax*1.06, label = stars, size = 5)

save_and_print(p_near, "57_nearest_staph_distance_MARCOpos_violin.png", 4.5, 4)

# 5) Save summary + Wilcoxon
sum_tbl <- d_marco %>%
  dplyr::group_by(niche) %>%
  dplyr::summarise(n = dplyr::n(),
                   mean = mean(nearest_sa_um, na.rm = TRUE),
                   median = median(nearest_sa_um, na.rm = TRUE),
                   .groups = "drop")
w_tbl <- tibble::tibble(test = "Wilcoxon (two-sided)", p_value = pval)

writexl::write_xlsx(
  list(summary = sum_tbl, wilcoxon = w_tbl),
  file.path(tables_dir, "57_nearest_staph_distance_MARCOpos.xlsx")
)


# --- put this once after you create `mph` ---
LOW_HIGH <- c("low","high")
mph$niche <- factor(mph$niche, levels = LOW_HIGH)

# ==========================
# Efferocytosis module score (robust AddModuleScore)
# ==========================
# Curated list (uppercase to match Xenium feature names)
eff_genes <- toupper(c(
  "MERTK","AXL","TYRO3","GAS6","PROS1","MRC1","LGALS3","TIMD4","HAVCR2",
  "STAB1","STAB2","ITGAV","ITGB5","ITGB3","CD36","MFGE8","ANXA1","ADAM17","S1PR1"
))
have_features <- rownames(mph[["Xenium"]])
eff_avail    <- intersect(eff_genes, have_features)
eff_missing  <- setdiff(eff_genes, have_features)
message("Efferocytosis genes used: ", paste(eff_avail, collapse=", "))
if (length(eff_missing)) message("Missing (not in panel): ", paste(eff_missing, collapse=", "))

if (length(eff_avail) > 0) {
  # Use a realistic control pool and small ctrl to avoid sampling errors
  var_pool <- VariableFeatures(mph)
  if (length(var_pool) < 100) var_pool <- have_features
  ctrl_per_gene <- 10  # small and safe for targeted panels
  nbin_use      <- 12  # fewer bins → more genes per bin
  
  mph <- AddModuleScore(
    object  = mph,
    features= list(eff_avail),
    name    = "EfferocytosisScore",
    pool    = var_pool,
    ctrl    = ctrl_per_gene,
    nbin    = nbin_use,
    assay   = "Xenium"
  )
  score_col <- "EfferocytosisScore1"
  
  # Violin (no dots) + slim box; force low→high order + consistent colors
  p_eff_vln <- VlnPlot(
    mph, features = score_col, group.by = "niche", pt.size = 0
  ) +
    scale_x_discrete(limits = LOW_HIGH) +
    scale_fill_manual(breaks = LOW_HIGH, values = c(low="#2ca02c", high="#d62728")) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.35) +
    labs(title = "Efferocytosis module score by niche", x = NULL, y = "Module score")
  save_and_print(p_eff_vln, "56_efferocytosis_score_violin.png", 6, 4)
  
  # UMAP split by niche (facet order follows factor levels)
  p_eff_umap <- FeaturePlot(mph, features = score_col, split.by = "niche") +
    ggtitle("Efferocytosis module score (UMAP, split by niche)")
  save_and_print(p_eff_umap, "56B_efferocytosis_score_umap_split.png", 10, 5)
  
  # Summary + Wilcoxon (low vs high)
  eff_df <- FetchData(mph, vars = c("niche", score_col)) |> dplyr::as_tibble()
  eff_df$niche <- factor(eff_df$niche, levels = LOW_HIGH)
  eff_sum <- eff_df |>
    dplyr::group_by(niche) |>
    dplyr::summarise(
      n      = dplyr::n(),
      mean   = mean(.data[[score_col]], na.rm = TRUE),
      median = median(.data[[score_col]], na.rm = TRUE),
      .groups="drop"
    )
  w <- wilcox.test(eff_df[[score_col]] ~ eff_df$niche)
  eff_stats <- tibble::tibble(test = "Wilcoxon (two-sided)", p_value = w$p.value)
  
  writexl::write_xlsx(
    list("efferocytosis_summary" = eff_sum,
         "efferocytosis_wilcoxon" = eff_stats),
    file.path(tables_dir, "56_efferocytosis_module_score.xlsx")
  )
} else {
  message("No efferocytosis genes available on this panel — skipping module score.")
}

# ---- constants (keep colors/order consistent) ----
if (!exists("LOW_HIGH", inherits = FALSE))   LOW_HIGH   <- c("low","high")
if (!exists("SCALE_NICHE", inherits = FALSE)) SCALE_NICHE <- c(low="#2ca02c", high="#d62728")

# ---- ensure module score exists (compute once if missing) ----
score_col <- "EfferocytosisScore1"
if (!(score_col %in% colnames(mph@meta.data))) {
  eff_genes <- toupper(c(
    "MERTK","AXL","TYRO3","GAS6","PROS1","MRC1","LGALS3","TIMD4","HAVCR2",
    "STAB1","STAB2","ITGAV","ITGB5","ITGB3","CD36","MFGE8","ANXA1","ADAM17","S1PR1"
  ))
  have_features <- rownames(mph[["Xenium"]])
  eff_avail <- intersect(eff_genes, have_features)
  
  # small, safe control set to avoid the sampling error
  ctrl_pool <- setdiff(have_features, eff_avail)
  ctrl_n    <- max(1, min(25, length(ctrl_pool)))
  
  mph <- AddModuleScore(mph, features = list(eff_avail),
                        name = "EfferocytosisScore", ctrl = ctrl_n)
}

# ---- subset to MARCO+ macrophages ----
mph_marco <- subset(mph, subset = MARCO_pos)
mph_marco$niche <- factor(mph_marco$niche, levels = LOW_HIGH)

# ---- violin (no dots) + slim box, MARCO+ only ----
eff_df_m <- FetchData(mph_marco, vars = c("niche", score_col)) |>
  dplyr::as_tibble() |>
  dplyr::mutate(niche = factor(niche, levels = LOW_HIGH))

p_eff_vln_marco <- ggplot(eff_df_m, aes(x = niche, y = .data[[score_col]], fill = niche)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.35) +
  scale_x_discrete(limits = LOW_HIGH) +
  scale_fill_manual(values = SCALE_NICHE, breaks = LOW_HIGH) +
  labs(title = "Efferocytosis module score (MARCO+ only)", x = NULL, y = "Module score") +
  theme_bw() + guides(fill = "none")
save_and_print(p_eff_vln_marco, "56C_efferocytosis_score_violin_MARCOpos.png", 6, 4)

# ---- UMAP feature plot, MARCO+ only, split by niche ----
p_eff_umap_marco <- FeaturePlot(mph_marco, features = score_col, split.by = "niche") +
  ggtitle("Efferocytosis score (UMAP), MARCO+ only")
save_and_print(p_eff_umap_marco, "56D_efferocytosis_score_umap_split_MARCOpos.png", 10, 5)

# ---- summary + Wilcoxon (MARCO+ only) ----
eff_sum_m <- eff_df_m |>
  dplyr::group_by(niche) |>
  dplyr::summarize(
    n      = dplyr::n(),
    mean   = mean(.data[[score_col]], na.rm = TRUE),
    median = median(.data[[score_col]], na.rm = TRUE),
    .groups = "drop"
  )
w_m <- wilcox.test(.data[[score_col]] ~ niche, data = eff_df_m)

writexl::write_xlsx(
  list(
    efferocytosis_summary_MARCOpos  = eff_sum_m,
    efferocytosis_wilcoxon_MARCOpos = tibble::tibble(
      test = "Wilcoxon (two-sided, MARCO+ only)", p_value = w_m$p.value
    )
  ),
  file.path(tables_dir, "56C_efferocytosis_module_score_MARCOpos.xlsx")
)
message("MARCO+ efferocytosis score Wilcoxon p = ", signif(w_m$p.value, 3))


# ==========================
# DEG hi vs lo within macrophages + optional FGSEA
# ==========================
Idents(mph) <- factor(mph$niche, levels = LOW_HIGH)
deg <- FindMarkers(mph, ident.1="high", ident.2="low", logfc.threshold=0, min.pct=0.05)
deg$gene <- rownames(deg)
deg <- deg[!grepl("^[0-9]+-", rownames(deg)), , drop = FALSE]
deg <- deg %>% arrange(desc(avg_log2FC))
write_xlsx(list("DEG_hi_vs_lo"=deg),
           file.path(tables_dir, "51_mph_DEG_hi_vs_lo.xlsx"))

top_up   <- deg %>% filter(avg_log2FC > 0) %>% slice_head(n=15) %>% pull(gene)
top_down <- deg %>% filter(avg_log2FC < 0) %>% slice_head(n=15) %>% pull(gene)
top_genes <- unique(c(top_up, top_down))
if (length(top_genes) > 0) {
  p_heat <- DoHeatmap(mph, features=top_genes, group.by="niche", raster=TRUE) +
    ggtitle("Macrophage DEG heatmap (hi vs lo)")
  save_and_print(p_heat, "52_mph_DEG_heatmap.png", 8, 6, 300)
}

# Optional FGSEA (Hallmark) — safely stringify leadingEdge for Excel
if (!requireNamespace("msigdbr", quietly=TRUE)) install.packages("msigdbr")
if (!requireNamespace("fgsea", quietly=TRUE))   install.packages("fgsea")
suppressPackageStartupMessages({ library(msigdbr); library(fgsea) })
ranked <- deg %>% transmute(gene, stat = avg_log2FC) %>% filter(!is.na(stat)) %>% arrange(desc(stat))
ranks <- setNames(ranked$stat, ranked$gene)
msig <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
path_list <- split(msig$gene_symbol, msig$gs_name)
fg <- fgseaMultilevel(pathways = path_list, stats = ranks, minSize = 10, maxSize = 500)
fg <- dplyr::arrange(as_tibble(fg), padj)
fg$leadingEdge <- vapply(fg$leadingEdge, function(v) paste(v, collapse=";"), "")
write_xlsx(list("fgsea_hallmark_hi_vs_lo"=fg),
           file.path(tables_dir, "55_fgsea_hallmark_hi_vs_lo.xlsx"))

p_fg <- fg %>%
  slice_head(n=20) %>%
  mutate(pathway = factor(pathway, levels=rev(pathway))) %>%
  ggplot(aes(x=pathway, y=NES, fill=NES>0)) +
  geom_col() + coord_flip() + theme_bw() +
  labs(title="FGSEA Hallmark (hi vs lo)", x=NULL, y="Normalized Enrichment Score") +
  scale_fill_manual(values=c("TRUE"="#d62728","FALSE"="#1f77b4")) +
  theme(legend.position="none")
save_and_print(p_fg, "55_fgsea_hallmark_bar.png", 7, 6, 300)

# ====================================================
# DEG / FGSEA: MARCO+ macrophages only (low → high)
# ====================================================

# 0) Subset to MARCO+ and lock plotting/order
mph_marco <- subset(mph, subset = MARCO_pos)
mph_marco$niche <- factor(mph_marco$niche, levels = LOW_HIGH)

# quick sanity check
print(table(mph_marco$niche, useNA = "ifany"))
if (any(table(mph_marco$niche) < 10)) {
  message("Warning: fewer than 10 MARCO+ cells in one niche; DEG may be underpowered.")
}

# 1) DEG (hi vs lo) within MARCO+ macrophages
Idents(mph_marco) <- mph_marco$niche
deg_marco <- FindMarkers(
  mph_marco, ident.1 = "high", ident.2 = "low",
  logfc.threshold = 0, min.pct = 0.05
)
deg_marco$gene <- rownames(deg_marco)

# remove 10x/bacterial probe IDs (e.g., "400-SA-16s")
deg_marco <- deg_marco[!grepl("^[0-9]+-", deg_marco$gene), , drop = FALSE]
deg_marco <- dplyr::arrange(deg_marco, dplyr::desc(avg_log2FC))

# save table
write_xlsx(
  list("DEG_hi_vs_lo_MARCOpos" = deg_marco),
  file.path(tables_dir, "51_mphMARCOpos_DEG_hi_vs_lo.xlsx")
)

# 2) Heatmap of top genes (balanced up/down if available)
top_up   <- deg_marco |> dplyr::filter(avg_log2FC > 0) |> dplyr::slice_head(n = 40) |> dplyr::pull(gene)
top_down <- deg_marco |> dplyr::filter(avg_log2FC < 0) |> dplyr::slice_head(n = 40) |> dplyr::pull(gene)
top_genes <- unique(c(top_up, top_down))

if (length(top_genes) > 0) {
  p_heat_marco <- DoHeatmap(mph_marco, features = top_genes, group.by = "niche", raster = TRUE) +
    ggtitle("Macrophage MARCO+ DEG heatmap (hi vs lo)")
  save_and_print(p_heat_marco, "52_mphMARCOpos_DEG_heatmap.png", 8, 6, 300)
}

# 3) Optional FGSEA (Hallmark) for MARCO+ only
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
if (!requireNamespace("fgsea",   quietly = TRUE)) install.packages("fgsea")
suppressPackageStartupMessages({ library(msigdbr); library(fgsea) })

ranked_m <- deg_marco |>
  dplyr::transmute(gene, stat = avg_log2FC) |>
  dplyr::filter(!is.na(stat)) |>
  dplyr::arrange(dplyr::desc(stat))
ranks_m <- setNames(ranked_m$stat, ranked_m$gene)

msig_h <- msigdbr(species = "Homo sapiens", category = "H") |>
  dplyr::select(gs_name, gene_symbol) |> as.data.frame()
path_list <- split(msig_h$gene_symbol, msig_h$gs_name)

fg_m <- fgseaMultilevel(pathways = path_list, stats = ranks_m, minSize = 10, maxSize = 500)
fg_m  <- dplyr::arrange(as_tibble(fg_m), padj)
fg_m$leadingEdge <- vapply(fg_m$leadingEdge, function(v) paste(v, collapse = ";"), "")

write_xlsx(
  list("fgsea_hallmark_hi_vs_lo_MARCOpos" = fg_m),
  file.path(tables_dir, "55_fgsea_hallmark_hi_vs_lo_MARCOpos.xlsx")
)

p_fg_m <- fg_m |>
  dplyr::slice_head(n = 20) |>
  dplyr::mutate(pathway = factor(pathway, levels = rev(pathway))) |>
  ggplot(aes(x = pathway, y = NES, fill = NES > 0)) +
  geom_col() + coord_flip() + theme_bw() +
  labs(title = "FGSEA Hallmark (MARCO+ hi vs lo)", x = NULL, y = "Normalized Enrichment Score") +
  scale_fill_manual(values = c("TRUE" = "#d62728", "FALSE" = "#1f77b4")) +
  theme(legend.position = "none")
save_and_print(p_fg_m, "55_fgsea_hallmark_bar_MARCOpos.png", 7, 6, 300)

# ==========================================================
# Focused macrophage-function heatmap (low → high)
# ==========================================================

if (!exists("LOW_HIGH")) LOW_HIGH <- c("low","high")

# Choose which object to use
obj <- if (exists("mph_marco")) mph_marco else mph
DefaultAssay(obj) <- "Xenium"
obj$niche <- factor(obj$niche, levels = LOW_HIGH)

# ---- curated macrophage modules (UPPERCASE to match collapsed symbols)
mac_sets <- list(
  Scavenger_Receptors = c("MARCO","MSR1","CD36","MRC1","STAB1","STAB2","CD163"),
  Fc_Receptors        = c("FCGR1A","FCGR2A","FCGR3A","CD64","CD32","CD16"),
  TLR_Signaling       = c("TLR1","TLR2","TLR4","TLR5","MYD88","TICAM1","IRAK1","NFKB1"),
  Antigen_Presentation= c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74","CIITA"),
  IFN_Response        = c("STAT1","IRF7","ISG15","IFITM3","HERC5","GBP1"),
  Cytokines_Chemokines= c("IL1B","TNF","CCL2","CCL3","CCL4","CXCL9","CXCL10"),
  Lipid_Metabolism    = c("APOE","ABCA1","SREBF2","PLIN2","LPL"),
  Efferocytosis       = c("MERTK","AXL","TYRO3","MFGE8","GAS6","PROS1","ITGAV","CD36","TIMD4")
)

# collapse to mapping (gene -> category), keep unique one-to-one
map_df <- do.call(rbind, lapply(names(mac_sets), function(cat) {
  data.frame(gene = toupper(mac_sets[[cat]]), category = cat, stringsAsFactors = FALSE)
})) |> distinct()

available <- rownames(obj[["Xenium"]])
map_df <- dplyr::filter(map_df, gene %in% available)

# Use DEG table from the object you chose
deg_tbl <- if (exists("deg_marco") && identical(colnames(GetAssayData(obj))[1], colnames(GetAssayData(mph_marco))[1])) {
  deg_marco
} else {
  deg
}

# Keep significant DEGs that are in our modules
sig <- deg_tbl |>
  dplyr::filter(!is.na(p_val_adj) & p_val_adj < 0.05) |>
  dplyr::mutate(gene = toupper(gene)) |>
  dplyr::inner_join(map_df, by = "gene")

# If nothing passes padj, relax to top by |log2FC|
if (nrow(sig) == 0) {
  sig <- deg_tbl |>
    dplyr::mutate(gene = toupper(gene)) |>
    dplyr::inner_join(map_df, by = "gene")
  message("No macrophage-module genes with padj<0.05; showing top by |log2FC|.")
}

# pick top N per category by |log2FC|
genes_per_cat <- 6  # <— tune
pick <- sig |>
  dplyr::mutate(rank = dplyr::dense_rank(-abs(avg_log2FC))) |>
  dplyr::group_by(category) |>
  dplyr::slice_max(order_by = abs(avg_log2FC), n = genes_per_cat, with_ties = FALSE) |>
  dplyr::ungroup()

# order rows: by category, then effect (up first), then logFC
pick <- pick |>
  dplyr::mutate(effect = ifelse(avg_log2FC > 0, "high>low", "low>high")) |>
  dplyr::arrange(category, desc(avg_log2FC))

sel_genes  <- unique(pick$gene)
row_cat    <- setNames(pick$category, pick$gene)[sel_genes]

if (length(sel_genes) == 0) {
  stop("No curated macrophage-function genes found among DEGs & panel.")
}

# -------------------------------
# Heatmap with row/column splits
# -------------------------------
if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
    !requi
    
# ==========================
# Save objects
# ==========================
saveRDS(x,   file.path(rds_dir, "xenium_with_cell_based_niches_r20k5.rds"))
saveRDS(mph, file.path(rds_dir, "macrophages_CD68_LYZ_umap_res0p2_cleaned.rds"))

message("All done. Plots -> ", plots_dir, " | Tables -> ", tables_dir, " | RDS -> ", rds_dir)
```


Figure 2: Spatial mapping of macrophage transcripts and S. aureus niche in human cutaneous wounds by Xenium in situ. Xenium analysis

Figure 3: Neutralizing MARCO significantly improved biofilm-mediated efferocytosis dysfunction in MDM exposed to BFhi ∆rexB CM.

Figure 4: MARCO overexpression in myeloid cells impaired dead cell clearance. scRNA analysis

Figure S2: MARCO+ macrophage subsets showed the presence of both pro-inflammatory M1 and pro-resolution M2 markers. scRNA

Figure S3: Elevated MARCO expression in human patients with non-healing diabetic wounds. scRNA(Bhasin -RS)

Figure S5: Lyz2-promoter–driven MARCO overexpression via TNT alters macrophage polarization and reshapes WE proteomic signatures in murine skin wounds. Akoya analysis
