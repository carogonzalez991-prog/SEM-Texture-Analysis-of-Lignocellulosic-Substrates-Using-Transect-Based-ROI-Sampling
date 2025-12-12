# ============================================================
# SEM Transect ROI Pipeline
# Initial vs Residual substrate comparison
# ============================================================

# ---- Packages ----
pkgs <- c("EBImage","ggplot2","dplyr","tidyr","readr","boot","patchwork")
if (!requireNamespace("EBImage", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("EBImage")
}
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
if(length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only=TRUE)

# ================= USER INPUT =================

# --- Image paths (EDIT BY USER) ---
path_ini <- "data/SEM/initial_sample.tif"
path_res <- "data/SEM/residual_sample.tif"

# --- Scale bar (µm) ---
bar_um_ini <- 100
bar_um_res <- 100

# --- ROI settings ---
N_ROIS   <- 5
ROI_W_UM <- 90
ROI_H_UM <- 90
OUT_PX   <- 384

# --- Output files ---
csv_rois    <- "outputs/sem_rois.csv"
csv_summary <- "outputs/sem_summary.csv"
csv_stats   <- "outputs/sem_stats.csv"
png_boxes   <- "figures/SEM_boxes.png"
png_rois    <- "figures/SEM_rois.png"
png_plots   <- "figures/SEM_metrics.png"

# ================= FUNCTIONS =================

im_read_gray <- function(path){
  img <- EBImage::readImage(path)
  if (EBImage::colorMode(img) != "Grayscale")
    img <- EBImage::channel(img, "gray")
  rng <- range(img@.Data, na.rm=TRUE)
  if (diff(rng)>0) img <- (img-rng[1])/diff(rng)
  img
}

calibrate_um_px <- function(path, bar_um){
  img <- im_read_gray(path)
  a <- as.array(img); h <- dim(a)[1]; w <- dim(a)[2]
  plot(NA, xlim=c(0,w), ylim=c(h,0), asp=1, axes=FALSE)
  rasterImage(as.raster(img), 0, h, w, 0)
  title("Click both ends of the scale bar")
  pts <- locator(2)
  bar_px <- sqrt(diff(pts$x)^2 + diff(pts$y)^2)
  list(img=img, um_px=bar_um/bar_px)
}

var_laplacian <- function(img){
  K <- matrix(c(0,1,0,1,-4,1,0,1,0),3,3,byrow=TRUE)
  var(as.numeric(EBImage::filter2(img,K)))
}

otsu_manual <- function(img, nbins=256){
  x <- as.numeric(img); x <- x[is.finite(x)]
  h <- hist(x, breaks=nbins, plot=FALSE)
  p <- h$counts/sum(h$counts)
  w <- cumsum(p); m <- cumsum(p*(1:nbins)); mT <- m[nbins]
  s <- (mT*w - m)^2/(w*(1-w)+1e-12)
  i <- which.max(s)
  mean(h$breaks[i:(i+1)])
}

pct_bright <- function(img){
  eq <- EBImage::equalize(img)
  thr <- otsu_manual(eq)
  mean(eq >= thr)*100
}

glcm_feats <- function(img, n=32){
  M <- as.array(img)
  M <- floor((M-min(M))/(max(M)-min(M)+1e-12)*(n-1))
  G <- matrix(0,n,n)
  for(i in 1:(nrow(M)-1))
    for(j in 1:(ncol(M)-1))
      G[M[i,j]+1,M[i,j+1]+1] <- G[M[i,j]+1,M[i,j+1]+1]+1
  P <- G/sum(G)
  i <- row(P); j <- col(P)
  c(
    contrast=sum((i-j)^2*P),
    homog=sum(P/(1+(i-j)^2)),
    energy=sum(P^2),
    entropy=-sum(P[P>0]*log(P[P>0]))
  )
}

extract_rois <- function(img, um_px){
  a <- as.array(img); h <- dim(a)[1]; w <- dim(a)[2]
  px_w <- round(ROI_W_UM/um_px); px_h <- round(ROI_H_UM/um_px)
  plot(NA, xlim=c(0,w), ylim=c(h,0), asp=1, axes=FALSE)
  rasterImage(as.raster(img), 0, h, w, 0)
  title("Click transect start and end")
  pts <- locator(2)
  t <- seq(1,N_ROIS)/(N_ROIS+1)
  lapply(t,function(tt){
    cx <- round(pts$x[1]+tt*(pts$x[2]-pts$x[1]))
    cy <- round(pts$y[1]+tt*(pts$y[2]-pts$y[1]))
    roi <- img[(cy-px_h/2):(cy+px_h/2),
               (cx-px_w/2):(cx+px_w/2)]
    EBImage::resize(roi,OUT_PX,OUT_PX)
  })
}

metrics_rois <- function(rois,label){
  do.call(rbind,lapply(seq_along(rois),function(i){
    g <- glcm_feats(rois[[i]])
    data.frame(
      roi=i, condition=label,
      lap_var=var_laplacian(rois[[i]]),
      pct_bright=pct_bright(rois[[i]]),
      contrast=g["contrast"],
      homogeneity=g["homog"],
      energy=g["energy"],
      entropy=g["entropy"]
    )
  }))
}

# ================= PIPELINE =================

cal_ini <- calibrate_um_px(path_ini, bar_um_ini)
cal_res <- calibrate_um_px(path_res, bar_um_res)

rois_ini <- extract_rois(cal_ini$img, cal_ini$um_px)
rois_res <- extract_rois(cal_res$img, cal_res$um_px)

df <- rbind(
  metrics_rois(rois_ini,"Initial"),
  metrics_rois(rois_res,"Residual")
)

write_csv(df, csv_rois)

summary_tbl <- df |>
  pivot_longer(-c(roi,condition)) |>
  group_by(condition,name) |>
  summarise(mean=mean(value), sd=sd(value), .groups="drop")

write_csv(summary_tbl, csv_summary)

cat("\nSEM pipeline completed successfully\n")
