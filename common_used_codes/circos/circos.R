library(circlize)
library(graphics)

setwd('/Users/stead/Desktop/circos')
load("chromatin_transition.RData")

all_states = rownames(mat) 
n_states = nrow(mat)

rownames(mat) = paste0("R_", seq_len(n_states)) 
colnames(mat) = paste0("C_", seq_len(n_states))


state_col = c("TssA" = "#E41A1C",    "TssAFlnk" = "#E41A1C",
              "TxFlnk" = "#E41A1C",  "Tx" = "#E41A1C",
              "TxWk" = "#E41A1C",    "EnhG" = "#E41A1C",
              "Enh" = "#E41A1C",     "ZNF/Rpts" = "#E41A1C",
              "Het" = "#377EB8",     "TssBiv" = "#377EB8",
              "BivFlnk" = "#377EB8", "EnhBiv" = "#377EB8",
              "ReprPC" = "#377EB8",  "ReprPCWk" = "#377EB8",
              "Quies" = "black")

# one for rows and one for columns
state_col2 = c(state_col, state_col)
names(state_col2) = c(rownames(mat), colnames(mat))

colmat = rep(state_col2[rownames(mat)], n_states)
colmat = rgb(t(col2rgb(colmat)), maxColorValue = 255)

qati = quantile(mat, 0.7)
colmat[mat > qati] = paste0(colmat[mat > qati], "A0")
colmat[mat <= qati] = paste0(colmat[mat <= qati], "20")
dim(colmat) = dim(mat)

circos.par(start.degree = -5, gap.after = c(rep(1, n_states-1), 10, rep(1, n_states-1), 10), 
           cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE) ###split 

cdm_res = chordDiagram(mat, col = colmat, grid.col = state_col2, directional = TRUE, 
                       annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))

head(cdm_res)

###Create plotting regions for a whole track
circos.track(track.index = 2, panel.fun = function(x, y) {
  if(abs(CELL_META$cell.start.degree - CELL_META$cell.end.degree) > 3) {
    sn = CELL_META$sector.index
    i_state = as.numeric(gsub("(C|R)_", "", sn))
    circos.text(CELL_META$xcenter, CELL_META$ycenter, i_state, col = "white", 
                font = 2, cex = 0.7, adj = c(0.5, 0.5), niceFacing = TRUE)
    xlim = CELL_META$xlim
    breaks = seq(0, xlim[2], by = 4e5)
    circos.axis(major.at = breaks, labels = paste0(breaks/1000, "KB"), labels.cex = 0.5)
  }
}, bg.border = NA)


for(i in seq_len(nrow(cdm_res))) { 
  if(cdm_res$value[i] > 0) {
    circos.rect(cdm_res[i, "x1"], -uy(1, "mm"),
                cdm_res[i, "x1"] - abs(cdm_res[i, "value"]), -uy(2, "mm"),
                col = state_col2[cdm_res$cn[i]], border = state_col2[cdm_res$cn[i]], sector.index = cdm_res$rn[i], track.index = 2)
  }
}

abs_max = quantile(abs(c(meth_mat_1, meth_mat_2) - 0.5), 0.95, na.rm = TRUE)
col_fun = colorRamp2(c(0.5 - abs_max, 0.5, 0.5 + abs_max), c("blue", "white", "red"))
col_fun2 = colorRamp2(c(-abs_max, 0, abs_max), c("green", "white", "orange"))

ylim = get.cell.meta.data("ylim", sector.index = rownames(mat)[1], track.index = 1)
y1 = ylim[1] + (ylim[2] - ylim[1])*0.4
y2 = ylim[2]


for(i in seq_len(nrow(cdm_res))) {
  if(cdm_res$value[i] > 0) {
    circos.rect(cdm_res[i, "x1"], y1, cdm_res[i, "x1"] - abs(cdm_res[i, "value"]), y1 + (y2-y1)*0.45, 
                col = col_fun(meth_mat_1[cdm_res$rn[i], cdm_res$cn[i]]), 
                border = col_fun(meth_mat_1[cdm_res$rn[i], cdm_res$cn[i]]),
                sector.index = cdm_res$rn[i], track.index = 1)
    
    circos.rect(cdm_res[i, "x1"], y1 + (y2-y1)*0.55, cdm_res[i, "x1"] - abs(cdm_res[i, "value"]), y2, 
                col = col_fun2(meth_mat_2[cdm_res$rn[i], cdm_res$cn[i]] - meth_mat_1[cdm_res$rn[i], cdm_res$cn[i]]), 
                border = col_fun2(meth_mat_2[cdm_res$rn[i], cdm_res$cn[i]] - meth_mat_1[cdm_res$rn[i], cdm_res$cn[i]]),
                sector.index = cdm_res$rn[i], track.index = 1)
    
    circos.rect(cdm_res[i, "x2"], y1, cdm_res[i, "x2"] - abs(cdm_res[i, "value"]), y1 + (y2-y1)*0.45, 
                col = col_fun(meth_mat_2[cdm_res$rn[i], cdm_res$cn[i]]), 
                border = col_fun(meth_mat_2[cdm_res$rn[i], cdm_res$cn[i]]),
                sector.index = cdm_res$cn[i], track.index = 1)
    
    circos.rect(cdm_res[i, "x2"], y1 + (y2-y1)*0.55, cdm_res[i, "x2"] - abs(cdm_res[i, "value"]), y2, 
                col = col_fun2(meth_mat_1[cdm_res$rn[i], cdm_res$cn[i]] - meth_mat_2[cdm_res$rn[i], cdm_res$cn[i]]), 
                border = col_fun2(meth_mat_1[cdm_res$rn[i], cdm_res$cn[i]] - meth_mat_2[cdm_res$rn[i], cdm_res$cn[i]]),
                sector.index = cdm_res$cn[i], track.index = 1)
  }
}
circos.clear()



###basic used
set.seed(999)
mat = matrix(sample(18, 18), 3, 6) 
rownames(mat) = paste0("S", 1:3) 
colnames(mat) = paste0("E", 1:6) 
mat

df = data.frame(from = rep(rownames(mat), times = ncol(mat)), to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)
df

chordDiagram(mat)
chordDiagram(df)
circos.clear()


circos.par(gap.after = c(rep(5, nrow(mat)-1), 15, rep(5, ncol(mat)-1), 15)) 
chordDiagram(mat)
circos.clear()

circos.par(gap.after = c(rep(5, length(unique(df[[1]]))-1), 15, rep(5, length(unique(df[[2]]))-1), 15))
chordDiagram(df) 
circos.clear()

##colors
grid.col = c(S1 = "red", S2 = "green", S3 = "blue",
             E1 = "grey", E2 = "grey", E3 = "grey", E4 = "grey", E5 = "grey", E6 = "grey")
chordDiagram(mat, grid.col = grid.col) 
chordDiagram(t(mat), grid.col = grid.col)
chordDiagram(mat, grid.col = grid.col, transparency = 0)
col_mat = rand_color(length(mat), transparency = 0.5)
dim(col_mat) = dim(mat)
chordDiagram(mat, grid.col = grid.col, col = col_mat)

###
mat2 = matrix(sample(100, 35), nrow = 5) 
rownames(mat2) = letters[1:5] 
colnames(mat2) = letters[1:7]
mat2
chordDiagram(mat2, grid.col = 1:7, directional = 1, row.col = 1:5)

mat3 = mat2
for(cn in intersect(rownames(mat3), colnames(mat3))) {
  mat3[cn, cn] = 0
}
mat3

###Advanced usage
chordDiagram(mat3, grid.col = 1:7, directional = 1, row.col = 1:5)

chordDiagram(mat, grid.col = grid.col, annotationTrack = "grid") 
chordDiagram(mat, grid.col = grid.col, annotationTrack = c("name", "grid"), annotationTrackHeight = c(0.03, 0.01))
chordDiagram(mat, grid.col = grid.col, annotationTrack = NULL)

list(ylim = c(0, 1),
     track.height = circos.par("track.height"), bg.col = NA,
     bg.border = NA,
     bg.lty = par("lty"),
     bg.lwd = par("lwd"))

chordDiagram(mat, annotationTrack = NULL,
             preAllocateTracks = list(track.height = 0.3))
chordDiagram(mat, annotationTrack = NULL, preAllocateTracks = list(list(track.height = 0.1),
                                                                   list(bg.border = "black")))

chordDiagram(mat, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
circos.track(track.index = 1, panel.fun = function(x, y) { circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                                                                       facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5)) }, bg.border = NA)



factors = 1:20 # just indicate there are 20 sectors 
circos.par(gap.degree = 0, cell.padding = c(0, 0, 0, 0),
           start.degree = 360/20/2, track.margin = c(0, 0), clock.wise = FALSE)
circos.initialize(factors = factors, xlim = c(0, 1))

circos.track(ylim = c(0, 1), factors = factors, bg.col = "black", track.height = 0.15) 
circos.trackText(x = rep(0.5, 20), y = rep(0.5, 20), 
                 labels = c(13, 4, 18, 1, 20, 5, 12, 9, 14, 11, 8, 16, 7, 19, 3, 17, 2, 15, 10, 6), 
                 cex = 0.8, factors = factors, col = "#EEEEEE", font = 2, facing = "downward") 

circos.track(ylim = c(0, 1), factors = factors, 
             bg.col = rep(c("#E41A1C", "#4DAF4A"), 10), bg.border = "#EEEEEE", track.height = 0.05) 

circos.track(ylim = c(0, 1), factors = factors,
             bg.col = rep(c("black", "white"), 10), bg.border = "#EEEEEE", track.height = 0.275)

circos.track(ylim = c(0, 1), factors = factors,
             bg.col = rep(c("#E41A1C", "#4DAF4A"), 10), bg.border = "#EEEEEE", track.height = 0.05)



