library(sequenza)



ratio.bayes <- function(depth.ratio
                        , cellularity
                        , ploidy
                        , avg.depth.ratio
                        ,weight.ratio = 100
                        , CNt.min = 0
                        ,CNt.max = 7
                        , CNn = 2
                        , priors.table = data.frame(CN = CNt.min:CNt.max,
                       value = 1)) {

   mufreq.tab <- data.frame(ratio = depth.ratio,
                            weight.ratio = weight.ratio)
   mufreq.depth.ratio <- model.points(cellularity = cellularity, ploidy = ploidy,
                                      types = cbind(CNn = CNn, CNt = CNt.min:CNt.max, Mt = 0),
                                      avg.depth.ratio = avg.depth.ratio)
   model.d.ratio      <- cbind(CNt = CNt.min:CNt.max, depth.ratio = mufreq.depth.ratio[, 2])

   model.pts          <- as.data.frame(model.d.ratio)
   rows.x             <- 1:nrow(mufreq.tab)

   priors <- rep(1, nrow(model.pts))
   for (i in 1:nrow(priors.table)) {
      priors[model.pts$CNt == priors.table$CN[i]] <- priors.table$value[i]
   }
   priors <- priors / sum(priors)

   bayes.fit <- function (x, mat, model.pts, priors) {
      test.ratio <- model.pts$depth.ratio
      min.offset <- 1e-323
      score.r    <- sequenza:::depth.ratio.dbinom(size = mat[x,]$weight.ratio, depth.ratio = mat[x,]$ratio, test.ratio)

      score.r    <- score.r * priors

      post.model <- score.r

      post.model[post.model == 0] <- min.offset

         max.lik <-  which.max(post.model)
         max.post <- c(as.numeric(model.pts[max.lik,1]), log2(post.model[max.lik]))

      max.post
   }
   bafs.L           <- mapply(FUN = bayes.fit, rows.x,
                         MoreArgs = list(mat = mufreq.tab,
                                         model.pts = model.pts,
                                         priors = priors),
                                         SIMPLIFY = FALSE)
   bafs.L           <- do.call(rbind, bafs.L)
   colnames(bafs.L) <- c("CNt", "L")
   bafs.L
}

ratio.model.fit <- function(cellularity = seq(0, 1, by = 0.01),
                          ploidy = seq(1, 7, by = 0.1),
                          mc.cores = getOption("mc.cores", 2L), ...) {

   result <- expand.grid(ploidy = ploidy, cellularity = cellularity,
                         KEEP.OUT.ATTRS = FALSE)

   fit.cp <- function(ii) {
      L.model <- ratio.bayes(cellularity = result$cellularity[ii],
                           ploidy = result$ploidy[ii], ...)
      sum(L.model[,2])
   }
   bayes.res <- sequenza:::mclapplyPb(X = 1:nrow(result), FUN = fit.cp, mc.cores = mc.cores)
   #bayes.res <- lapply(X = 1:nrow(result), FUN = fit.cp)
   result$L <- unlist(bayes.res)
   z <- tapply(result$L, list(result$ploidy, result$cellularity), mean)
   x <- as.numeric(rownames(z))
   y <- as.numeric(colnames(z))
   max.lik <- max(result$L)
   LogSumLik <- log2(sum(2^(result$L - max.lik))) + max.lik
   znorm <- 2^(z - LogSumLik)
  
   list(ploidy = x, cellularity = y, lpp = znorm)
}


genome.view.logR <- function(seg.cn, ...) {
   chr.order <- unique(seg.cn$chromosome)
   #seg.list  <- split(x = seg.cn[,c("chromosome", "start.pos", "end.pos", "mean")],
   #                   f = seg.cn$chromosome)
   seg.list  <- split(x = seg.cn[,c("chromosome", "start.pos", "end.pos", "log2")],
                      f = seg.cn$chromosome)
   seg.list  <- seg.list[order(order(chr.order))]
   seg.max   <- lapply(X = seg.list, FUN = function(x) x[nrow(x), "end.pos" ])
   seg.pos   <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
   seg.max   <- cumsum(as.numeric(do.call(rbind, seg.max)))
   chr.offset <- 0
   for (i in 1:length(seg.pos)){
      seg.pos[[i]] <- seg.pos[[i]] + chr.offset
      colnames(seg.pos[[i]]) <- c("abs.start","abs.end")
      chr.offset   <- seg.max[i]
   }
   seg.max      <- sapply(X = seg.pos, FUN = function(x) x[nrow(x), "abs.end" ])
   abs.list     <- mapply(cbind, seg.list, seg.pos, SIMPLIFY = FALSE)
   abs.segments <- do.call(rbind, abs.list)
      plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
           y = c(min(abs.segments$log2), max(abs.segments$log2)), type = "n",
           ylab = "logR", xlab = "Position (Mb)",
           xaxt='n', yaxt = 'n', xaxs = "i", ...)
      axis(labels = min(abs.segments$log2):max(abs.segments$log2),
           at = min(abs.segments$log2):max(abs.segments$log2),
           side = 2, line = 0, las = 1)
      #abline(h = c(min(abs.segments$mean):max(abs.segments$mean)), lty = 2)
      segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
               y0 = abs.segments$log2, y1= abs.segments$log2, col="red", lwd = 5, lend = 1)

   abline(v = c(0, seg.max), lty = 3)
   for (i in 1:length(abs.list)){
      max.pos <- nrow(abs.list[[i]])
      mtext(chr.order[i], side = 3, line = 0,
            at = sum(abs.list[[i]]$abs.start[1], abs.list[[i]]$abs.end[max.pos])/2)
      #axis(labels = as.character(round(seq(abs.list[[i]]$start.pos[1]/1e6, abs.list[[i]]$end.pos[max.pos]/1e6, by = 20), 0)),
      #     at = seq(abs.list[[i]]$abs.start[1], abs.list[[i]]$abs.end[max.pos], by = 2e7), outer = FALSE, cex = par("cex.axis")*par("cex"),
      #     side = 1 , line = 0)
   }
   axis(labels = as.character(round(seq(abs.list[[1]]$start.pos[1]/1e6, abs.list[[1]]$end.pos[nrow(abs.list[[1]])]/1e6, by = 50), 0)),
        at = seq(abs.list[[1]]$abs.start[1], abs.list[[1]]$abs.end[nrow(abs.list[[1]])], by = 5e7), outer = FALSE, cex = par("cex.axis")*par("cex"),
        side = 1 , line = 1)
}

copynumber.view <- function(seg.cn, info.type = "AB", ...) {
   chr.order <- unique(seg.cn$chromosome)
   seg.list  <- split(x = seg.cn[,c("chromosome", "start.pos", "end.pos", "A", "B", "CNt")],
                      f = seg.cn$chromosome)
   seg.list  <- seg.list[order(order(chr.order))]
   seg.max   <- lapply(X = seg.list, FUN = function(x) x[nrow(x), "end.pos" ])
   seg.pos   <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
   seg.max   <- cumsum(as.numeric(do.call(rbind, seg.max)))
   chr.offset <- 0
   for (i in 1:length(seg.pos)){
      seg.pos[[i]] <- seg.pos[[i]] + chr.offset
      colnames(seg.pos[[i]]) <- c("abs.start","abs.end")
      chr.offset   <- seg.max[i]
   }
   seg.max      <- sapply(X = seg.pos, FUN = function(x) x[nrow(x), "abs.end" ])
   abs.list     <- mapply(cbind, seg.list, seg.pos, SIMPLIFY = FALSE)
   abs.segments <- do.call(rbind, abs.list)
   if (info.type == "AB") {
      abs.segments <- na.exclude(abs.segments)
      plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
           y = c(-0.1, (max(abs.segments$A)+0.1)), type = "n",
           ylab = "Copy number", xlab = "Position (Mb)",
           xaxt='n',  yaxt='n', xaxs = "i", ...)
      axis(labels = 0:max(abs.segments$A),
           at = 0:max(abs.segments$A),
           side = 2, line = 0, las = 1)
      #abline(h = c(0:max(abs.segments$A)), lty = 2)
      segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
               y0 = (abs.segments$B-0.1), y1 = (abs.segments$B-0.1), col="blue", lwd = 5, lend = 1)
      segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
               y0 = (abs.segments$A+0.1), y1 = (abs.segments$A+0.1), col="red", lwd = 5, lend = 1)
   } else {
      abs.segments <- abs.segments[!is.na(abs.segments$CNt), ]
      plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
           y = c(min(abs.segments$CNt), max(abs.segments$CNt)), type = "n",
           ylab = "Copy number", xlab = "Position (Mb)",
           xaxt='n', yaxt = 'n', xaxs = "i", ...)
      axis(labels = min(abs.segments$CNt):max(abs.segments$CNt),
           at = min(abs.segments$CNt):max(abs.segments$CNt),
           side = 2, line = 0, las = 1)
      #abline(h = c(min(abs.segments$CNt):max(abs.segments$CNt)), lty = 2)
      segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
               y0 = abs.segments$CNt, y1= abs.segments$CNt, col="red", lwd = 5, lend = 1)
   }
   abline(v = c(0, seg.max), lty = 3)
   for (i in 1:length(abs.list)){
      max.pos <- nrow(abs.list[[i]])
      mtext(chr.order[i], side = 3, line = 0, las=2,
            at = sum(abs.list[[i]]$abs.start[1], abs.list[[i]]$abs.end[max.pos])/2)
      #axis(labels = as.character(round(seq(abs.list[[i]]$start.pos[1]/1e6, abs.list[[i]]$end.pos[max.pos]/1e6, by = 20), 0)),
      #     at = seq(abs.list[[i]]$abs.start[1], abs.list[[i]]$abs.end[max.pos], by = 2e7), outer = FALSE, cex = par("cex.axis")*par("cex"),
      #     side = 1 , line = 0)
   }
   axis(labels = as.character(round(seq(abs.list[[1]]$start.pos[1]/1e6, abs.list[[1]]$end.pos[nrow(abs.list[[1]])]/1e6, by = 50), 0)),
        at = seq(abs.list[[1]]$abs.start[1], abs.list[[1]]$abs.end[nrow(abs.list[[1]])], by = 5e7), outer = FALSE, cex = par("cex.axis")*par("cex"),
        side = 1 , line = 1)
}

