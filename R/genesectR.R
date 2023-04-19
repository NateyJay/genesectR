
# require(stringr)



# Example -----------------------------------------------------------------


# master_set <- c(str_glue("Gene_{1:1000}"))
#
# ls <- list(Set_A= sample(master_set, 300),
#            Set_B= sample(master_set, 27),
#            Set_C= sample(master_set, 99))
#
#
# gs <- gs_import(ls, master_set)
# gs <- gs_compute_matricies(gs, mc=T)
# gs_plot_fischer(gs, breaks=3)



# Functions ---------------------------------------------------------------


#' Import function for gene list data
#'
#' @description Normalizes all input gene names to match the master list and produces a base gs object
#'
#' @param set_list List object of vectors containing gene names.
#' @param master_set Vector containing all gene names in the transcriptome.
#'
#' @return A list object which contains coded set information.
#' @export
#'
#' @examples
gs_import <- function(set_list, master_set) {
  library(stringr)

  message(str_glue("{length(set_list)} sets imported"))
  message(str_glue("{length(master_set)} entries in master set"))

  message("Removing entries from sets not in master")

  for (i in 1:length(set_list)) {
    s = set_list[[i]]
    before = length(s)
    s = match(s ,master_set)
    s = s[!is.na(s)]
    after = length(s)
    message(paste("     ", names(set_list)[i], " ",
                  before, " -> ", after, sep=''))
    set_list[[i]] = s
  }
  return(list(set_list=set_list, master_set=master_set))
}

#' fisher test
#'
#' @param setA
#' @param setB
#' @param master_len
#'
#' @return
#' @export
#'
#' @examples
gs_fisher <- function(setA, setB, master_len, master_set) {
  # Returns a list of stats between 2 sets of genes considering the length of all genes
  # This includes 1-tailed p-values for greater and less than, as well as the odd-ratio and size of overlap


  length(setA %in% master_set)
  length(setB %in% master_set)


  union = length(which(setA %in% setB))
  AnotB = length(which(setA %notin% setB))
  BnotA = length(which(setB %notin% setA))
  nots  = master_len - union - AnotB - BnotA

  df <- data.frame(A = c(union, AnotB), NotA=c(BnotA, nots))
  row.names(df) <- c("B", "NotB")

  g_test = fisher.test(df, alternative='g')
  g_test$p.value
  g_test$estimate

  l_test = fisher.test(df, alternative='l')
  l_test$p.value
  l_test$estimate

  ls = list(greater_p = g_test$p.value,
            lesser_p = l_test$p.value,
            odds_rat = as.numeric(g_test$estimate),
            overlap = df[1,1])
  return(ls)
}

#' montecarlo test
#'
#' @param setA
#' @param setB
#' @param master_set
#' @param n
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
gs_monte_carlo <- function(setA, setB, master_set, n= 1000, verbose=F) {
  # This randomly samples genes for our two sets and reports the size of their overlap (n=1000).
  # P-values are derived empirically as the proportion of random selections which where greater, or less than the known overlap.
  # Then, it calculates the sd and mean for the random-overlaps to estimate the Z-statistic for our known overlap relative to this random distribution.
  # It returns a list of the greater/lesser p-values as well as the z-stat.

  lenO = length(which(setA %in% setB))

  counter = 0

  ovs = c()

  for (i in 1:n) {
    if (i %% 100 == 0 & verbose) {
      message(i)
    }
    ranA = sample(master_set, length(setA))
    ranB = sample(master_set, length(setB))

    overlap = length(which(ranA %in% ranB))
    ovs = c(ovs, overlap)
    if (overlap < lenO) {
      counter = counter + 1
    }
  }

  x = ovs
  # qqnorm(x)
  # qqline(x)

  pop_sd <- sd(x)*sqrt((length(x)-1)/(length(x)))
  pop_mean <- mean(x)

  z_stat <- (lenO - pop_mean) / pop_sd

  # hist(ovs, xlim=c(0,700))
  # abline(v=lenO, col='red')

  # ?ztest

  greater_p = 1 - (counter / n )
  lesser_p = counter / n

  # g_test = binom.test(counter,n,alternative='g')
  # l_test = binom.test(counter,n,alternative='l')

  # z_score = corpora::z.score(counter, n)

  ls = list(greater_p = greater_p,
            lesser_p = lesser_p,
            z_score = z_stat)


}

#' Wrapper to compute all matricies
#'
#' @param gs
#' @param mc
#'
#' @return
#' @export
#'
#' @examples
gs_compute_matricies <- function(gs, mc=F) {
  # pval : FET p-values
  # padj : FET adjusted p-values (BH)
  # odds : FET odds-ratio
  # overlap : size of known overlap
  # mc_pval : montecarlo p-values
  # m_z : montecarlo z-stats

  comp.df <- data.frame()

  message("  -> performing pairwise tests")
  pairs = as.data.frame(t(combn(1:length(gs$set_list),2)))
  for (i in 1:nrow(pairs)) {
    message(".", appendLF=F)
    # print(pairs[i,])
    a_i = pairs[i,1]
    b_i = pairs[i,2]
    setA = gs$set_list[[a_i]]
    setB = gs$set_list[[b_i]]

    fi_test = gs_fisher(setA, setB, length(gs$master_set), gs$master_set)


    add.df <- data.frame(A=names(gs$set_list[a_i]),
                         B=names(gs$set_list[b_i]),
                         over=     fi_test$greater_p,
                         under=    fi_test$lesser_p,
                         odds=     fi_test$odds_rat,
                         overlap=  fi_test$overlap)

    if (mc) {
      mc_test = gs_monte_carlo(setA, setB, gs$master_set, n=1000)

      add.df <- cbind(add.df,
                      data.frame(
                        mc_over=  mc_test$greater_p,
                        mc_under= mc_test$lesser_p,
                        z_score=  mc_test$z_score
                        )
                      )

    }
    comp.df <- rbind(comp.df, add.df)

  }
  message(" done")
  #
  # # gtools::combinations(4,2,names(set_list))
  # m = matrix(data = rep(NA, length(gs$set_list)^2), ncol=length(gs$set_list))
  # dimnames(m) <- list(names(set_list), names(set_list))
  #
  #

  message("  -> forming matricies")
  ## now building a beautiful matrix
  i = 0
  vec = c()
  odds_vec = c()
  lap_vec = c()
  mc_vec = c()
  z_vec = c()
  for (col in names(gs$set_list)) {
    tri='upper'
    for (row in names(gs$set_list)) {
      i = i + 1

      # print(paste(i,col,row))

      # message(paste(i,col,row, sep=' -- '))
      if (col == row) {
        vec = c(vec, 0)
        odds_vec = c(odds_vec, 1)
        lap_vec = c(lap_vec, length(gs$set_list[[col]]))
        mc_vec = c(mc_vec, 0)
        z_vec = c(z_vec, 1)
        tri='lower'

      } else {

        comp = comp.df[comp.df$A == col & comp.df$B == row,]

        if (nrow(comp) == 0){
          comp = comp.df[comp.df$B == col & comp.df$A == row,]
        }

        if (nrow(comp) > 1) {
          stop("ERROR")
        }

        if (tri == 'upper') {
          vec = c(vec, comp$over)
          mc_vec = c(mc_vec, comp$mc_over)
        } else {
          vec = c(vec, comp$under)
          mc_vec = c(mc_vec, comp$mc_under)
        }
        odds_vec = c(odds_vec, comp$odds)
        lap_vec = c(lap_vec, comp$overlap)
        z_vec = c(z_vec, comp$z_score)
      }

    }
  }

  message("  -> cleaning up and exporting")
  m = matrix(vec, byrow=F, ncol=length(gs$set_list))
  m_adj <- matrix(p.adjust(vec, 'BH'), byrow=F, ncol=length(gs$set_list))
  m_odds = matrix(odds_vec, byrow=F, ncol=length(gs$set_list))
  m_lap  = matrix(lap_vec, byrow=F, ncol=length(gs$set_list))
  if (mc) {
    m_mc  <- matrix(mc_vec, byrow=F, ncol=length(gs$set_list))
    m_z    = matrix(z_vec, byrow=F, ncol=length(gs$set_list))
  }

  namify <- function(m) {
    rownames(m) <- names(gs$set_list)
    colnames(m) <- names(gs$set_list)
    diag(m) <- NA
    return(m)
  }

  m <- namify(m)
  m_odds <- namify(m_odds)
  m_lap <- namify(m_lap)
  m_adj <- namify(m_adj)
  if (mc) {
    m_mc <- namify(m_mc)
    m_z <- namify(m_z)
  }

  ls <- list(fet_pval=m,
             fet_padj=m_adj,
             fet_odds=m_odds,
             overlap=m_lap)
  if (mc) {
    ls$mc_pval= m_mc
    ls$mc_z=m_z
  }

  gs$mats <- ls
  return(gs)
}

#' Basic plotting structure
#'
#' @param gs
#' @param breaks
#' @param plot_type
#'
#' @return
#' @export
#'
#' @examples
gs_make_plot.df <- function(gs, breaks, plot_type="pval") {
  # This function produces a dataframe of mapping dimensions and colors for all of the resulting plots.
  # Plots are produced in a blank frame with the base dimension of a box determined by "box_size".
  # The plotting frames dimensions are determined by the box_size and number of comparisons

  set_list = gs$set_list
  mats = gs$mats
  box_size = 1

  tot_len = length(set_list)

  dim = tot_len
  gs$dim = dim

  plot.df <- data.frame(x_name= rep(names(set_list), each=tot_len),
                        y_name= rep(names(set_list), times=tot_len),
                        x1 = rep(0:(tot_len-1), each=tot_len)*box_size,
                        y1 = rep((tot_len-1):0, times=tot_len)*box_size)


  plot.df$x2 = plot.df$x1 + box_size
  plot.df$y2 = plot.df$y1 + box_size

  plot.df$text.x = (plot.df$x2 + plot.df$x1) / 2
  plot.df$text.y = (plot.df$y2 + plot.df$y1) / 2
  plot.df$pval = as.vector(mats$fet_pval) #round(as.vector(mats$fet_pval),3)
  filter = (plot.df$pval >= 0.01 & !is.na(plot.df$pval))
  plot.df$pval[filter] <- round(plot.df$pval[filter],2)
  filter = (plot.df$pval < 0.01 & !is.na(plot.df$pval))
  plot.df$pval[filter] <- formatC(plot.df$pval[filter], format = "E", digits = 1)
  # plot.df$pval[plot.df$pval < 0.001] <- "<0.001"

  plot.df$mc <- as.vector(mats$mc_pval) #round(as.vector(mats$mc_pval),3)
  filter = (plot.df$mc >=0.001 & !is.na(plot.df$mc))
  plot.df$mc[filter] <- round(plot.df$mc[filter], 2)
  plot.df$mc[plot.df$mc < 0.001] <- "<0.001"

  plot.df$z <- as.vector(mats$mc_z) #round(as.vector(mats$mc_z), 1)
  plot.df$z <- round(plot.df$z, 2)

  plot.df$odds = as.vector(mats$fet_odds) #round(as.vector(mats$fet_odds),2)
  plot.df$odds <- round(plot.df$odds, 2)

  plot.df$overlap <- as.vector(mats$overlap)


  ### Making text block fields
  make_text <- function(v) {
    out <- paste(v, "\n(",plot.df$overlap, ")", sep='')
    out[is.na(v)] <- NA
    return(out)
  }

  plot.df$text_pval <- make_text(plot.df$pval)
  plot.df$text_odds <- make_text(plot.df$odds)
  plot.df$text_mc <- make_text(plot.df$mc)
  plot.df$text_z <- make_text(plot.df$z)

  ## Filtering double-shown data for odds-ratios
  filter = (!is.na(plot.df$odds) &
              plot.df$odds < 1 &
              as.vector(upper.tri(mats$fet_odds)))
  plot.df$text_odds[filter] <- NA

  filter = (!is.na(plot.df$odds) &
              plot.df$odds > 1 &
              as.vector(!upper.tri(mats$fet_odds)))
  plot.df$text_odds[filter] <- NA

  ## Filtering double-shown data for z-statistic
  filter = (!is.na(plot.df$z) &
              plot.df$z < 0 &
              as.vector(upper.tri(mats$mc_z)))
  plot.df$text_z[filter] <- NA
  plot.df$z[filter] <- NA

  filter = (!is.na(plot.df$z) &
              plot.df$z > 0 &
              as.vector(!upper.tri(mats$mc_z)))
  plot.df$text_z[filter] <- NA
  plot.df$z[filter] <- NA


  ## Making colors for fisher pvalue text
  plot.df$col_pval = "white"
  plot.df$col_pval[as.vector(mats$fet_pval) < 0.01 &
                     as.vector(upper.tri(mats$fet_pval))] <- "yellow"
  plot.df$col_pval[as.vector(mats$fet_pval) < 0.01 &
                     as.vector(!upper.tri(mats$fet_pval))] <- "darkolivegreen1"
  plot.df$col_pval[is.na(as.vector(mats$fet_pval))] <- "grey60"


  ## Making colors for montecarlo pvalue text
  plot.df$col_mc = "white"
  plot.df$col_mc[as.vector(mats$mc_pval) < 0.01 &
                   as.vector(upper.tri(mats$mc_pval))] <- "coral"
  plot.df$col_mc[as.vector(mats$mc_pval) < 0.01 &
                   as.vector(!upper.tri(mats$mc_pval))] <- "cadetblue"
  plot.df$col_mc[is.na(as.vector(mats$mc_pval))] <- "grey60"


  shadify <- function(od, pal='Reds 3', log=F) {
    if (length(od) == 0) {
      return("grey60")
    }
    od = abs(od)
    if (log) {
      od = log(od,2)
    }
    limit = max(od, na.rm=T)
    scale = hcl.colors(10, palette=pal, rev=T)
    ab = abs(od)
    i_vec = round(ab / max(ab, na.rm=T) * 9 ) + 1
    shade = scale[i_vec]
    return(shade)
  }

  # shade_and_legend_z <- function(x, top_pal, bot_pal) {
  #
  #   maximum = max(abs(x), na.rm=T)
  #
  #
  #   filter = (x > 0 & !is.na(x))
  #   log10(x[filter]) hcl.colors(10, top_pal)
  #   axisTicks(c(0,1900), log=F)
  #
  # }


  ### Making colors for z-scores
  plot.df$col_z = "white"
  filter = (!is.na(plot.df$z) &
              plot.df$z > 0 &
              as.vector(upper.tri(mats$mc_z)))
  plot.df$col_z[filter] = shadify(plot.df$z[filter], log=T, pal="Blues 2")
  filter = (!is.na(plot.df$z) &
              plot.df$z < 0 &
              as.vector(!upper.tri(mats$mc_z)))
  plot.df$col_z[filter] = shadify(plot.df$z[filter], log=T, pal= 'Purples 2')
  plot.df$col_z[is.na(as.vector(mats$mc_z))] <- "grey60"


  ### Making colors for odds ratios
  plot.df$col_odds = "white"
  filter = (!is.na(plot.df$odds) &
              plot.df$odds > 1 &
              as.vector(upper.tri(mats$fet_odds)))
  plot.df$col_odds[filter] = shadify(plot.df$odds[filter], log=T, pal="Greens 2")
  filter = (!is.na(plot.df$odds) &
              plot.df$odds < 1 &
              as.vector(!upper.tri(mats$fet_odds)))
  plot.df$col_odds[filter] = shadify(plot.df$odds[filter], log=T, pal= 'Reds 3')

  plot.df$col_odds[is.na(as.vector(mats$fet_odds))] <- "grey60"
  gs$plot.df <- plot.df
  return(gs)
}


#' Wrapper to plot the fisher exact test analysis
#'
#' @param gs
#' @param breaks
#' @param gap
#' @param margin
#'
#' @return
#' @export
#'
#' @examples
gs_plot_fischer <- function(gs, breaks = c(), gap=0.2, margin=5) {
  set_list = gs$set_list
  mats = gs$mats

  round_and_na_fix <- function(x,r=2) {
    x = round(x,r)
    x[is.na(x)] <- ''
    return(x)
  }

  rm_non_plots <- function(df) {
    n = unique(df$x_name)
    c <- gtools::combinations(length(n),2,n)
    for (i in 1:nrow(c)) {
      a = df$padj[df$x_name == c[i,1] & df$y_name == c[i,2]]
      b = df$padj[df$x_name == c[i,2] & df$y_name == c[i,1]]

      # if (df$x_name == c[i,1] == df$y_name == c[i,2]) {
      #   df$to_print[df$x_name == c[i,1] & df$y_name == c[i,2]] <- F

      if (a > b) {
        df$to_print[df$x_name == c[i,1] & df$y_name == c[i,2]] <- F
      } else if (b > a) {
        df$to_print[df$x_name == c[i,2] & df$y_name == c[i,1]] <- F
      } else {
        df$to_print[df$x_name == c[i,2] & df$y_name == c[i,1]] <- T
      }
    }
    return(df)
  }

  p.df <- data.frame(x_name=rep(names(set_list), each=length(set_list)),
                     y_name=rep(names(set_list), length(set_list)))

  p.df$padj <- as.vector(mats$fet_padj)
  p.df$odds <- as.vector(mats$fet_odds)
  p.df$lodds  <- log2(p.df$odds)
  p.df$inter <- as.vector(mats$overlap)
  p.df$to_print <- T
  p.df <- rm_non_plots(p.df)



  p.df$xi <- match(p.df$x_name, names(set_list)) - 1
  p.df$yi <- match(rev(p.df$y_name), names(set_list)) - 1


  p.df$x <- p.df$xi + (cumsum(1:length(set_list) %in% breaks) * gap)[match(p.df$x_name, names(set_list))]
  p.df$y <- p.df$yi - (cumsum(1:length(set_list) %in% breaks) * gap)[match(p.df$y_name, names(set_list))]
  p.df$y <- p.df$y - min(p.df$y, na.rm=T) + 0.5


  p.df$col = 'white'
  p.df$col[p.df$x_name == p.df$y_name] <- 'grey50'
  # hcl.pals(type='diverging')
  max_odds=max(abs(p.df$lodds[is.finite(p.df$lodds)]), na.rm=T)
  palette = hcl.colors(100, 'Purple-Green')[20:80]
  filter = (!is.na(p.df$lodds) & p.df$to_print)
  p.df$col[filter] <- range_to_color(p.df$lodds[filter],
                                     pal=palette,
                                     minmax=c(-max_odds, max_odds))

  par(mar=c(2,margin,margin,2), xpd=T)



  plot(1,type='n', xlab='', ylab='', xlim=c(0,max(p.df$x+1)), ylim=c(0,max(p.df$x+1)), axes=F)

  rect(p.df$x, p.df$y, p.df$x + 1, p.df$y+1, col=p.df$col)

  t.df <- p.df[p.df$to_print & !is.na(p.df$inter),]
  # t.df <- p.df

  text(t.df$x+.5, t.df$y+.5, round_and_na_fix(t.df$lodds), pos=3, font=2)
  text(t.df$x+.5, t.df$y+.5, round_and_na_fix(log10(t.df$padj)),adj=c(0.5,0.5))
  text(t.df$x+.5, t.df$y+.5, str_glue("({round_and_na_fix(t.df$inter)})"), pos=1, font=3)

  text(unique(p.df$x)+.5, max(p.df$y,na.rm=T)+1, unique(p.df$x_name),pos=3, font=3)

  text(0, unique(p.df$y)+.5, unique(p.df$y_name),pos=2, font=3)

  lens = unlist(lapply(gs$set_list, length))
  filter = (p.df$x_name == p.df$y_name)
  text(p.df$x[filter]+.5, p.df$y[filter]+.5, str_glue("({lens})"), font=3)

  # text(unique(p.df$x)+.5, max(p.df$y,na.rm=T)+1, str_glue("{unique(p.df$x_name)}\n({unlist(lapply(gs$set_list, length))})"),pos=3, font=3)
  #
  # text(0, unique(p.df$y)+.5, str_glue("{unique(p.df$y_name)}\n({unlist(lapply(gs$set_list, length))})"),pos=2, font=3)


  usr = par()$usr

  color_legend(0, 0,palette, dim.x = 1, .25, minmax=c(-round(max_odds,1),round(max_odds,1)),
               main='L2-fold odds-ratio')

}


#' Multi-comparison plot
#'
#' @param gs
#' @param expand
#' @param value.cex
#' @param label.cex
#'
#' @return
#' @export
#'
#' @examples
gs_multi_plot <- function(gs, expand, value.cex=0.8, label.cex=0.8) {
  set_list = gs$set_list

  x_len = length(set_list) - length(expand)

  if (x_len < 1) {
    stop("error: set_list - expand must be > 0")
  }

  mat = gs$mats$fet_padj
  m.df = reshape2::melt(mat)
  names(m.df) <- c("rows", "cols", "padj")

  order_paste <- function(x) {
    x = x[order(x)]
    return(paste(x, collapse='-'))
  }

  m.df$comp = apply(m.df[,c(1,2)], 1, order_paste)

  direction <- function(x) {
    n = colnames(mat)
    if (x[1] == x[2]) {
      return(NA)
    }

    i1 = match(x[1],n)
    i2 = match(x[2],n)
    if (i2 < i1) {
      return('lower')
    } else {
      return('upper')
    }
  }

  m.df$tri = apply(m.df[,c(1,2)], 1, direction)

  m.df <- m.df[order(m.df$padj),]
  m.df <- m.df[!duplicated(m.df$comp),]
  m.df <- m.df[complete.cases(m.df),]

  x_names = colnames(mat)[colnames(mat) %notin% expand]
  y_names = expand

  p.df <- data.frame(x_name=rep(x_names, length(y_names)),
                     y_name=rep(y_names, each=length(x_names)))
  p.df$comp <- apply(p.df[,c(1,2)], 1, order_paste)
  p.df$pval <- m.df$padj[match(p.df$comp, m.df$comp)]

  mat_to_val = function(x, mat) {
    return(mat[x[1], x[2]])
  }
  p.df$odds <- apply(p.df[c(1,2)], 1, mat_to_val, gs$mats$fet_odds)
  p.df$padj <- m.df$padj[match(p.df$comp, m.df$comp)]
  p.df$overlap <- apply(p.df[c(1,2)], 1, mat_to_val, gs$mats$overlap)
  p.df$x <- match(p.df$x_name, unique(p.df$x_name))
  p.df$y <- match(p.df$y_name, unique(p.df$y_name))

  max_dim = max(c(p.df$y,p.df$x)) + 1

  p.df$x = max_dim - max(p.df$x) + p.df$x - 1
  p.df$y = max_dim - max(p.df$y) + p.df$y - 1

  p.df$log_odds <- log2(p.df$odds)
  p.df$log_padj <- log10(p.df$padj)


  palette = hcl.colors(100, 'Purple-Green')[20:80]


  minmax=abs(c(max(p.df$log_odds, na.rm=T), min(p.df$log_odds, na.rm=T)))
  minmax = c(-1*minmax, minmax)
  p.df$fill_color <- range_to_color(p.df$log_odds, palette, minmax=minmax)



  par(xpd=T, mar=rep(5,4))
  plot(1,1, type='n', xlim=c(0,max_dim), ylim=c(0,max_dim), xlab='', ylab='', axes=F)
  rect(p.df$x, p.df$y, p.df$x+1, p.df$y+1, col=p.df$fill_color)

  text(p.df$x+.5, p.df$y+.5, round(p.df$log_odds,1), pos=3,# adj=c(0.5,-0.3),
       cex=value.cex, font=2)
  text(p.df$x+.5, p.df$y+.5, round(p.df$log_padj,1),# adj=c(0.5,1.3),
       cex=value.cex, font=1)
  text(p.df$x+.5, p.df$y+.5, str_glue("({round(p.df$overlap,1)})"), pos=1, # adj=c(0.5,1.3),
       cex=value.cex, font=3)

  text(min(p.df$x) - 0.5, unique(p.df$y)+0.5, unique(p.df$y_name), adj=c(1,0.5), cex=label.cex, font=3)
  text(unique(p.df$x)+0.5, max(p.df$y) + 1.5, unique(p.df$x_name), adj=c(0,0.5), srt=90, cex=label.cex, font=3)

  x1 = unique(p.df$x); x2 = x1 + 1; y1 = max(p.df$y)+1.1; y2=y1 + .25
  rect(x1,y1,x2,y2, col='grey85', border='black')
  text((x1+x2)/2, (y1+y2)/2, lengths(gs$set_list)[x_names], cex=label.cex, font=3)


  x1 = min(p.df$x) - .35; x2 = x1 + .25; y1 = unique(p.df$y); y2=y1 + 1
  rect(x1,y1,x2,y2, col='grey85', border='black')
  text((x1+x2)/2, (y1+y2)/2, lengths(gs$set_list)[y_names], cex=label.cex, font=3, srt=90)



  break_n = 21

  minmax=round(max(abs(p.df$log_odds[is.finite(p.df$log_odds)]), na.rm=T),1)


  color_legend(min(p.df$x), 0,palette, dim.x = 1, .25, minmax=c(-round(minmax,1),round(minmax,1)),
               main='L2-fold odds-ratio')



}




# master_set <- str_glue("Gene_{1:1000}")
#
# ls <- list(Set_A= sample(master_set, 300),
#            Set_B= sample(master_set[1:100], 27),
#            Set_C= sample(master_set, 99),
#            Set_V= sample(master_set[1:100], 15),
#            Set_W= sample(master_set, 201),
#            Set_X= sample(master_set, 500),
#            Set_Y= sample(master_set, 44),
#            Set_Z= sample(master_set, 766))
#
#
# gs <- gs_import(ls, master_set)
# gs <- gs_compute_matricies(gs)
#
#
# gs_multi_plot(gs, c("Set_V","Set_W","Set_X","Set_Y","Set_Z"))




#' palette function
#'
#' @param values
#' @param pal
#' @param minmax
#' @param log
#'
#' @return
#' @export
#'
#' @examples
range_to_color <- function(values, pal, minmax=NULL, log=F) {
  values = unlist(values)

  if (log) {
    values = log10(values)
  }

  if (is.null(minmax)) {
    minmax = c(min(values[is.finite(values)]), max(values[is.finite(values)]))
  }

  range = minmax[2] - minmax[1]

  norm = (values - minmax[1]) / range
  norm[norm < 0] <- 0
  norm[norm > 1] <- 1

  idx  = norm * (length(pal) - 1) + 1
  idx = round(idx)

  return(pal[idx])
}


#' Legend plotting function
#'
#' @param ori.x
#' @param ori.y
#' @param colors
#' @param minmax
#' @param dim.x
#' @param dim.y
#' @param main
#' @param cex
#'
#' @return
#' @export
#'
#' @examples
color_legend <- function(ori.x,ori.y, colors, minmax, dim.x, dim.y, main='', cex=.8) {
  x2 = ori.x + dim.x
  for (i in 1:length(colors)) {
    x1 = ori.x + ((i-1)/(length(colors)-1)) * dim.x
    rect(x1,ori.y,x2,ori.y+dim.y,border=NA,col=colors[i])
  }
  rect(ori.x, ori.y, ori.x + dim.x, ori.y + dim.y)
  text(ori.x, ori.y + dim.y/2, minmax[1], pos=2, cex=cex)#adj=c(1.2,0.5))
  text(ori.x + dim.x, ori.y + dim.y/2, minmax[2], pos=4, cex=cex)# adj=c(-.2,0.5))
  text(ori.x + dim.x/2, ori.y, main, pos=1, cex=cex)


}


`%notin%` <- Negate(`%in%`)

