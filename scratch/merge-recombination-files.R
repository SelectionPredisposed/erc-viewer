
setwd('/home/oyvind/gitprojects/erc-viewer/resources/CEU/')
rec <- read.table(gzfile('CEU-1-final.txt.gz'), header = T, stringsAsFactors = F)

ggplot(rec) + geom_line(aes(Position.bp., Rate.cM.Mb.))

# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/
dflist = list()
for(i in seq(1,22,1)){
  df <- read.table(gzfile(paste0('CEU-',i,'-final.txt.gz')), header = T, stringsAsFactors = F)
  df$chr <- i
  assign(x = paste0("rr",i), value = df)
}

rrlist <- sprintf("rr%s",seq(1,22,1))
x.list <- lapply(rrlist, get)

rr.all <- do.call(rbind, x.list)

ggplot(subset(rec, Position.bp. < 1000000)) + geom_line(aes(Position.bp., Rate.cM.Mb.))
