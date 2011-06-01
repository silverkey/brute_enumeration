source('analyze_initiator.R')
library(RMySQL)

org.code = 'MM';
wl = 7

db = paste(sep='',org.code,'_promoter_20_20_NSH_100_WL_',wl)
usr = 'mysql_dev'
pwd = 'riiGbs'
zcutoff = 2.5
stat = paste(sep='','stat_',org.code,'_promoter_20_20.fa_100_',wl,'_counts.xls.csv')
pdfname = paste(sep='',org.code,'_',wl,'bp.pdf')
drv = dbDriver("MySQL")
con = dbConnect(drv,db,usr,pwd)
t = read.table(file = stat, sep = ',', header = T)
t = na.omit(t[t$z.seq > zcutoff,])
pos = analyze.words(t$word,con)
pos = na.omit(pos[pos$z.score > zcutoff,])
pos = pos[order(pos$median.diff,decreasing=T),]
draw.words(pos$word,con,pdfname)
