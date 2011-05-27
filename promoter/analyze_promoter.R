# This script works on the output files from the
# brute_enumeration_mysql.pl script

analyze.words = function(w,con,name) {
  tab = data.frame(row.names=seq(1:length(w)))
  zscore = c()
  topcount = c()
  mediandiff = c()
  for(i in 1:length(w)) {
  
    word = w[i]
    query = paste(sep='',"SELECT id FROM word WHERE word = '",word,"'")
    res = dbSendQuery(con,query)
    id = as.numeric(fetch(res))

    query = paste(sep='','SELECT pos FROM seq_occurrence WHERE word_id = ',id)
    res = dbSendQuery(con,query)
    pos = fetch(res,n=-1)
    pos = pos$pos

    # WORK ON COUNTS OF THE HISTOGRAM BREAKS
    h = hist(pos,breaks=seq(1,100,1),plot=F)
    # GET THE AVERAGE Z.SCORE FROM THE MAXIMUM AND THE LEFT AND THE RIGTH BREAKS
    # THIS SEEMS TO BE MORE EFFECTIVE THAN TO WORK ONLY ON: 
    # z.score = (max(h$counts)-mean(h$counts))/sd(h$counts)
    z.score = supported.z(h)
    # THE POSITION OF THE PEAK
    top.count = h$breaks[h$counts==max(h$counts)][1]
    # OTHER EFFECTIVE MEASURE: THE DIFFERENCE BETWEEN THE COUNTS OF THE MAXIMUM
    # AND THE MEDIAN OF THE DISTRIBUTION OF THE COUNTS FOR THIS WORD 
    median.diff = max(h$counts)-median(h$counts)
    
    zscore = c(zscore,z.score)
    topcount = c(topcount,top.count)
    mediandiff = c(mediandiff,median.diff)
  }
  tab$word = w
  tab$z.score = zscore
  tab$position = topcount
  tab$median.diff = mediandiff
  tab
}

supported.z = function(h) {
  z = (max(h$counts)-mean(h$counts))/sd(h$counts)
  p.max = h$breaks[h$counts == max(h$counts)]
  l.z = (h$counts[p.max-1]-mean(h$counts))/sd(h$counts)
  r.z = (h$counts[p.max+1]-mean(h$counts))/sd(h$counts)
  mean(c(z,l.z,r.z))
}

draw.words = function(w,con,pdfname) {
  pdf(file=paste(sep='',pdfname,'.pdf'),paper='a4r',width=8.3,height=11.7,pointsize=8)
  par(mfrow=c(3,2))
  for(i in 1:length(w)) {  
    word = w[i]
    query = paste(sep='',"SELECT id FROM word WHERE word = '",word,"'")
    res = dbSendQuery(con,query)
    id = as.numeric(fetch(res))
    query = paste(sep='','SELECT pos FROM seq_occurrence WHERE word_id = ',id)
    res = dbSendQuery(con,query)
    pos = fetch(res,n=-1)
    pos = pos$pos
    hist(pos,breaks=seq(1,100,1),col=8,main=word)
  }
  dev.off()
}


library(RMySQL)
db = 'MM_upstream_100_1_NSH_100_WL_6'
usr = 'mysql_dev'
pwd = 'riiGbs'
zcutoff = 3
stat = 'stat_MM_upstream_100_1.fa_100_6_counts.xls.csv'
pdfname = 'MM'
drv = dbDriver("MySQL")
con = dbConnect(drv,db,usr,pwd)
t = read.table(file = stat, sep = ',', header = T)
t = na.omit(t[t$z.seq > zcutoff,])
pos = analyze.words(t$word,con)
pos = na.omit(pos[pos$z.score > zcutoff,])
pos = pos[order(pos$median.diff,decreasing=T),]
draw.words(pos$word,con,pdfname)


library(RMySQL)
db = 'FC_upstream_100_1_NSH_100_WL_6'
usr = 'mysql_dev'
pwd = 'riiGbs'
zcutoff = 2.5
stat = 'stat_FC_upstream_100_1.fa_100_6_counts.xls.csv'
pdfname = 'FC'
drv = dbDriver("MySQL")
con = dbConnect(drv,db,usr,pwd)
t = read.table(file = stat, sep = ',', header = T)
t = na.omit(t[t$z.seq > zcutoff,])
pos = analyze.words(t$word,con)
pos = na.omit(pos[pos$z.score > zcutoff,])
pos = pos[order(pos$median.diff,decreasing=T),]
draw.words(pos$word,con,pdfname)


library(RMySQL)
db = 'DM_upstream_100_1_NSH_100_WL_6'
usr = 'mysql_dev'
pwd = 'riiGbs'
zcutoff = 3
stat = 'stat_DM_upstream_100_1.fa_100_6_counts.xls.csv'
pdfname = 'DM'
drv = dbDriver("MySQL")
con = dbConnect(drv,db,usr,pwd)
t = read.table(file = stat, sep = ',', header = T)
t = na.omit(t[t$z.seq > zcutoff,])
pos = analyze.words(t$word,con)
pos = na.omit(pos[pos$z.score > zcutoff,])
pos = pos[order(pos$median.diff,decreasing=T),]
draw.words(pos$word,con,pdfname)
