library(RMySQL)
library(beeswarm)
org.code = 'DR';
word = 'TATAAT'
wl = 6
spacing = 1
db = paste(sep='',org.code,'_upstream_100_1_NSH_100_WL_',wl)
usr = 'mysql_dev'
pwd = 'riiGbs'
pdfname = paste(sep='',org.code,'_',word,'.pdf')
drv = dbDriver("MySQL")
con = dbConnect(drv,db,usr,pwd)

query = paste(sep='',"SELECT id FROM word WHERE word = '",word,"'")
res = dbSendQuery(con,query)
id = as.numeric(fetch(res))

query = paste(sep='','SELECT COUNT(DISTINCT seq_id) FROM seq_occurrence WHERE word_id = ',id)
res = dbSendQuery(con,query)
seq = fetch(res,n=-1)
seq = seq[,1]

query = paste(sep='','SELECT sequences FROM shuffled_counts WHERE word_id = ',id)
res = dbSendQuery(con,query)
shu = fetch(res,n=-1)
shu = shu[,1]

query = paste(sep='','SELECT pos FROM seq_occurrence WHERE word_id = ',id)
res = dbSendQuery(con,query)
pos = fetch(res,n=-1)
pos = pos$pos

pdf(file=pdfname,paper='a4r',width=9,height=12,pointsize=8)
par(mfrow=c(1,2))

hist(pos-100,breaks=50,col=8,main=paste(sep=' ',org.code,'-',word,'positions in core promoters'),xlab='bp relative to putative TSS',ylab='number of promoters')

offset = 150
min = min(c(seq,min(shu))) - offset
max = max(c(seq,max(shu))) + offset
beeswarm(c(shu,seq),method='center',pwcol=c(rep(2,length(shu)),1),pch=19,horiz=F,ylim=c(min,max),main=paste(sep=' ',org.code,'-','promoters containing',word),ylab='number of promoters',labels='',spacing=spacing)
legend('topright',legend=c('shuffled','real'),pch=19,cex=1.5,col=c(2,1))

dev.off()

