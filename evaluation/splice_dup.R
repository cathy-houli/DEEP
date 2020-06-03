#input format: chr start end
#output format chr start end cov
ReadFile = "test.sort.splice";

WriteFile = "test.sort.cov.splice";

splice = read.table(ReadFile, header = FALSE, sep = '\t');
dup = duplicated(splice);
dup_start = which(dup==FALSE);

Unique = splice[!dup,];
Cov = rep(0,nrow(Unique));

Cov[1:(nrow(Unique)-1)] = dup_start[2:nrow(Unique)]-dup_start[1:(nrow(Unique)-1)];
Cov[nrow(Unique)] = nrow(splice)-dup_start[nrow(Unique)]+1;

Unique_Cov = data.frame(Unique,Cov);
options(scipen=200);

write.table(Unique_Cov, file = WriteFile, sep = "\t", row.names = FALSE,col.names = FALSE, quote = FALSE);
