
# 
# $1=/home/ubuntu/SHARE/ISseqOutput/INSPIIREDData/I1_fq_matched.fq
# $2=/home/ubuntu/SHARE/ISseqOutput/INSPIIREDData/R1_fq_matched.fq
# $3=/home/ubuntu/SHARE/ISseqOutput/INSPIIREDData/R2.fq


paste -d '\n' $1 $2 | sed -n 'p;n;n;N;s/\n//p' > $3
gzip $3
