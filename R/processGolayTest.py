from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import *
from golay import *
import commands
import sys


print str(sys.argv[1])
print str(sys.argv[2])
print str(sys.argv[3])


out_handle = open(str(sys.argv[2]),"w")

for seq_record in SeqIO.parse(str(sys.argv[1]), "fastq"):  
  
  #print str(seq_record.seq)
  res = decode(str(seq_record.seq))
  print str(seq_record.seq)
  print res

  if res[0] != None:
    SeqIO.write(SeqRecord(Seq(res[0], SingleLetterAlphabet()), id=seq_record.id, description=""), out_handle, "fasta")

done_handle = open(str(sys.argv[3]),"w")
done_handle.close()
