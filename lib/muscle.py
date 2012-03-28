import      utils       as      u
import      os 
from        Bio         import  AlignIO
import      itertools   as      it

class Muscle:
    """ Perform progressive multiple alignment using muscle """
    def __init__(self, log, rts):
        u.check_cmd('muscle')
        self.log        = log
        self.rts        = rts.subdir('muscle')

    def align_contigs(self, contigs):
        """ Calculate multiple alignment of contigs """
        infile   =   self.rts.tempfile('contigs.fas')
        outfile  =   self.rts.tempfile('contigs_aln.fas')
        u.write_fasta(contigs, infile)

        # Run muscle:
        cmd = "muscle -in %s -out %s -maxiters 1 -diags 2>/dev/null" % (infile, outfile)
        if os.system(cmd) != 0:
            return None

        # Read in alignment:
        aln = AlignIO.read(open(outfile), "fasta") 
        return aln
