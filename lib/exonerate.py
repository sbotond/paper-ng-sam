import      utils       as      u
import      os 
from        Bio         import  SeqIO
import      itertools   as      it

class Exonerate:
    """ Perform pairwise alignments using exonerate  """
    def __init__(self, log, rts):
        self.log        = log
        self.rts        = rts
        u.check_cmd('exonerate')

    def _prepare_input(self, s1, s2):
        """ Save sequences as temporary files """
        rts = self.rts.subdir('exonerate')

        f1  = rts.tempfile('s1.fas')
        fh1 = open(f1,'w')
        fh1.write(">s1\n%s\n" % s1)
        fh1.flush()
        fh1.close()

        f2  = rts.tempfile('s2.fas')
        fh2 = open(f2,'w')
        fh2.write(">s2\n%s\n" % s2)
        fh2.flush()
        fh2.close()

        r   = rts.tempfile('out.txt')
        return f1, f2, r


    def check_strands(self, s1, s2):
        """ Check whether two sequences align with the same strand """
        if len(s1) == 0 or len(s2) == 0:
            return None
        f1, f2, r   = self._prepare_input(s1, s2)

        cmd  = """exonerate --verbose 0 --showalignment no --showvulgar no --ryo '%qS|%tS\\n' """
        cmd += "-m affine:local -e -100 -o -100 -n 1 --target %s --query %s > %s" % (f1, f2, r)

        if os.system(cmd) != 0:
            return None

        ryo = file(r).readlines()
        if len(ryo) != 1:
            return None

        ryo = ryo[0].rstrip()
        if ryo == '':
            return None

        # Parse strands:
        strands = ryo.split('|')
        if strands[0] == strands[1]:
            return False
        else:
            return True
        
    def seq_cmp(self, s1, s2):
        """ Compare two sequences by pairwise alignment """
        f1, f2, r   = self._prepare_input(s1, s2)

        cmd  = "exonerate --verbose 0 --showalignment no --showvulgar no -m affine:local -e -100 -o -100 --ryo '%tab|%tae|%pi\\n' -n 1 "
        cmd += "--target %s --query %s > %s " % (f1, f2, r)

        if os.system(cmd) != 0:
            return None

        ryo = file(r).readlines()
        if len(ryo) != 1:
            return None

        ryo = ryo[0].rstrip()
        if ryo == '':
            return None

        qb, qe, pi  = ryo.split('|')
        qb  =   int(qb)
        qe  =   int(qe)
        pi  =   float(pi)

        aln_len = abs(qe - qb)

        res = {
            'percent_identity':         pi,
            'alignment_length':         aln_len,
        }
        return res
