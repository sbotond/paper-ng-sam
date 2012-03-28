import      utils       as      u
import      os 
from        Bio         import  SeqIO
import      itertools   as      it

class Velvet:
    """ Assemble reads using velvet """
    def __init__(self, fqs, kmer_length, min_ctgl, ins_len, max_div, rts, log, clean=False, asm_dir="assembly"):
        self.log        = log
        self.rts        = rts
        self.fqs        = fqs
        self.kmer_length= kmer_length
        self.min_ctgl   = min_ctgl
        self.ins_len    = ins_len
        self.max_div    = max_div
        self.clean      = clean
        self.asm_rts    = self.rts.subdir(asm_dir)
        self.asm_dir    = self.asm_rts.base
        u.check_cmd('velveth')
        u.check_cmd('velvetg')

    def velveth(self):
        """ Build k-mer hash """
        self._prepare_fasta()
        cmd = "velveth %s %d -shortPaired -fasta %s >/dev/null" % (self.asm_dir, self.kmer_length, self.fasta)
        if os.system(cmd) != 0:
            return None
        for f in ('Log','Roadmaps', 'Sequences'):
            self.asm_rts.register(os.path.join(self.asm_dir, f))
        return True

    def _prepare_fasta(self):
        """ Prepare input for velveth """
        output   = self.rts.tempfile("reads.fas")
        ofh      = open(output, 'w')
        stream1  = SeqIO.parse(self.fqs[0],'fastq')
        stream2  = SeqIO.parse(self.fqs[1],'fastq')

        for record in it.chain( it.izip(stream1, stream2) ):
            SeqIO.write(record, ofh, 'fasta')    

        ofh.flush()
        ofh.close()
        self.fasta  = output

    def velvetg(self):
        """ Assemble short reads using velvetg """
        cmd  = "velvetg %s -scaffolding no -read_trkg yes -max_divergence %f "  % (self.asm_dir, self.max_div)
        cmd += "-min_contig_lgth %d -exp_cov auto -ins_length %d >/dev/null"    % (self.min_ctgl, self.ins_len)
        if os.system(cmd) != 0:
            return None
        for f in ('contigs.fa','Graph2', 'LastGraph', 'PreGraph', 'stats.txt'):
            self.asm_rts.register(os.path.join(self.asm_dir, f))
        return True

    def parse_contigs(self):
        """ Parse contigs """
        cf  = os.path.join(self.asm_dir, 'contigs.fa')
        if not os.path.exists(cf):
            return None
        contigs = u.Fasta(cf).slurp()
        if len(contigs) == 0:
            return None
        return contigs

