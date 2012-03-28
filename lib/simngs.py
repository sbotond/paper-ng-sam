
import      utils       as      u
import      os 

class SimNGS:
    """ simulate Illumina sequencing using simNGS """
    def __init__(self, run_file, read_length, insert_size, log, rts):
        self.fq_tmp     = 'tmp'
        self.log        = log
        self.run_file   = run_file
        self.read_length= read_length
        self.insert_size= insert_size
        self.isize      = insert_size - 2 * read_length
        self.sim_cmd    = self._build_sim_cmd()
        self.rts        = rts
        u.check_cmd('simLibrary')
        u.check_cmd('simNGS')

    def sim(self, name, seq, cov):
        """ Simulate from one sequence """
        ref     = self.rts.tempfile("ref_tmp.fas")
        fh      = open(ref,"w")
        fh.write(">%s\n%s\n" % (name, seq))
        fh.flush()
        fh.close()
        cmd = self._build_lib_cmd(cov, ref) + " | " + self.sim_cmd
        ret = os.system(cmd)
        if ret != 0:
            self.log.fatal("Failed to simulate sequencing for %s" % name)
        os.unlink(ref)

    def _build_sim_cmd(self):
        cmd = "simNGS -n %d -p paired -o fastq -O %s %s 2>/dev/null; cat tmp_end1.fq >> end1.fq; cat tmp_end2.fq >> end2.fq; rm tmp_end?.fq" % (self.read_length, self.fq_tmp, self.run_file)
        return cmd
   
    def _build_lib_cmd(self, cov, ref):
        cmd = "cd %s; simLibrary -r %d -i %d -x %d %s 2>/dev/null" % (self.rts.base, self.read_length, self.isize, cov, ref)
        return cmd

