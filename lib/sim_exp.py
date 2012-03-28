import      utils       as      u
import      os 
import      dendropy

class SimPcrDil:
    """ Simulate PCR amplifications and dilutions """
    def __init__(self, name, init_popsize, pcr_eff, nr_cycles_mut, dilf_after_mut, nr_cycles_cln, dilf_after_cln, nr_cycles_cov, total_cov, sample_size_mut, path, log, rdir='.', clean=True, mut_only=False):
        self.log            = log
        self.name           = name
        self.init_popsize   = init_popsize
        self.pcr_eff        = pcr_eff
        self.nr_cycles_mut  = nr_cycles_mut
        self.dilf_after_mut = dilf_after_mut
        self.nr_cycles_cln  = nr_cycles_cln
        self.dilf_after_cln = dilf_after_cln
        self.nr_cycles_cov  = nr_cycles_cov
        self.total_cov      = total_cov
        self.sample_size_mut= sample_size_mut
        self.path           = path
        self.script         = os.path.join(self.path,'pcr_coal.R')
        self.rdir           = rdir
        self.clean          = clean
        self.mut_only       = 0
        if mut_only:
            self.mut_only   = 1
        self.out_nwk        = os.path.join(self.rdir, self.name + ".nwk")
        self.out_cov        = os.path.join(self.rdir, self.name + ".cov")
        self._check_script()

    def simulate(self):
        """ Run NG-SAM simulation """
        ret     = 256
        count   = 0 
        while True:
            if ret == 0:
                break
            elif count < 6:
                self._cleanup()
                cmd = self._construct_cmd()
                ret = self._run_cmd(cmd)
                self.sample_size_mut *= 2
                count   += 1
            else:
                self.log.fatal("Failed to simulate experiment!")

        res = self._harness_results()
        if self.clean:
            self._cleanup()
        return res
           
    def _cleanup(self):
        """ Remove pcr_coal.R output """
        if os.path.exists(self.out_nwk):
            os.remove(self.out_nwk)
        if os.path.exists(self.out_cov):
            os.remove(self.out_cov)

    def _harness_results(self):
        res = {}
        cov = {}
        for line in file(self.out_cov):
            name, cv    = line.split()
            cov[name]   = int(float(cv) * self.total_cov)
        res['cov']  = cov
        if len(res['cov']) > 0:
            res['tree'] = dendropy.Tree.get_from_path(self.out_nwk, schema='newick', as_rooted=True)
        return res

    def _check_script(self):
        if not os.path.exists(self.script):
            self.log.fatal('Simulation script is not in sepcified path!')

    def _construct_cmd(self):
        args    = [
            self.name,
            self.init_popsize,
            self.pcr_eff,
            self.nr_cycles_mut,
            self.dilf_after_mut,
            self.nr_cycles_cln,
            self.dilf_after_cln,
            self.nr_cycles_cov,
            self.sample_size_mut,
            self.mut_only,
        ]
        cmd = self.script
        for arg in args:
            cmd += " %s" %arg
        cmd += " 2>/dev/null"
        return cmd

    def _run_cmd(self, cmd):
        old_wd  = os.getcwd()
        os.chdir(self.rdir)
        ret = os.system(cmd)
        os.chdir(old_wd)
        return ret
