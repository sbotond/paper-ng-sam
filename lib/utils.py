import      os
import      sys
import      time
import      copy
import      subprocess                      as      sp
from        Bio                             import  SeqIO
from        Bio                             import  Seq
from        Bio.Alphabet                    import  generic_dna
import      matplotlib
matplotlib.use('Agg')
from        matplotlib                      import  pyplot          as  plt
from        collections                     import  defaultdict
from        matplotlib.backends.backend_pdf import  PdfPages
import      numpy                           as      np

class Log:
    """ Logging utility class """
    def __init__(self, fname=None, level=0):
        self.level = level
        if fname == None:
            self.fname  = "<sys.stderr>"     
            self.file   = sys.stderr
        else:
            self.file   = open(fname, "w")
            self.fname  = fname

    def close(self):
        """ Close log file """
        self.file.flush()
        self.file.close()

    def log(self, message):
        """ Log message """
        if self.level < 0:
            return
        self.file.write("[%s] %s\n" % (time.strftime("%y-%m-%d %H:%M:%s"), message) )

    def fatal(self, message):
        """ Log message and exit """
        self.file.write("[%s] %s\n" % (time.strftime("%y-%m-%d %H:%M:%s"), message) )
        sys.exit(1)

class Rtemp:
    """ Utility class handling temporary storage """
    def __init__(self, base, log, autoclean=False):
        self.log        = log
        self.autoclean  = autoclean
        self.parent     = None
        self.children   = [ ]
        self.files      = [ ]
        if os.path.isdir(base) != True:
            log.fatal("The base must be a directory: %s" % base)
        self.base   = os.path.abspath(base)

    def exists(self, fname):
        """ Check wheteher a file exists """
        return os.path.exists(fname) 

    def _iterate_fname(self, fname):
        """ Iterate until we don't have a name clash """
        i           = 0
        orig_fn     = fname # basename
        while True:
            if (fname in self.files) or self.exists(fname):
                i       += 1
                fname   = orig_fn + ("_%03d" % i)
            else:
                break
        return fname

    def tempfile(self, name):
        """ Get a temporary file """
        fname       = self._iterate_fname(os.path.join(self.base, name))
        self.register(fname)            
        return fname

    def temp_fh(self, name):
        """ Get a temporary file handle """
        fname   = self.tempfile(name)
        f       = open(fname, "w")
        return f

    def clean(self):
        """ Remove registered temporary files """
        for child in self.children:
            child.clean()   # call cleanup on children.
        tmp = list(self.files)
        for f in tmp:
            self.remove(f)
        # Delete the directory if children is  a subdir:
        if self.parent != None:
            os.rmdir(self.base)

    def remove(self, fname):
        """ Remove a temporary file """
        if not (fname in self.files):
            self.log.fatal("The file %s is not mannaged by this object!" % fname)
        if os.path.exists(fname):
            os.remove(fname)
        self.files.remove(fname)

    def subdir(self, dname):
        """ Get a mannaged temporary subdirectory """
        clone           = copy.copy(self)        
        dname           = os.path.join(self.base, dname)
        clone_base      = self._iterate_fname(dname)
        os.mkdir(clone_base)
        clone.base      = clone_base
        clone.parent    = self
        clone.children  = []
        clone.files     = []
        self.children.append(clone)
        return clone

    def register(self, fname):
        """ Register temporary file """
        self.files.append(fname)
    
    def unregister(self, fname):
        """ Unregister temporary file """
        self.files.remove(fname)

    def __del__(self):
        if self.autoclean:
            self.clean()

class Fasta:
    """ Fasta parsing class """
    def __init__(self, infile):
        self.infile     = infile
        self.in_fh      = open(infile, "r")
        self.iter       = SeqIO.parse(self.in_fh,'fasta')

    def __iter__(self):
        """ Return iterator """
        return iter(self.iter)

    def slurp(self):
        """ Slurp sequences """
        records = { }
        for s in iter(self):
            records[s.name] = str(s.seq)
        return records

def parse_target_seq(f):
    """ Read in and parse target sequence """
    tmp = Fasta(f).slurp()
    res = { }
    if len(tmp) != 1:
        raise ValueError("Too many sequences in file %s!" % f)
    name    = tmp.keys()[0]
    seq     = tmp.values()[0]
    res['full_name'] = name
    res['seq']       = seq.upper()
    return res

def parse_bl_file(blf):
    """ Read in branch length scaler file """
    return float(file(blf).readlines()[0])

class Res:
    """ Class for saving results """
    def __init__(self, f, name):
        self.f      = f
        self.name   = name

    def save(self,status,nmut=-10, targ_len=-10, cons_len=-10, seq_ident=-10):
        if status not in (0,-1,-2,-3,-4,-5,-6, -7):
            raise ValueError("Invalid result status: %s" % status)
        fh  = open(self.f,"w")
        fh.write("%s\t%d\t%d\t%d\t%d\t%6f\n" %(self.name,status, nmut, targ_len, cons_len, seq_ident) )
        fh.flush()
        fh.close()


def check_cmd(name):
    """ Check whether a command is in the path """
    cmd = "which %s > /dev/null" % name
    if os.system(cmd) != 0:
        raise ValueError("The command \"%s\" is not present in the path!")

def revcomp(seq):
    """ Reverse complement sequence using Biopython """
    tmp = Seq.Seq(seq, generic_dna) 
    tmp = tmp.reverse_complement()
    return str(tmp)

def write_fasta(sequences, fname):
    """ Write sequence to fasta file """
    fh  = open(fname, 'w')
    for name, seq in sequences.iteritems():
        fh.write(">%s\n%s\n" % (name, seq))
    fh.flush()
    fh.close()

def consensus(aln):
    """ Calculate majority-rule consensus """
    aln_len = aln.get_alignment_length()
    cons    = ''
    for i in xrange(aln_len):
        counts  = defaultdict(int)
        col = aln[:, i]
        for base in col:
            counts[base] += 1
        swapped = dict (zip(counts.values(), counts.keys()))
        char    = swapped[max(swapped.keys())] 
        if char != '-':
            cons += char
    return cons 

def gen_target(ulen, unr):
    """ Generate a target sequence with unit number and length specified """
    # Generate random unit:  
    alphabet    = ('A','T','G','C')
    unit        = ''
    for i in xrange(ulen):
        unit += alphabet[ np.random.randint(len(alphabet)) ]
    # Return sequence:
    return unit * unr

class Report:
    """ Class for plotting reports """
    def __init__(self, pdf):
        self.pdf    = pdf
        self.pages  = PdfPages(pdf)

    def hexbin(self, x, y, z,title="", xlab="", ylab="",f=np.mean, grids=100, xs='linear',ys='linear', cmap=None,tfs=11, point=None,cline=None):
        """ Simple hexbin plot """
        fig = plt.figure()
        plt.hexbin(x=x, y=y, C=z, reduce_C_function=f, gridsize=grids, xscale=xs, yscale=ys, cmap=cmap)
        xv  = [ min(x), max(x) ]
        yv  = [ min(y), max(y) ]
        if point != None:
            plt.plot(
                point['x'],
                point['y'],
                color=point['color'],
                marker=point['marker'],
                ms=point['ms'],
                hold=True,
            ) 
            xv.append(point['x'])
            yv.append(point['y'])
        if cline != None:
           self.plot_cline(plt, cline, x, y) 
        plt.xlim( (min(xv),max(xv)) )
        plt.ylim( (min(yv), max(yv)) )
        plt.colorbar()
        plt.xlabel(xlab, fontsize=tfs)
        plt.ylabel(ylab, fontsize=tfs)
        plt.title(title, fontsize=tfs)
        self.pages.savefig(fig)
        plt.close(fig)

    def plot_cline(self, plt, cline, x, y):
        step = 0.125
        x = np.arange(min(x), max(x),step,dtype=float)
        y = np.arange(min(y), max(y),step,dtype=float)
        xx = []
        yy = []
        for a in x:
            for b in y:
                if a * b == float(cline):
                    xx.append(a)
                    yy.append(b)
        plt.plot(xx,yy,'--',lw=1,color='black')

    def plot_hash(self, h, title="", xlab="", ylab=""):
        """ Visualise hash as a bar plot """
        fig = plt.figure()
        plt.bar(h.keys(), h.values(), width=0.1)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)

    def close(self):
        self.pages.close()

    def __del__(self):
        self.close()
