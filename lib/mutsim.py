import      utils       as      u
import      os 
import      dendropy
import      numpy       as      np

class MutSim:
    """ Simulate mutations along a tree using a continous-time Markov process """
    def __init__(self, model_file, bl_scaler=1.0, alphabet=('A', 'T','G','C')):
        self.model_file = model_file
        self.alphabet   = alphabet
        self.epsilon    = 10.0**-10 # A small value.
        self.model      = self.read_model(model_file)
        self.bl_scaler  = bl_scaler
        u.check_cmd('R')

    def read_model(self, model_file):
        """ Read mutation model from file """
        alphabet_length = len(self.alphabet)
        # Construct unsclaed Q matrix: 
        m   = np.matrix( np.zeros((alphabet_length, alphabet_length), dtype=float) )
        for l in file(self.model_file):
            name, rate  = l.split()
            rate        = float(rate)
            # Set zero rates to a small value:
            if rate == 0.0:
                rate = self.epsilon
            f, t        = name.split("->")
            m[ self._idx(f,t) ]   = rate
        # Build scaled Q matrix:
        self.m          = self._build_model(m)
        # Precalculate eigen-decomposition:
        self.u, self.v  = np.linalg.eig(self.m)
        self.v_inv      = np.linalg.inv(self.v)

    def _check_root_seq(self, root_seq):
        """ Check the validity of the root sequence """
        for s in root_seq:
            if s not in self.alphabet:
                raise ValueError("Root sequence symbol %s not in alphabet!" % s)

    def _idx(self,f, t):
        """ Get tuple idices for a pair of bases """
        for s in (f, t):
            if s not in self.alphabet:
               raise ValueError("Symbol %f not in alphabet!" % s)
        return (self.alphabet.index(f), self.alphabet.index(t))

    def _build_model(self, m):
        for i in xrange(m.shape[1]):
            m[i,i]  = -np.sum(m[i,])
        m = self._scale_q(m)
        return m

    def _scale_q(self, q):
        """ Scale Q matrix """
        self.eq_dist    = self._calc_eq_dist(q)
        # Calculate expected number of substitutions at equilibrium:
        exp_subst       = 0
        for i in xrange(len(self.alphabet)):
            for j in xrange(len(self.alphabet)):
                if i == j: continue
                exp_subst += self.eq_dist[i] * q[i,j]
        # Calculate scaling factor and return scaled
        # Q matrix:
        self.scaling_factor = 1.0/exp_subst
        return q * self.scaling_factor

    def calc_exp_subst(self, bf):
        """ Calculate the expected number of substitutions at equilibrium """
        exp_subst       = 0
        q               = self.m
        for i in xrange(len(self.alphabet)):
            for j in xrange(len(self.alphabet)):
                if i == j: continue
                exp_subst += bf[i] * q[i,j]
        return exp_subst


    def _calc_eq_dist(self, q):
        """ Get equilibrium distribution using eigen-decomposition """
        epsilon         = np.finfo(np.double).eps * 10**2
        u,v             = np.linalg.eig(q.T) 
        i               = np.where(np.abs(u) < epsilon)[0]
        eq_dist         = v[:,i]
        eq_dist         = eq_dist/np.sum(eq_dist)
        eq_dist         = np.squeeze(np.asarray(eq_dist))
        return eq_dist

    def calc_P(self, t):
        """ Calculate transition probabilities for time t """
        u, v, v_inv  = self.u, self.v, self.v_inv
        U    = np.diag(np.exp(u * t))
        P    = np.dot( np.dot(v, U), v_inv)
        return P

    def sim(self, tree, root_seq):
        """ Simulate substitutions along a tree """
        if not tree.is_rooted:
            raise ValueError("Cannot simulate on unrooted tree!")
        self._check_root_seq(root_seq)
        self.root_seq   = root_seq
        self.sequences  = { }
        count           = 0

        for edge in tree.preorder_edge_iter():
            # Discarding edges with illegal lengths:
            if edge.length is None:
                continue
            tail_seq    =   self._get_seq(edge.tail_node, count)
            head_seq    =   self._evolve_branch(seq=tail_seq, length=(edge.length * self.bl_scaler))
            self.sequences[edge.head_node]  = head_seq

    def get_tips(self):
        """ Get the tip labels from a dendropy object """
        res = { }
        for node, seq in self.sequences.iteritems():
            if node.is_leaf():
               res[node.taxon.label] = seq 
        return res

    def _get_seq(self, node, count):
        """ Get the sequence associated with a node """
        if count == 0:
           self.sequences[node] =   self.root_seq
        return self.sequences[node]

    def _evolve_branch(self, seq, length):
        """ Simulate substitutions along a branch """
        if length == 0:
            return seq
        P       = self.calc_P(length) 
        sites   = list(seq) 
        for i in xrange(len(sites)):
            symbols     = self.alphabet
            probs       = np.array(P[self.alphabet.index(sites[i]),])
            probs.shape = 4
            new_idx     = np.random.multinomial(1, probs, 1)[0]
            new_idx     = np.where(new_idx == 1)[0][0]
            sites[i]    = self.alphabet[ new_idx ]
        return ''.join(sites) 

def calc_basefreq(seq):
    """ Calculate base frequencies of a sequence """ 
    alphabet    = ('A','T','G','C')
    res = np.zeros((len(alphabet)), dtype=float)
    for s in seq:
        if s not in alphabet:
            raise ValueError("calc_basefreq: symbol %s not in alphabet!" % s)
        i   = alphabet.index(s)
        res[i] += 1
    return res/np.sum(res)

def hm_dist(s1, s2):
    """ Calculate the hamming distance for a pair of sequences """
    if len(s1) != len(s2):
        raise ValueError('Size mismatch')
    a1  = np.asarray(list(s1))
    a2  = np.asarray(list(s2))
    match   = float(np.sum(a1 == a2))
    l       = float(len(a1))
    return (l-match)/l

