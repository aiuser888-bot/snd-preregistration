"""SND reference implementation (protocol V2.0.13). See ../../README.md."""
__version__ = "2.0.13"
from .design import make_block_rng, block_angles
from .stats import center_two_way, block_scores, lam
from .randtest import randomization_test
from .fcrice import build_table, load_table, ucb95, TABLE_CHECKPOINTS
from .kappa import rho_x_bernoulli_copula, kappa_exact, coverage_prediction
