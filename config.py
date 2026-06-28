# ─────────────────────────────────────────────
# PhyloGenie Configuration
# ─────────────────────────────────────────────
# 1. Sign up free at https://www.ncbi.nlm.nih.gov/account/
# 2. Get API key at https://www.ncbi.nlm.nih.gov/account/settings/
# 3. Paste your email and API key below

NCBI_EMAIL = "f20230843@pilani.bits-pilani.ac.in"   # REQUIRED — put your email here
NCBI_API_KEY = "5b385ee9fa60f950abefddcae7462963b109"                        # OPTIONAL but recommended (10x faster)

# Output directory for job files and PDFs
OUTPUT_DIR = "outputs"

# Max sequences per job
MAX_SEQUENCES = 10

# BLAST settings
BLAST_HITLIST_SIZE = 20
BLAST_EXPECT = 0.001

# Augustus web API
AUGUSTUS_URL = "https://bioinf.uni-greifswald.de/augustus/submission.php"
