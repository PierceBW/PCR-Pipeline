# Virtual environment (venv)

This project uses a **venv** (no conda). Use it so scripts like `species_filter_fast.py` have `beautifulsoup4`, `tqdm`, and `biopython` available.

## Setup (already done once)

```bash
cd /path/to/PCR-Pipeline
python3 -m venv venv
source venv/bin/activate    # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Use it

**Option A – Activate, then run anything**

```bash
source venv/bin/activate
python species_filter_fast.py --cutoff 95
# or
jupyter notebook
```

**Option B – Run without activating**

```bash
./venv/bin/python species_filter_fast.py --cutoff 95
./venv/bin/jupyter notebook
```

## What’s installed

- **beautifulsoup4** – parse HTML reference files
- **tqdm** – progress bars
- **biopython** – sequence/tree handling in notebooks and `reroot_trees.py`
