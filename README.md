# Direct-fit SVI and restricted W-SVI

Research discussion material on weighted algebraic-conic fitting and global butterfly-arbitrage checks for standard SVI and a restricted quartic W-SVI family.

## Research note

Open [Direct-fit WSVI discussion note](Direct_fit_WSVI_discussion_note.html). This is a self-contained HTML document with embedded equations and comparison figures. Download the file to view it in a browser.

## Reproduce the examples

```bash
python -m pip install -r requirements.txt
python svi_wsvi_nobfly.py --model wsvi --outdir results/wsvi
python svi_wsvi_nobfly.py --model svi --outdir results/svi
```

The script produces smile, Durrleman and combined comparison PNG plots, and supports ATM weighting with `--boost` and `--width`.

The W-SVI model here is a **restricted quartic** with fixed effective scale, not the complete Clevenger–Wan WSVI family. The optimizer is local. Polynomial tests certify butterfly-arbitrage feasibility for rounded parameters, not global optimality. Numerical examples use synthetic targets; they are not empirical market fits.

## Contents

- `Direct_fit_WSVI_discussion_note.html` — self-contained explanation and figures
- `svi_wsvi_nobfly.py` — executable examples in both model families
- `verify_vogt_algebra.py` — symbolic verification utility
- `figures/` — representative comparison plots

## References

- Clevenger and Wan (2026), *WSVI: A Dimensionless Shape Family for Implied Volatility and Its Static No-Arbitrage Structure*, https://arxiv.org/abs/2608.22620
- Gatheral and Jacquier (2014), *Arbitrage-free SVI volatility surfaces*
