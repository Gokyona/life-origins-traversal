# Stage 4c-1 specification: polymer chemistry (levels E1 and E2)

**Status:** specification, 26 October 2026 (00:20 Toronto time) · **Plan:** [stage4_plan.md](stage4_plan.md), step 4c-1 · **Scope:** liquid water, carbon chemistry, Earth-like surface or sea floor

This is the first layer of the bottom-up simulator. It contains only feedstock, monomers and polymers: no catalysts, templates or membranes yet. Later layers add to it without changing it.

## 1. What the layer must answer

For each environment (a point in layer A), the layer must answer two questions:
- **E1:** Do activated monomers reach a working concentration?
- **E2:** Do chains of length ≥ L_min (default 8; 4 and 16 as sensitivity settings) persist, that is, form faster than they hydrolyse or wash out?

## 2. State

The simulation has N reaction sites (default 256). One site is a small volume: a pore, a droplet or a pool. All sites receive the same environment; they differ only by chance. Each site holds:
- counts of feedstock F, activated monomers M\*_a and inactive monomers M_a, for an alphabet of size A (default 2, later 4);
- polymers stored explicitly as sequences up to L_max = 32. Each site has a fixed number of slots (default 2,048), and a polymer is a length plus a bit-packed sequence. A polymer that reaches L_max stops growing in this layer.

Sequences are already stored so that later layers can give them functions (catalysis, templating) without changing the format.

## 3. Processes and rate laws

The time step is adaptive (tau-leap), with time measured in hours.

| process | rate | environment dependence | status |
|---|---|---|---|
| feedstock supply | s_F | feedstock supply rate (layer A); composition sets the monomer alphabet weights | input |
| monomer synthesis F → M | k_syn · F | energy source: UV flux × UV yield (layer C), or redox (H₂) × redox yield (layer C) | yields are layer C |
| activation M → M\* | k_act · M | same energy terms; dehydration during dry phases | layer C |
| deactivation M\* → M | k_deact · M\* | temperature (Arrhenius, Ea layer C) | layer C |
| chain growth (condensation) | k_c · c_eff · M\* · (chain ends) | concentration factor (evaporation × thermal gradient × mineral adsorption); dry phase multiplies by f_dry (layer C) | k_c and f_dry are layer C |
| chain initiation M\* + M\* → dimer | k_c · c_eff · M\*² | as above | |
| hydrolysis per bond | k_h(T, pH, Mg²⁺, K⁺) | Li & Breaker (1999) form, see below; dry phase multiplies by water activity a_w (≈ 0 when dry) | form verified; overall scale shared with RNA |
| washout | 1/τ for all free species | residence time τ (dilution and flow) | input |
| UV damage | k_UV · UV per bond | UV flux | layer C |
| mineral protection | adsorbed fraction θ hydrolyses at k_h · (1 − p_min) | mineral surface availability and type | θ and p_min are layer C |

**Hydrolysis.** Li & Breaker (1999) measured RNA phosphodiester cleavage and gave an empirical rate law in pH, temperature, K⁺ and Mg²⁺:
- the rate rises about tenfold per pH unit;
- it rises by a factor 10^0.07 per °C, consistent with their activation energy of 29 kcal/mol;
- it rises with Mg²⁺ to about the power 0.8.

At 37 °C, pH 7.4, 0.25 M K⁺ and 5 mM Mg²⁺ their projected rate is about 1.4 × 10⁻⁷ per minute.

We use their pH, temperature and Mg²⁺ dependences relative to that reference point. The K⁺ term is replaced by a salinity term with a layer C exponent, because the printed equation is ambiguous in its K⁺ dependence, as our own recalculation showed. Outside their measured ion ranges (Mg²⁺ 5–50 mM, K⁺ 0.03–3.16 M) the law is extrapolated, and this is flagged.

Peptide-like bonds hydrolyse far more slowly: the half-life is about 400 years at 25 °C and about 5.5 weeks at 100 °C (Radzicka & Wolfenden 1996, as tabulated by Wolfenden & Snider 2001). Peptide-like chains (alphabet from amino-acid-like feedstock) therefore use their own Arrhenius line through these two points.

**Wet–dry and hot–cold cycling.** Cycles are square waves with the period and amplitude from layer A:
- **Dry phase:** a_w → a_w,dry (layer C, near 0). Concentration is multiplied by the evaporation factor, and growth is multiplied by f_dry.
- **Wet phase:** volume is restored, so all concentrations are divided by the evaporation factor.
- **Hot–cold cycling:** the temperature of the Arrhenius terms alternates.

**Time.** Real "available time" (up to millions of years) cannot be simulated step by step. Each run therefore lasts until the chain-length distribution is stationary (at most T_max steps), and is scored after that. Available time enters only where a process is slower than the run (it caps what slow steps can achieve, through the ratio of available time to the slowest time scale). This approximation is stated in every report.

## 4. Validation before any map

1. **Reversible polymerisation without cycling or flow.** At steady state the chain-length distribution must follow the most-probable (Flory) distribution P(L) ∝ p^(L−1), where p is set by the ratio of growth to hydrolysis. Compare the simulated distribution with the analytic one over at least three values of p.
2. **Arrhenius check.** The simulated hydrolysis of a fixed pool of dimers must reproduce the 10^0.07 per °C slope and the reference rate above.
3. **Cycling.** At the same mean water activity, wet–dry cycling must lengthen chains compared with constant wet conditions. This is the qualitative result that motivates the hot-spring scenario (Damer & Deamer 2020). We do not fit a number to it.
4. **Washout.** With τ much shorter than the time needed to grow chains, E2 must fail. With no feedstock, E1 must fail.
5. **Distractors.** Layer A variables that act only in later layers (for example ion effects on membranes) must have zero importance in this layer's map.

## 5. Output and map

Per run, the layer records:
- the concentration time series of M\*;
- the final chain-length histogram;
- the fraction of sites with any chain ≥ L_min;
- E1 and E2 as yes/no outcomes.
- **Composition record** (for the composition map; see the plan, Section 2.3): monomer shares by alphabet type, chain-length histogram, sequence diversity (number of distinct sequences and Shannon entropy per length class), share of chains in adsorbed (mineral) state, and the time at which E1 and E2 were first reached. Later layers add their own fields to the same record.

The map is a COMPASS joint nested map of E1 and E2 over all of layer A plus the layer C constants of this layer (about 25 dimensions). It uses the dimension handling chosen from the 4a scaling test. The step finishes with a cost measurement per run, importance, partial dependence, and two slices: temperature × cycling and pH × Mg²⁺.

## 6. Checkpoints and resuming (all stage 4 programs)

The owner requires that a power cut loses at most a few minutes of work, and that computation continues afterwards (26 October 2026).

- **Atomic saves.** Every result file is written to a temporary file and then renamed, so a cut never leaves a half-written file. Finished work is appended as one line per completed batch (JSON Lines) and flushed to disk (fsync).
- **Checkpoints at least every 5 minutes**, and at the end of every COMPASS round. A checkpoint holds the engine state: points simulated so far and their outcomes, the used-point mask, the round number and the random-generator states.
- **Resume.** Every script takes `--resume`. On restart it reads the last checkpoint and the result lines, skips finished work, and continues. Seeds are derived from (config, repetition, round, batch), so a resumed run gives the same numbers as an uninterrupted one. A small test checks this: kill a run mid-way, resume it, and compare with an uninterrupted run.
- **One command to restart everything.** `lot-s4/resume_all.sh` restarts every unfinished job with `--resume`. Starting it automatically at boot is a change to the machine, and needs the owner's approval first.
- **Off-machine copies.** Finished result files are copied to the repository after each step, so a disk failure on the benchmark machine cannot lose published results.

## 7. Cost target

Batch the runs on the GPU like the Part I models. The target is ≤ 50 ms per run, so that a 100,000-run map takes about 1.5 hours on the B70. If it is much slower, reduce N or the slot count for screening runs (multi-fidelity), and keep full runs for points near the boundary.

## References

- Damer, B. & Deamer, D. (2020). The hot spring hypothesis for an origin of life. *Astrobiology* 20(4), 429–452. (Checked for the [condition review](condition_review.md).)
- Li, Y. & Breaker, R. R. (1999). Kinetics of RNA degradation by specific base catalysis of transesterification involving the 2′-hydroxyl group. *Journal of the American Chemical Society* 121(23), 5364–5372. We checked a copy of the paper on 26 October 2026 for the volume, issue, pages, the rate law, the reference rate and the activation energy. The first page with the full title was not in the copy; the title as given here is to be confirmed.
- Radzicka, A. & Wolfenden, R. (1996). Rates of uncatalyzed peptide bond hydrolysis in neutral solution and the transition state affinities of proteases. *Journal of the American Chemical Society* 118, 6105–6109. Values as tabulated in Wolfenden, R. & Snider, M. J. (2001), *Accounts of Chemical Research* 34(12), 938–945, p. 940, via the BioNumbers database (entry 105352), checked 26 October 2026.
