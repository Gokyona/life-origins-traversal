# Stage 4c-2 specification: catalysis and self-sustaining networks (level L1)

**Status:** draft specification, 26 October 2026 (01:15 Toronto time) · **Plan:** [stage4_plan.md](stage4_plan.md), step 4c-2 · **Builds on:** [4c-1](stage4c1_spec.md) (polymer chemistry), [stage 1e](../results/stage1e_raf_benchmark_b70.md) (autocatalytic sets)

The second layer adds catalysis to the polymer chemistry of 4c-1. Nothing else changes. The question is whether, and in which environments, a set of polymers arises that catalyses its own formation from the available food: a self-sustaining chemical network. This is level L1, "self-maintenance".

## 1. Catalysis rule

- **Who catalyses.** Any polymer of length ≥ L_cat (layer C, 3–8) may catalyse ligation reactions A + B → AB between two other polymers or monomers.
- **Which reactions.** Whether sequence C catalyses reaction (A, B) is decided by a fixed hash of the three sequences: catalysis occurs if hash(C, A, B) < p_cat(L_C). The assignment is therefore random but reproducible: the same molecules always have the same chemistry, in every site and every run. This is the polymer model of Kauffman and of Hordijk & Steel, applied to the molecules that actually form.
- **How likely.** p_cat(L) = p₀ · (1 − e^(−(L − L_cat + 1)/λ)), so longer chains are more likely to catalyse, up to p₀. Both p₀ and λ are layer C dimensions over wide logarithmic ranges, because the catalytic probability of prebiotic polymers has not been measured in this form.
- **How much.** A catalysed ligation runs κ times faster than the uncatalysed one (layer C, 10–10⁶), in proportion to the catalyst's concentration.
- **Minerals.** Mineral surfaces add generic, non-specific catalysis of condensation: a factor scaled by mineral availability and type (layer C per type). This lets the map show whether mineral catalysis can replace or seed polymer catalysis.
- **Environment.** Catalysed rates follow the same temperature dependence as uncatalysed ones, with their own activation energy (layer C, lower than uncatalysed). Catalysts are hydrolysed like any polymer, which links catalysis back to pH, temperature and Mg²⁺.

## 2. Detector for L1

At checkpoints during the run, and at the end:

1. **Build the network.** Collect the realised reaction network in each site: species present above a threshold count, and the ligation reactions among them with their catalysts.
2. **Find autocatalytic sets.** Run the RAF algorithm of Hordijk & Steel on that network. The food set is the monomers and the feedstock.
3. **Score L1.** L1 holds if all three conditions are met:
   - **Size:** a RAF set with at least three reactions exists.
   - **Persistence:** its members persist over at least five molecular turnover times (the mean lifetime of a member against hydrolysis and washout).
   - **Self-sustaining:** the members are more abundant than in the same run with catalysis switched off (a paired control run with the same seed). This rules out sets that exist only on paper.

Composition record fields added by this layer:
- RAF size and number of independent RAF subsets;
- share of mineral against polymer catalysis;
- catalyst lengths;
- monomer composition of RAF members;
- the time L1 was first reached, relative to E2.

## 3. Validation before any map

1. **Known threshold.** In a well-mixed site with constant concentrations, set the catalysis parameters directly. The probability that a RAF exists must rise steeply at about one to two reactions catalysed per molecule. This is the threshold reported by Hordijk, Kauffman & Steel (2011) and reproduced in stage 1e within 2%.
2. **Catalysis off.** With p₀ = 0, L1 must never occur, including with mineral catalysis alone. Mineral catalysis is non-specific, so it cannot form a RAF set by itself; it only speeds up E2.
3. **Shuffled control.** Re-hash catalysis each step instead of fixing it. Persistent RAF sets must disappear, because a set that is no longer the same chemistry cannot sustain itself.
4. **Hydrolysis link.** At the reference catalysis settings, raising the temperature or pH must shorten catalyst lifetimes and shrink the L1 region in a direction consistent with the hydrolysis law of 4c-1.

## 4. Map

A joint nested map of E1, E2 and L1 over layer A and the layer C constants of 4c-1 and 4c-2 (about 30 dimensions). It uses COMPASS as validated at that dimension in [4a](../results/stage4a_scaling_b70.md). Because L1 requires E2, simulations stop early where E2 fails, which saves compute in most of the space. Checkpoints and resuming are as in 4c-1, Section 6.

## 5. What this layer cannot show

Catalysis is assigned at random, not derived from chemistry. The map therefore says, for each environment, how much catalytic ability the polymers would need for a self-sustaining network to arise. It does not say whether real prebiotic polymers had it. That requirement becomes a layer C dimension. Its importance, and the threshold the map finds for it, are themselves the result, and they can be compared with experiments that measure the catalytic activity of random sequences.

## References

- Hordijk, W., Kauffman, S. A. & Steel, M. (2011). Required levels of catalysis for emergence of autocatalytic sets in models of chemical reaction systems. *International Journal of Molecular Sciences* 12(5), 3085–3101. (Checked for the stage 1 report.)
- Hordijk, W. & Steel, M. (2004). Detecting autocatalytic, self-sustaining sets in chemical reaction systems. *Journal of Theoretical Biology* 227, 451–461. (Checked for the stage 1 report and the homepage reference list.)
- Kauffman, S. A. (1986). Autocatalytic sets of proteins. *Journal of Theoretical Biology* 119, 1–24. (Checked for the homepage reference list.)
