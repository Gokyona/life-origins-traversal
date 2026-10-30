# About Gokyona Environment

**Version 1.0 · 2026 · Author: Gokyona Wahiro (pen name), independent researcher, Toronto**

- GitHub: <https://github.com/Gokyona/life-origins-traversal>
- Homepage: <https://gokyona.github.io>
- License: code MIT · text CC BY 4.0

### What it is

Gokyona Environment is an origin-of-life simulation platform that runs on your own computer. Starting from simple feedstock, it simulates a whole chemical chain:

environment → feedstock and activated monomers → chains (polymerisation, breakage) → template copying (with copy errors) → amphiphiles self-assembling into vesicles (enclosing, growing, dividing, washing out) → optional selection on a sequence "trait".

The page shows the run live, with the simulated time and a Darwin panel (real half vs control half, lineage and clade spread, ρ(w,z), the L4 signal and more).

The platform's engine is the Stage-4 engine of paper 2: with one 64-bit word per chain it is bit-identical to the paper's code. The platform extends it to long sequences.

### If you use this platform, please cite

1. **COMPASS: mapping where life-like behaviour appears in stochastic models, with an ordered, non-repeating active search.** Gokyona Wahiro, 2026 (draft; preprint link to follow).
   COMPASS is an ordered, never-repeating active-learning engine that maps where levels of life-like behaviour appear in the condition space of stochastic models. On six models of up to 29 dimensions it needed 1.6–8× fewer simulations than Sobol or random sampling for the same map accuracy; on a 54-dimensional model it found environments a ground-truth sample never hit.

2. **How far up the ladder of life can environments carry a bottom-up chemistry? An ordered search over 13 environmental variables.** Gokyona Wahiro, 2026 (draft; preprint link to follow).
   A fully bottom-up artificial chemistry, searched layer by layer over 13 environmental variables. It climbs to coupled protocells only in cool, slowly flushed, enclosed water bodies with temperature cycles of days to weeks, and the ladder stops before Darwinian evolution. The engine of this platform is the Stage-4 engine of this paper.

### Limitations (please read)

- This is a model, not life.
- Heredity relies on very fast copying. In the model, template-copy extension is about 10⁹–10¹¹× faster than measured non-enzymatic chemistry. At measured rates, heredity in the model cannot start from scratch.
- The selection "trait" is an imposed score table, not derived from the chemistry.
- No single lineage took over from the others in any test.
- The time scale follows the model's rate constants; it is not a calibrated Earth clock. Lengths beyond 64 letters are unvalidated.
- The new "experimental options" are all off by default and have only been smoke-tested and screened.

### AI disclosure

This platform was built with the AI assistant Claude (Anthropic), which wrote code and ran tests under the author's direction.
