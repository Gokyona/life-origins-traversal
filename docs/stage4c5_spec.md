# Stage 4c-5 specification: coupling, L3 and L4 detection

**Status:** specification, drafted 28 October 2026, 22:40 Toronto time · **Plan:** [stage4_plan.md](stage4_plan.md), step 4c-5 · **Builds on:** [4c-4](stage4c4_spec.md) (vesicles, E3), with all of 4c-1 to 4c-4 unchanged

4c-5 adds no new chemistry or physics. The coupling between content and growth is already in 4c-4: osmotic growth (β), division, and copying inside vesicles. This step only records what vesicles do and detects two levels:

- **L3 coupled protocell:** vesicles whose contents help them grow and divide, and which pass that content on.
- **L4 Darwinian evolution:** a heritable vesicle property rises beyond what drift would give.

The dynamics must stay bit-identical to chem4c4 (validation test 1).

## 1. Vesicle lineage record

The record deferred in 4c-4 (§8 row 11) is now required. It is kept on the host (CPU memory, then disk), not on the GPU, and it is written at every checkpoint.

**Per vesicle per checkpoint:**
- vesicle id (int32), parent id, generation;
- membrane a_j and volume v_j;
- membrane gained since the previous checkpoint, as a cumulative uptake counter (osmotic and area-driven uptake, not counting the halving at division);
- divisions since the previous checkpoint;
- enclosed chain letters and enclosed letters of length ≥ L_H;
- osmotic load Π_j as used by the dynamics;
- the two dominant enclosed lineage keys, with their enclosed letters. A lineage key is as in the 4c-3 L2 detector. It is stored as a 64-bit hash, with a host-side table from hash to sequence.

**Per division event,** aggregated to checkpoint resolution: parent id, both daughter ids, and the checkpoint interval.

Vesicles lost by washout, breakdown or a missing daughter slot are marked with their cause and time.

The record draws no random numbers from the simulation stream.

## 2. Detectors

All detectors use the **joint window**. This is the longest run of consecutive checkpoints ending at the last checkpoint in which both of these hold:
- E3 strict holds (vesicles hold chains in at least half of the sites);
- the L2 lineage is present (the L2 window condition of 4c-3).

L3 and L4 are evaluated only in runs where L2_long and E3 strict both hold. A run is **evaluable** if the joint window holds at least 30 (vesicle, interval) pairs and at least 3 vesicle generations. Otherwise it is marked not evaluable, and that share is reported.

**L3 loose (content helps growth).** For every vesicle and checkpoint interval in the window, take:
- the relative growth g = (membrane gained in the interval) / (mean a_j in the interval);
- the load Π_j at the start of the interval.

Dividing by a_j removes the trivial size effect, since uptake is proportional to a_j. Compute the Spearman correlation ρ(g, Π) pooled over sites. The null is shuffled contents: 999 permutations of Π among the vesicles of the same site and checkpoint. L3 loose holds if ρ ≥ 0.1 and the permutation p-value is < 0.01.

**L3 strict (coupled protocell).** L3 loose, and content is passed on.
- Take the vesicles present at the start of the window as ancestors.
- Content is passed on if, at the last checkpoint, at least half of the sites hold a descendant at least 3 generations below its ancestor.
- That descendant must still enclose the ancestor's dominant lineage key in **chains made by copying after the window start**, with those new letters at least 0.5 times the ancestor's lineage letters at the window start. Every chain carries the step at which it was made, as a record that does not change the dynamics. Chains that existed at the window start do not count, however often they were released and re-enclosed.
- **Copied, not re-enclosed:** the share of sites passing must exceed that of the paired random-copy run (same environment, constants and seed) by at least 0.25. (Revised twice on 28 October, 23:00 and 23:30, see §6.)

**L4 (Darwinian evolution).** L3 strict, plus a heritable trait whose rise beats drift. The trait z of a vesicle is its enclosed L2 lineage letters per volume.
- **Heritability:** the parent–daughter regression of z, over the division events in the window, has a slope > 0 with p < 0.01.
- **Selection beyond drift:** for the ancestors at the start of the window, the fitness w is the number of descendants alive at the last checkpoint. The selection statistic is the Spearman correlation ρ(w, z).
  - The null is the neutral-marker control. Neutral markers assigned to the ancestors cannot change the dynamics. So a neutral marker's correlation with w has the same distribution as ρ(w, z) after permuting z among the ancestors. 999 such permutations stand in for 999 independent neutral markers, equivalent to the neutral-sequence-position markers of stage 3c.
  - L4 needs ρ(w, z) above the null's 99th percentile, with at least 10 ancestors.

**L4 sequence version (reported, not a gate).** Whether specific lineage keys rise in frequency among vesicles beyond the same permutation null. It is reported as which keys rise and by how much.

**Paired controls,** with the same environment, constants and seed:
- **β = 0:** no osmotic coupling. Over many runs, L3 loose should not pass more often than the false-positive rate.
- **Random-copy:** the copy control of 4c-3. Copies have random sequences, so lineage keys cannot persist and L3 strict should fail.

At every map point, the rates of both controls are reported next to the main rates.

**Composition fields added** (the "content of the first protocells" of the plan):
- for L3 vesicles: number of enclosed chains, length distribution, dominant lineage keys and their sequences;
- membrane size at division, generation time, and maximum generation;
- ρ(g, Π) and, for L4, ρ(w, z) as effect sizes.

## 3. Validation before any map

The criteria are fixed now, before any result.

1. **Bit-identity.** chem4c5 with the record on equals chem4c4 (b79360f6) bit for bit. This is checked on state hashes for 3 seeds at 3 points (#3709, engine12538, one point without amphiphiles).
2. **Record consistency.**
   - Generation = parent's generation + 1.
   - Daughters' membrane sums to the parent's at division, within float tolerance.
   - Enclosed letters in all vesicles plus outside letters equal the site's chain letters at each checkpoint.
   - Every vesicle id has exactly one birth and at most one end.
3. **L3 positive and negative control (fixture).** A closed site holds 8 vesicles: 4 loaded with non-copying chains (k_t = 0) and 4 empty, with amphiphile supply and shear.
   - With β = 5: ρ(g, Π) ≥ 0.3 and p < 0.01 in at least 15 of 16 seeds.
   - With β = 0: L3 loose passes in at most 1 of 16 seeds.
4. **Inheritance control (fixture).** Vesicles start with a self-copying lineage (copying on, at #3709's copying constants), with division by shear.
   - The L3 strict inheritance test passes, including the paired random-copy condition, in at least 14 of 16 seeds.
   - With k_t = 0 as the main run (paired with its own random-copy run), it passes in at most 1 of 16 seeds.
   - The raw persistence share without the paired condition is reported for the main run, the k_t = 0 run and the random-copy run. It shows how much re-enclosure alone keeps.
5. **False-positive rate at real points.** β = 0 runs at the 4 membrane-compatible 9a points and engine12538, 16 seeds each (80 runs). L3 loose and L4 pass in at most 4 of the 80 runs (5%).
6. **L4 positive control (fixture).** Two vesicle types: a self-copying lineage, and non-copying chains of the same initial mass, with β = 5 and division. L4 passes in at least 14 of 16 seeds. With β = 0, it passes in at most 1 of 16.
7. **Kill and resume** give bit-identical results, including the record.
8. **Memory pre-flight** under the standing rule:
   - worst case at batch 256: the joint points with all vesicle slots used, and the record on;
   - peak at most 80% of the 28 GB VRAM cap and of the 24 G scope;
   - no growth over at least 5 batches;
   - host RAM and disk use of the record measured and extrapolated to the map size.

**Comparison with Part I (reported, not a gate).** The plan asks whether emergent protocells behave as in Part I (stage 3a knockouts). At the joint points, report:
- L2_long with vesicles against without (amphiphile 0, the "no compartments" knockout);
- L3 with β against β = 0 (the "no coupling" knockout).

The 4c-4 validation already found the L2 margin at #3709 lower with vesicles than without (0.34 against 0.61). If enclosure hurts heredity here while compartments helped in Part I, that is a finding to report and explain, not a test failure.

**Rule.** If tests 1 to 8 all pass by these criteria, the map steps below run without a confirmation stop. If any test fails, the run stops and the failure is reported. Nothing is redesigned without review, and any redesign is recorded in an implementation table like 4c-4 §8 before any map.

## 4. Maps

The full 4c-4 map found no overlap of heredity and vesicles in 32,000 points. A full-space L3 map would therefore find nothing new. The joint nested map over the whole space belongs to 4d. 4c-5 maps where L3 and L4 can occur at all.

**5a. Point study.** 16 seeds each, with the paired β = 0 and random-copy controls, at:
- engine12538;
- the 4 membrane-compatible 9a points, with their successful membrane sets (up to 8 per point);
- every confirmed joint hit of the 9c map.

Report L3 loose, L3 strict and L4 rates per point, the effect sizes, the evaluable share, and the composition fields.

**5b. Focused map.** The box3 v2 region (4c-4 §9c), with the same masking and conditional membrane constants.
- Ground truth 500 × 16, COMPASS 16,000, batch 256.
- Levels: E1, E2, L1, L2_long, E3 strict, and joint, as before, plus L3 loose, L3 strict and L4.
- The paired controls run at every ground truth point and at every new L3 or L4 hit.
- Every new L3 strict or L4 hit is confirmed with 16 runs.

**Horizon.** Horizon, step count and time step are the same as in chem4c4. If the joint window is too short for 3 vesicle generations at most joint points, this is reported as the main limit. A longer-horizon variant would be a separate decision, not a silent change.

## 5. What this step cannot show

- **The coupling is only osmotic.** No other content-dependent growth law exists in the model. Coupling by catalysed amphiphile synthesis or by membrane-binding chains is not modelled.
- **L4 here is selection among vesicles on enclosed lineage load.** It is not the evolution of new functions: chains have no catalytic phenotype beyond the 4c-2 catalysis rules.
- **The fixtures in tests 3, 4 and 6 are constructed.** They show that the detectors work, not that coupling arises. Only 5a and 5b say where it arises.
- **The vesicle constants are layer C ranges.** Where L3 depends on β or the division constants, the result is "it depends on X", which is a pointer for experiments.

## 6. Implementation choices and revision (28 October 2026, before any GPU result)

**Second revision (23:30).** A CPU smoke test of the first revision showed that the paired comparison does not separate copying from no copying. The original 63-letter templates (16 per vesicle) stay in the vesicle family through rehydration and do not decay within the fixture horizon. So the ancestor's key persists equally in the random-copy run, and the test measured the survival of the original molecules, not copies. The test now counts only chains made by copying after the window start (§2). The paired random-copy condition stays. Under random copy, new chains do not carry the key; with k_t = 0, there are no new chains. This was done before any GPU validation result.

**Revision of L3 strict and test 4 (23:00).** On a CPU smoke test of the test 4 fixture at #3709, content persisted over 3 generations even with k_t = 0. The cause is the cycle-averaged rehydration of 4c-4 (§8 row 14). When a daughter finds no free slot, its chains are released outside, and rehydration encloses them again. So the §2 assumption of dilution by about 1/8 does not hold in wet–dry environments. The fixed 0.5 threshold was therefore replaced by a comparison with the paired random-copy run (§2). This was done before any validation result on the GPU. L3 loose (growth against load) is not affected.

The spec leaves some details open. The implementation (chem4c5.py) takes the simplest reading:

| # | Open point | Choice |
|---|---|---|
| 1 | Membrane gained | per slot: the net change of a_j from uptake, exchange and dissolution; formation (a₀) and halving at division are excluded |
| 2 | (vesicle, interval) pairs for L3 loose | slot lineage: alive at both checkpoints and never emptied or refilled by a newborn. The first daughter keeps the slot, so pairs survive divisions within an interval |
| 3 | Lineage key | whole-sequence key of chains ≥ L_H. A sequence and its complement count as one lineage (the smaller key). The top 2 per vesicle by enclosed letters |
| 4 | Trait z | letters of the 2 dominant enclosed lineages divided by the volume share φ_j |
| 5 | Heritability pairs | parent at checkpoint k−1, each daughter at checkpoint k, for divisions in interval k |
| 6 | Gate | E1, L2 loose, copy length ≥ L_H and E3 strict in the run; the runners then require L2 strict from the paired random-copy run, so the gate is L2_long ∩ E3 strict |
| 7 | Ancestry table | at most 2,048 vesicle ids per site; an overflow is flagged per run |
| 8 | End causes | washout, breakdown outside the window, dissolution below the CVC, division. A division without a free slot is recorded with no second daughter id |
| 9 | Records on disk | only for runs with an L3 or L4 hit; the pre-flight measures the full size |
| 10 | Fixtures (tests 3, 4, 6) | evaluated on the whole record from checkpoint 0, without the L2/E3 gate. Test 3: s_shear 0.02 (slow division, so pairs exist). Tests 4 and 6: #3709 environment and constants, s_shear = s_cycle = 0.01, a 6,000 h horizon, 63-letter templates, 4 starting vesicles with free slots. Test 6: 2 vesicles with copying 63-mers and 2 with non-copying chains of the same letters |
| 11 | Test 5 membrane sets | the first window-compatible 9a set per point |
| 12 | L4 sequence version | the 3 keys with the largest rise in frequency from the window start to the last checkpoint; no permutation test |

## References

As in [4c-4](stage4c4_spec.md): Chen, Roberts & Szostak 2004 for osmotic competition, and Zhu & Szostak 2009 for growth and division. Stage 3c (neutral markers) and stage 3a (knockouts) are in this repository.
