#!/bin/bash
# Copy the P3-1 files to publish into publish/ and write publish/MD5SUMS.txt.
cd "$(dirname "$0")"
P=publish
rm -rf $P; mkdir -p $P/code $P/results $P/logs $P/ref $P/literature
cp chem_p31.py chem_p31_pre_inhibition.py chem_p31_8b6da3c1.py chem_p31_pre_coldgap.py chem4c3.py cat4c3.py raf4c2.py ckpt.py thermal.py \
   validate_p31.py validate_p31_sim.py run_p31.py kill_resume_p31.sh preflight_p31.py gpu_queue.sh gpu_queue_cg.sh gpu_queue_final.sh \
   reagent_fit.py reagent_simcheck.py cold_gap.py cold_accept.py cold_accept_pre_inhibition.py inhibit_fit.py eext_uncert.py \
   make_publish.sh $P/code/
cp REPORT_P31.md cold_gap_report.md inhibit_report.md $P/
cp validate_p31.json validate_p31_t4_emptycache.json validate_p31_coldgap_partial_8b6da3c1.json validate_p31_final.json \
   reagent_fit.json reagent_simcheck.json cold_gap.json cold_accept.json inhibit_fit.json eext_uncert.json \
   preflight_p31_AS2U.json preflight_p31_AS2U_final.json $P/results/
cp regress_4c3_stored.pkl $P/results/
mkdir -p $P/results/kill_test_p31; cp -r kill_test_p31/ka kill_test_p31/kb $P/results/kill_test_p31/; cp kill_test_p31/*.log $P/logs/ 2>/dev/null
cp gpu_temp.log gpu_queue_*.log validate_p31_*run.log kill_resume_p31_run.log preflight_p31_AS2U_run.log preflight_p31_AS2U_final_run.log $P/logs/ 2>/dev/null
cp ref/part3_p31_spec.md ref/part3_alphabet_transfer.md ref/alphabet_tm.py ref/part3_cold_gap.md $P/ref/
cp ../literature/README_aitken2023.md ../literature/MD5SUMS_aitken2023.txt ../literature/aitken2023_jacs_main.txt ../literature/aitken2023_jacs_SI.txt $P/literature/
(cd $P && find . -type f ! -name MD5SUMS.txt | sort | xargs md5sum > MD5SUMS.txt)
wc -l $P/MD5SUMS.txt; du -sh $P
