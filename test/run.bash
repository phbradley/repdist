#!/bin/bash

echo
echo "# Compute repdist scores for two small human alpha chain repertoires, specified at the allele level"

cmd="../bin/repdist --outprefix test1_ --tcrs_file1 human_alpha_tcrs_file1.txt --tcrs_file2 human_alpha_tcrs_file2.txt --organism human --database ../db"

echo
echo $cmd
$cmd

echo
echo "# Comparing output to expected output (this command should produce no output):"
cmd="diff expected_test1_repdist_scores.tsv test1_repdist_scores.tsv"
echo
echo $cmd
$cmd

echo
echo "# Compute repdist scores for two small human beta chain repertoires, specified at the allele level"

cmd="../bin/repdist --outprefix test2_ --tcrs_file1 human_beta_tcrs_file1.txt --tcrs_file2 human_beta_tcrs_file2.txt --organism human --database ../db"

echo
echo $cmd
$cmd

echo
echo "# Comparing output to expected output (this command should produce no output):"
cmd="diff expected_test2_repdist_scores.tsv test2_repdist_scores.tsv"
echo
echo $cmd
$cmd

echo
echo "# Compute repdist scores for two small mouse alpha chain repertoires, specified at the allele level"

cmd="../bin/repdist --outprefix test3_ --tcrs_file1 mouse_alpha_tcrs_file1.txt --tcrs_file2 mouse_alpha_tcrs_file2.txt --organism mouse --database ../db"

echo
echo $cmd
$cmd

echo
echo "# Comparing output to expected output (this command should produce no output):"
cmd="diff expected_test3_repdist_scores.tsv test3_repdist_scores.tsv"
echo
echo $cmd
$cmd

echo
echo "# Compute repdist scores for two small mouse beta chain repertoires, specified at the allele level"

cmd="../bin/repdist --outprefix test4_ --tcrs_file1 mouse_beta_tcrs_file1.txt --tcrs_file2 mouse_beta_tcrs_file2.txt --organism mouse --database ../db"

echo
echo $cmd
$cmd

echo
echo "# Comparing output to expected output (this command should produce no output):"
cmd="diff expected_test4_repdist_scores.tsv test4_repdist_scores.tsv"
echo
echo $cmd
$cmd

