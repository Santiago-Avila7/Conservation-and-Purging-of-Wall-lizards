#!/bin/bash

# Files containing samples and chromosomes
samples_file="Samples_Used"
chromosomes_file="Chromosome_names_temp.txt"
output_task_list="task_list.txt"

# Clear the task list file if it exists
> $output_task_list

# Generate the task list
while read sample; do
  while read chr; do
    echo "${sample} ${chr}" >> $output_task_list
  done < $chromosomes_file
done < $samples_file

echo "Task list created: $output_task_list"
