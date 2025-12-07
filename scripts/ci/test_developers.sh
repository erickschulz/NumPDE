#!/usr/bin/env bash
# This is a slightly modified version of compile_test.sh to fit gitlab ci workflow

task() {
  
    if [[ $d =~ /CMakeFiles/ ]]; then #should not be checked
      return
    fi

        cd "$d"
        [[ $d =~ ./developers/(.*)/ ]]
        cmd_test="./${BASH_REMATCH[1]}_test_mastersolution";
        cmd_mastersolution="./${BASH_REMATCH[1]}_mastersolution";

        if [[ -f "$cmd_test" ]]; then
          echo "Executing $cmd_test";
          eval "$cmd_test --gtest_color=yes"
        else
          echo "*** WARNING: No unit tests found in $d ***";
        fi
        # Run solution if exists:
        if [[ -f "$cmd_mastersolution" ]]; then
          echo "Executing $cmd_mastersolution";
          if ! output=$(eval $cmd_mastersolution 2>&1) ; then
            echo "ERROR: ";
            printf "$output";
            exit 1;
          fi
        else
          echo "*** WARNING: No mastersolution found in $d ***";
        fi
        cd ../..
}
set -e
for d in ./developers/*/ ;
do
  task $d #parallelizing did not work, can be modified at a later stage
done
