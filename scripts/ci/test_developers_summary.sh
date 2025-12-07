#!/usr/bin/env bash

# This script runs all tests and provides a summary of passed and failed tests.

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

passed_tests=0
failed_tests=0
failed_test_names=()

task() {
    if [[ $d =~ /CMakeFiles/ ]]; then #should not be checked
      return
    fi

    cd "$d"
    [[ $d =~ ./developers/(.*)/ ]]
    local problem_name=${BASH_REMATCH[1]}
    cmd_test="./${problem_name}_test_mastersolution"
    cmd_mastersolution="./${problem_name}_mastersolution"

    echo "-------------------------------------------------------------------"
    echo "Processing problem: $problem_name"
    echo "-------------------------------------------------------------------"

    local test_failed=0

    if [[ -f "$cmd_test" ]]; then
        echo "Executing $cmd_test"
        if eval "$cmd_test"; then
            echo -e "${GREEN}Test passed: $cmd_test${NC}"
        else
            echo -e "${RED}Test failed: $cmd_test${NC}"
            test_failed=1
        fi
    else
        echo "*** WARNING: No unit tests found in $d ***"
    fi

    # Run solution if exists:
    if [[ -f "$cmd_mastersolution" ]]; then
        echo "Executing $cmd_mastersolution"
        if output=$(eval "$cmd_mastersolution" 2>&1); then
            echo -e "${GREEN}Mastersolution executed successfully: $cmd_mastersolution${NC}"
        else
            echo -e "${RED}ERROR in mastersolution: $cmd_mastersolution${NC}"
            printf "$output\n"
            test_failed=1
        fi
    else
        echo "*** WARNING: No mastersolution found in $d ***"
    fi

    if [[ $test_failed -eq 1 ]]; then
        failed_tests=$((failed_tests + 1))
        failed_test_names+=("$problem_name")
    else
        passed_tests=$((passed_tests + 1))
    fi

    cd ../..
}

# Iterate over all subdirectories in the developers directory
for d in ./developers/*/ ; do
    # In order to not have to deal with file not found errors
    # we will first check if the directory exists
    [ -d "$d" ] || continue
    task "$d"
done

echo "==================================================================="
echo "Test Summary"
echo "==================================================================="
echo -e "${GREEN}Passed tests: $passed_tests${NC}"
echo -e "${RED}Failed tests: $failed_tests${NC}"

if [[ $failed_tests -gt 0 ]]; then
    echo "-------------------------------------------------------------------"
    echo "Failed problems:"
    for failed_test in "${failed_test_names[@]}"; do
        echo -e "${RED}- $failed_test${NC}"
    done
    echo "-------------------------------------------------------------------"
    exit 1
fi

exit 0
