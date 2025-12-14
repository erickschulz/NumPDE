#!/usr/bin/env bash
# Non-blocking test runner with summary reporting
# Tests continue even on failures to provide complete picture

# Arrays to track results
declare -a PASSED_TESTS
declare -a FAILED_TESTS
declare -a FAILED_MAINS
declare -a NO_TESTS

task() {
    if [[ $d =~ /CMakeFiles/ ]]; then
      return
    fi

    cd "$d"
    [[ $d =~ ./developers/(.*)/ ]]
    PROJECT_NAME="${BASH_REMATCH[1]}"
    cmd_test="./${PROJECT_NAME}_test_mastersolution"
    cmd_mastersolution="./${PROJECT_NAME}_mastersolution"

    # Run test executable if it exists
    if [[ -f "$cmd_test" ]]; then
      echo "Executing $cmd_test"
      if eval "$cmd_test --gtest_color=yes" ; then
        PASSED_TESTS+=("$PROJECT_NAME")
      else
        FAILED_TESTS+=("$PROJECT_NAME")
      fi
    else
      echo "*** WARNING: No unit tests found in $d ***"
      NO_TESTS+=("$PROJECT_NAME")
    fi

    # Run main executable if it exists
    if [[ -f "$cmd_mastersolution" ]]; then
      echo "Executing $cmd_mastersolution"
      if output=$(eval $cmd_mastersolution 2>&1) ; then
        : # Success, do nothing
      else
        echo "ERROR: $cmd_mastersolution failed"
        printf "$output"
        FAILED_MAINS+=("$PROJECT_NAME")
      fi
    fi

    cd ../..
}

# Process all developer projects
for d in ./developers/*/ ; do
  task $d
done

# Print summary
echo ""
echo "========================================="
echo "         TEST SUMMARY"
echo "========================================="
echo ""
echo "PASSED: ${#PASSED_TESTS[@]} projects"
for project in "${PASSED_TESTS[@]}"; do
  echo "  ✓ $project"
done

if [[ ${#FAILED_TESTS[@]} -gt 0 ]]; then
  echo ""
  echo "FAILED: ${#FAILED_TESTS[@]} projects"
  for project in "${FAILED_TESTS[@]}"; do
    echo "  ✗ $project"
  done
fi

if [[ ${#FAILED_MAINS[@]} -gt 0 ]]; then
  echo ""
  echo "MAIN EXECUTABLE FAILURES: ${#FAILED_MAINS[@]} projects"
  for project in "${FAILED_MAINS[@]}"; do
    echo "  ✗ $project"
  done
fi

if [[ ${#NO_TESTS[@]} -gt 0 ]]; then
  echo ""
  echo "NO TESTS: ${#NO_TESTS[@]} projects"
  for project in "${NO_TESTS[@]}"; do
    echo "  ⚠ $project"
  done
fi

echo ""
echo "========================================="
echo "Total: $((${#PASSED_TESTS[@]} + ${#FAILED_TESTS[@]})) tests executed"
echo "========================================="

# Exit with failure if any tests failed (GitHub Actions will show RED)
# But CI is configured to continue-on-error so pipeline won't stop
if [[ ${#FAILED_TESTS[@]} -gt 0 ]] || [[ ${#FAILED_MAINS[@]} -gt 0 ]]; then
  exit 1
else
  exit 0
fi
