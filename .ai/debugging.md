# Debugging

Never guess at a fix. Reproduce and understand the cause in isolation before changing package files, then apply the smallest change.

**Results first.** A suspected bug that could change numerical results (RoC values, peak points, smoothing, binning, dissimilarity, age handling) is confirmed with a test or a reproducible example and reported to the user with the evidence. Do not fix it until the user decides.

## Steps

1. **Reproduce in a throwaway script** outside the tracked package files: the session's temporary or scratch folder, or `tempdir()`. Never put debug files in `R/`, `tests/`, or the package root. The script calls `devtools::load_all()` with the package root as `path`, uses `example_data` or a minimal input, and prints enough to show what happens.
2. **Run it in a clean session** with `Rscript "<path>/debug_<topic>.R"`, never in an interactive session with leftover state. Test one hypothesis per script; iterate with `debug_<topic>2.R` if needed. On Windows, `Rscript` exits with code 1 after any warning; read the output.
3. **Probe the mechanism.** For scoping or data-masking problems, print environments, inspect internals (`body()`, `getAnywhere()`), and run the failing and working call patterns side by side.
4. **Write a test that reproduces the bug** in the function's test file and confirm it fails.
5. **Apply the minimal fix** in `R/`, with a short comment explaining why it is needed.
6. **Clean up** every debug script and output file.
7. **Run the focused test file,** then `devtools::test()`, then `devtools::check()`. A previously passing test that now fails is a new problem caused by the fix; go back to step 1.

## Reporting

State the cause, the evidence that confirmed it, the change made, and which checks in step 7 were run. If a step could not be run, say so.
