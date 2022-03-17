# The harmonic mean p-value for combining dependent tests 
## A python implementation of harmonicmeanp

The harmonic mean p-value (hmp) is a useful test for combining the p-values from multiple dependent tests, however, I could not find any python implementations so I wrote my own.

Included is 

- Basic hmp test (*hmp.stat(..)*)
- Worst-case upper bound (*hmp.upper_bound(..)*)
- Asymptotically exact hmp test (*hmp.hmp(..)*)

**hmp.hmp** should be preferred wherever possible

Example usage

```python
import hmp

L = 50 # Number of tests

# Something that should definetly be accepted, if all tests < alpha
p_values = [0.049] *  L

print(hmp.stat(p_values)) # Basic (0.049)
print(hmp.upper_bound(p_values)) # Overly conservative worst-case (0.521)
print(hmp.hmp(p_values)) # Exact (0.0468)

# Now lets try something that should definetly be rejected, all tests >> alpha
p_values = [0.5] *  L 

print(hmp.stat(p_values)) # Basic (0.5)
print(hmp.upper_bound(p_values)) # Over Conservative (>1)
print(hmp.hmp(p_values)) # Exact (0.654)
```

If you want to pass weights in, they can be passed to any of the methods. Weights must be non-negative and sum to 1.
```python
hmp.hmp(p_values, weights=[1/L] * L) # Uniform weights
```


The original reference can be viewed [here](https://www.pnas.org/doi/10.1073/pnas.1814092116) and the r-implementation is [here](https://cran.r-project.org/web/packages/harmonicmeanp/index.html)

Any citations should go to the original author, but feel free to use and reference this implementation as you see fit.

```
@article{wilson2019harmonic,
  title={The harmonic mean p-value for combining dependent tests},
  author={Wilson, Daniel J},
  journal={Proceedings of the National Academy of Sciences},
  volume={116},
  number={4},
  pages={1195--1200},
  year={2019},
  publisher={National Acad Sciences}
}
```
