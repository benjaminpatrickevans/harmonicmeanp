# The harmonic mean p-value for combining dependent tests 
## A python implementation of harmonicmeanp

This is a useful test for combining the p-values from multiple dependent tests, however, I could not find any python implementations so I wrote my own.

Included is 

- Basic hmp test (hmp.stat(..))
- Worst-case upper bound (hmp.upper_bound(..))
- Asymptotically exact hmp test (hmp.hmp(..))


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
