Python code to-do list
----------------------

- If the spectral range of interest is relatively narrow, so that it does not include one of the strong absorption peaks of the gas, then the cross-section within that range will be biased because it is missing the long-tail of the long peak that pushes up the minimum absorption level. For example, if you try to plot SO2 in the range 1.4-1.6um, it will say that the cross-section is zero. But if you try to plot it in the wide range of, say, 1.0-18um, then you will see that the 1.4-1.6um range has only the tail part and no independent absorption peaks inside it. Figure out a way to incorporate the tails that lie outside the spectral range of interest.
- Incorporate the "total partition sums" into the digitzation of the spectral lines.
