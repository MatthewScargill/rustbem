A simple but fast (ideally), multithreaded (we'll see) spectrum solver for 2D Helmholtz quantum billiard problems. Implements a Boundary Element Method (BEM) to construct boundary data, assemble boundary integral operators, and scan wavenumbers to approximate resonant modes. Includes simple spectrum and billiard visualisations. 

I first wrote this in Julia for my Master's thesis. It was very slow but it worked. As soon as there was never any reason to ever look at that code ever again (graduation), I was converted to the church of the Crab. Rewriting it in Rust has just been an excuse to tinker (but if you find yourself in a situation where this is actually useful to you then i pity you dearly)


