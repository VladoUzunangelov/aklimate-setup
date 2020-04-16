---
# Thu Apr 16 15:18:22 PDT 2020 chrisw


I tried to build a new docker environment from scratch. I encountered some `index out of bounds` errors with it. When I switched back to an older version to use as base container, there were no errors.

stuartlab/aklimate:dev	d65e1e4376c3	Old environment, no SPICER, no AKLIMATE, works
stuartlab/aklimate:dev_spicer	b7e3d2e0214e	Newer environment built on top of d65e1e4376c3, includes SPICER, no AKLIMATE, works
stuartlab/aklimate:dev_no_spicer	8533e26d7b3b	Newer environment built on top of ubuntu:18.04 (newer version of stuartlab/aklimate:dev), no SPICER, no AKLIMATE, out-of-bounds error

The R sessionInfo of each environment is saved.

