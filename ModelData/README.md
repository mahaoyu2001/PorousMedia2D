# Pore-centre data

Each text file contains two whitespace-separated columns, `x` and `y`, with
one circular-pore centre per row. The supplied RVE is `1 x 1`, and the pore
radius is `a = 0.005`. Some centres lie just outside the RVE so that pores
crossing an external boundary are clipped by the Boolean-cut operation; this
is intentional and those rows must not be discarded.
