# PageRank Algorithm: Theory & Efficient CSR Implementation

This repository contains a theoretical analysis and an optimized Python implementation of Google's **PageRank algorithm**. The implementation leverages the **Power Method** combined with the **Compressed Sparse Row (CSR)** matrix format to achieve high computational efficiency and a low memory footprint.

---

## Algorithm Overview

1. **Preprocessing & Dangling Nodes:** 
   The web graph is ingested into adjacency lists tracking outgoing links and backlinks. Dangling nodes (pages without outgoing links) are identified and recorded in a binary indicator vector $\mathbf{e}_{\text{dangling}}$.
   
2. **Sparse Format Encoding (CSR):** 
   The non-zero entries of the link matrix $A$ are stored using CSR arrays (`AA`, `JA`, `IA`).

3. **Power Iteration:** 
   The PageRank vector is updated iteratively without explicitly constructing the dense transition matrix $M = (1-m)A + mS$:
   $$\mathbf{v}^{(k+1)} = (1-m) \cdot \text{CSR\_matvec}(A, \mathbf{v}^{(k)}) + \text{dangling\_contrib} + \text{teleport}$$

4. **Convergence Check:** 
   Iterations terminate when $\Vert{}\mathbf{v}^{(k+1)} - \mathbf{v}^{(k)}\Vert{}_1 < \text{relTol}$.