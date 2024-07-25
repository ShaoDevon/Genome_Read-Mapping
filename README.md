# Read-Mapping

In human genome sequencing, it is often too difficult or expensive to sequence a whole genome. To solve this problem, next generation sequencing is implemented, where multiple short reads of the genome are sequenced and compared to an already-known reference genome. A problem arises: how can we recover the original sequence from the multiple short reads? Furthermore, the human genome is very long, and sequencers generate a large amount of reads. Solving the problem of recovering the original sequence is not only about completing it, but completing it in a space-efficient and time-efficient manner. 

This project aims to solve the problem in this manner. It provides an implementation of the minimizers method: from the genome a dictionary where the keys are minimizers and the values are their locations in the genome is generated, and using the minimizer table each read is mapped to its optimal location in the reference genome, thus reconstructing the donor genome. 
