# Genome Read Sequencing 
# Author: Devon Shao 
# Date: 2/10/2024

# This is a function that takes in a list of reads from a donor genome as well as a reference genome, and finds all the mutations between the original donor genome and the reference genome
def find_mutations(ref_genome, read_list, kmer_length, minimizer_length, subs_mismatch_threshold):
    mutations = []  # Define a list to hold the detected mutations
    indels = [] # Define a list to hold detected indels
    #matches_errors = [0,0]
    minimizer_table = make_minimizer_table(ref_genome, kmer_length, minimizer_length) # Create, from the ref_genome, a table of minimizers and corresponding locations
    for read in read_list: # Iterate through the reads, adding detected mutations to mutations list
        #print(len(read))
        if len(read) < minimizer_length:
            continue
        read_minimizer = find_read_minimizer(read, minimizer_length) 
        minimizer_sequence = read_minimizer[0] # Define the sequence of the minimizer 
        minimizer_position_in_read = read_minimizer[1] # Define where the minimizer is within the read (start position)
        if minimizer_sequence not in minimizer_table: # Edge case where minimizer from read is not in minimizer table, some mutation happened in the minimizer
            continue
        minimizer_positions = minimizer_table[minimizer_sequence] # Define a list of all the diff occurrences of the minimizer within the ref genome
        for position in minimizer_positions:
            genome_sequence = ref_genome[position-minimizer_position_in_read:position-minimizer_position_in_read+len(read)]
            mismatch_count = hamming_distance(genome_sequence, read)
            if mismatch_count <= subs_mismatch_threshold:
                for j in range(len(read)):
                    try:
                        if genome_sequence[j] != read[j]:
                            mutations.append(">S" + str(position-minimizer_position_in_read+j+1)+ " " + genome_sequence[j] + " " + read[j] + "\n")
                    except IndexError:
                        #print(f"genome sequence:{genome_sequence}, minimizer_positions:{minimizer_positions}, position:{position}, minimizer_position_in_read:{minimizer_position_in_read}, read_length:{len(read)}")
                        pass
            else:
                j = 0
                while genome_sequence[j] == read[j] and j < len(read):
                    j += 1
                if hamming_distance(genome_sequence[j+1:], read[j:len(read)-1]) <= 1:
                    indels.append(">D" + str(position-minimizer_position_in_read+j+1) + " " + genome_sequence[j] + "\n")
                if hamming_distance(genome_sequence[j:len(read)-1], read[j+1:]) <= 1:
                    indels.append(">I" + str(position-minimizer_position_in_read+j+1) + " " + read[j] + "\n")
                    #if genome_sequence[j] != read[j]:
                        #mutations.append(">S" + str(position-minimizer_position_in_read+j+1)+ " " + genome_sequence[j] + " " + read[j] + "\n")
                
    mutations = consensus(mutations, 2)
    indels = consensus(indels, 2)
    mutations.extend(indels)
    return(mutations)

# This function takes in a reference genome and specificed kmer and minimizer lenghts and makes a minimizer hashmap with keys as minimizers and values as the positions of the minimizers in the genome
def make_minimizer_table(ref_genome, kmer_length, minimizer_length): 
    minimizer_table = {}  # Define output table, should be a dictionary with keys = string of minimizer sequence and values = all positions of those sequences in ref_genome
    for outer_kmer_start_pos in range(len(ref_genome)-kmer_length + 1):  # Iterate through all k-mers of ref_genome 
        inner_kmers = {}
        current_outer_kmer = ref_genome[outer_kmer_start_pos:outer_kmer_start_pos+kmer_length]
        for inner_kmer_start_pos in range(kmer_length - minimizer_length + 1):    # Iterate through all k-mers of the current k-mer, adding each to a list and then sorting the list to take the smallest value as the minimizer 
            inner_kmers[current_outer_kmer[inner_kmer_start_pos:inner_kmer_start_pos + minimizer_length]] = inner_kmer_start_pos + outer_kmer_start_pos
        sorted_inner_kmers = sorted(inner_kmers.keys())
        minimizer = sorted_inner_kmers[0]
        if minimizer in minimizer_table:
            if inner_kmers[minimizer] not in minimizer_table[minimizer]:
                minimizer_table[minimizer].append(inner_kmers[minimizer])
        else:
            minimizer_table[minimizer] = [inner_kmers[minimizer]]
    return(minimizer_table)


# This function takes a read requence and specified minimizer length and returns the minimizer of the read sequence
def find_read_minimizer(read, minimizer_length):
    minimizer = []  # return a two-item list consisting of a string, the minimizer, and an int, the position of the first char of the kmer
    kmers = {}  # The kmers from which the lexicographically smallest one will become the minizer for the currrent read
    for kmer_start_pos in range(len(read) - minimizer_length + 1):
        kmers[read[kmer_start_pos:kmer_start_pos + minimizer_length]] = kmer_start_pos
    sorted_kmers = sorted(kmers.keys())
    minimizer.append(sorted_kmers[0])
    minimizer.append(kmers[sorted_kmers[0]])
    return minimizer

# This is a function that takes in two strings and returns their hamming distance
def hamming_distance(p:str, q:str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(p, q))

# This is a function that is ran to eliminate read errors. It takes all the detected mutations and only keeps the mutations that appear more than the specified threshold amount of times. 
def consensus(mutations, min_occurrences_threshold):
    counts = {}

    for mutation in mutations:
        if mutation in counts:
            counts[mutation] += 1
        else:
            counts[mutation] = 1

    filtered_mutations = [mutation for mutation, count in counts.items() if count >= min_occurrences_threshold]

    return filtered_mutations
